import math
import numpy as np
import cv2
import time

distancePower = 1
world = []
class camera:
	def facePoint(self, coords):
		theta, phi, _ = cartesian2Spherical(coords[0] - self.origin[0], coords[1] - self.origin[1], coords[2] - self.origin[2])
		self.rotate(phi - self.phi, theta - self.theta)

	def faceDir(self, theta, phi):
		self.rotate(theta - self.theta, phi - self.phi)

	def rotate(self, theta, phi):
		self.theta += theta
		self.phi += phi
		R1 = np.array([[math.cos(self.theta), -math.sin(self.theta),  0], [math.sin(self.theta), math.cos(self.theta), 0], [0, 0, 1]], np.float32)
		R2 = np.array([[math.cos(self.phi), 0, math.sin(self.phi)], [0, 1, 0], [-math.sin(self.phi), 0, math.cos(self.phi)]], np.float32)
		Matrix = R1@R2
		self.coords[:,:,0] = self.screen[:, :, 0] * Matrix[0, 0] + self.screen[:, :, 1] * Matrix[0, 1] + self.screen[:, :, 2] * Matrix[0, 2]
		self.coords[:,:,1] = self.screen[:, :, 0] * Matrix[1, 0] + self.screen[:, :, 1] * Matrix[1, 1] + self.screen[:, :, 2] * Matrix[1, 2]
		self.coords[:,:,2] = self.screen[:, :, 0] * Matrix[2, 0] + self.screen[:, :, 1] * Matrix[2, 1] + self.screen[:, :, 2] * Matrix[2, 2]
		self.coords = parallelNormalizeVector(self.coords)
		
	def __init__(self, theta, phi, fov, resolution, x, y, z):
		self.origin = np.array([x, y, z], np.float32)
		Y, X = np.mgrid[:resolution[0], :resolution[1]]
		coords = np.zeros((resolution[0], resolution[1], 3), np.float32)
		width, height = resolution[0], resolution[1]
		X = (X - ((X.shape[1] - 1) / 2))
		Y = (Y - ((Y.shape[0] - 1) / 2))
		coords[:,:,0] = X[:,:]
		coords[:,:,1] = Y[:,:]
		coords[:,:,2] = np.array(spherical2Cartesian(fov / 2, fov / 2, np.array(cartesian2Spherical(X[0, -1], Y[-1, 0], math.sqrt(((width / 2) ** 2 + (height / 2) ** 2))), np.float32)[2]), np.float32)[2]
		self.screen = coords * 1
		self.coords = coords * 1
		self.theta = 0
		self.phi = 0
		self.normals = coords * 0
		self.frame = np.zeros((resolution[0], resolution[1]), np.float32)
		self.faceDir(theta, phi)
		
	def FindSphereIntersections(self, Sphere):
		r = Sphere[3]
		
		Sphere = np.array(Sphere, np.float32)[:3] - self.origin
		# ~ print(Sphere)
		vect = (self.coords) * 1
		maggs = np.amax(np.absolute(vect), 2)
		
		vect[:,:,0] /= maggs
		vect[:,:,1] /= maggs
		vect[:,:,2] /= maggs
		
		output = vect*0
		
		a = np.sum(vect ** 2, axis=-1)
		b = np.sum(2 * (-Sphere) * vect, axis = -1)
		c = np.sum((vect * 0 - Sphere) ** 2, axis = -1) - r ** 2
		
		truths = b**2 - 4 * a * c < 0
		
		a[truths] = 1
		b[truths] = 0
		c[truths] = 0
		
		ans = (-b - np.sqrt(b**2 - (4 * a * c))) / (2 * a)
		ans2 = (-b + np.sqrt(b**2 - (4 * a * c))) / (2 * a)
		ans[ans < 0] = ans2[ans < 0]
		ans[ans < 0] = 0
		
		inters = self.coords * 0
		
		inters[:,:,0] = ans * vect[:, :, 0] #+ self.origin[0]
		inters[:,:,1] = ans * vect[:, :, 1] #+ self.origin[1]
		inters[:,:,2] = ans * vect[:, :, 2] #+ self.origin[2]
		# ~ start = time.time()
		conda = np.logical_xor(self.coords < 0, inters > 0)
		cond = 1 - np.logical_and(conda[:,:,0], np.logical_and(conda[:,:,1], conda[:,:,2]))
		truths = np.logical_or(cond, truths)
		return truths, inters
		
	def RenderSpheres(self, objects):
		norm = self.normals * 0
		intersections = norm * 1
		Distances = norm[:, :, 0] * 1
		AllTruths = norm[:, :, 0] + 1
		for a in objects:
			# ~ print(a.shape)
			truths, newInters = self.FindSphereIntersections(a)
			LocalNorms = newInters - np.array(a, np.float32)[:3]
			AllTruths = np.logical_and(truths, AllTruths)
			dists = np.sqrt(np.sum(newInters**2, axis=-1))
			cond = np.logical_and(np.logical_or(dists < Distances, Distances == 0), dists != 0)
			Distances[cond] = dists[cond]
			norm[cond] = LocalNorms[cond]
		bright = np.sum(parallelNormalizeVector(norm) * self.coords, axis = -1) - Distances
		bright -= np.min(bright)
		norm -= np.min(norm)
		bright[AllTruths] = 0
		norm[AllTruths] = self.coords[AllTruths] - 100
		frame = norm * 1
		frame[:,:,0] *= bright
		frame[:,:,1] *= bright
		frame[:,:,2] *= bright
		self.frame = normalize(frame)
		show(self.frame, 1)
			

elasticity = 0.99




GravityConstant = 6.67408e-11
DTime = 1 / 5000
class HeavenlyBody:
	def __init__(self, location, radius, density = 1, velocityVect = [0, 0, 0]):
		self.mass = 4/3 * math.pi * radius**3 * density
		self.radius = radius
		self.location = location
		self.velocity = velocityVect
		self.FrameForce = self.velocity*0
		
	def updateFrameForce(self, others):
		Force = self.velocity * 0
		for a in others:
			massRatio = self.mass / (self.mass + a.mass)
			r = ((self.location[0] - a.location[0]) ** 2 + (self.location[1] - a.location[1]) ** 2 + (self.location[2] - a.location[2]) ** 2)**0.5
			F = (a.mass * self.mass * GravityConstant) / r**2
			T, P, _ = cartesian2Spherical(a.location[0] - self.location[0], a.location[1] - self.location[1], a.location[2] - self.location[2])
			Xm, Ym, Zm = spherical2Cartesian(T, P, F)
			if(r < self.radius + a.radius):
				n = normalizeVector(self.location - a.location)
				dot = np.sum(self.velocity * n)
				if(dot < 0):
					Force -= self.velocity
					newDirection = normalizeVector(self.velocity - 2 * dot * n)
					Pself = self.velocity * self.mass
					Pother = a.velocity * a.mass
					TotalMagnitude = getVectorMagnitude(Pself) + getVectorMagnitude(Pother)
					Force += newDirection * (TotalMagnitude * massRatio) / self.mass * elasticity
			else:
				Force += (np.array([Xm, Ym, Zm], np.float64) * massRatio * a.mass / self.mass) * DTime
		self.FrameForce = Force

	def updatePosition(self):
		self.velocity += self.FrameForce
		self.location += self.velocity * DTime

	def asSphere(self):
		return np.array([self.location[0], self.location[1], self.location[2], self.radius], np.float32)







def show(img, time = 0):
	cv2.imshow("image", img)
	cv2.waitKey(time)
	if time == 0:
		cv2.destroyAllWindows()

def normalize(img):
	img = img * 1.0
	img -= np.min(img)
	img /= np.max(img)
	img *= 255.999
	return np.uint8(img)

def cartesian2Spherical(x, y, z):
	return np.arctan2(y,x), np.arctan2(np.sqrt(np.power(x, 2) + np.power(y, 2)), z),  np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
	
def spherical2Cartesian(theta, phi, rho = 1):
	return rho * np.cos(theta) * np.sin(phi), rho * np.sin(theta) * np.sin(phi), rho * np.cos(phi)
	
def getVectorMagnitude(Vect):
	return np.sum(np.power(Vect, 2)) ** 0.5
	
def normalizeVector(Vect):
	t = np.arctan2(Vect[1], Vect[0])
	p = np.arctan2(np.sqrt(np.power(Vect[0], 2) + np.power(Vect[1], 2)), Vect[2])
	return np.array([np.cos(t) * np.sin(p), np.sin(t) * np.sin(p), np.cos(p)], np.float32)

def parallelNormalizeVector(Vect):
	t = np.arctan2(Vect[:,:,1], Vect[:,:,0])
	p = np.arctan2(np.sqrt(np.power(Vect[:,:,0], 2) + np.power(Vect[:,:,1], 2)), Vect[:,:,2])
	Vect[:,:,0] = np.cos(t) * np.sin(p)
	Vect[:,:,1] = np.sin(t) * np.sin(p)
	Vect[:,:,2] = np.cos(p)
	return Vect


XT = 0
YT = 0
ZT = 0

neptune = HeavenlyBody(np.array([-20, -10, 40], np.float32), 12, 5, np.array([400, 0, 0], np.float32))
saturn = HeavenlyBody(np.array([20, 0, 40], np.float32), 10, 5, np.array([-400, 0, 0], np.float32))
# ~ uranus = HeavenlyBody(np.array([0, 0, 100], np.float32), 2, 1, np.array([0, 0, 0], np.float32))


# ~ theta, phi, _ = cartesian2Spherical(neptune.location[0] - XT, neptune.location[1] - YT, neptune.location[2] - ZT)
theta = 0
phi = 0
fov = 120
fov = fov / 180 * math.pi
w = 600
h = 600
cam = camera(theta, phi, fov, [w,h], XT, YT, ZT)
x = 0
import time
# ~ frames = []
while True:
	start = time.time()
	cam.RenderSpheres([neptune.asSphere(), saturn.asSphere()])
	# ~ cam.facePoint(saturn.location)
	neptune.updateFrameForce([saturn])
	saturn.updateFrameForce([neptune])
	# ~ uranus.updateVelocity([saturn, neptune])
	neptune.updatePosition()
	saturn.updatePosition()
	# ~ print(getVectorMagnitude(neptune.velocity * neptune.mass) + getVectorMagnitude(saturn.velocity * saturn.mass), end = "\r")
	end = time.time()
	print(f"{int(1 / (end - start))}fps    ", end = "\r")
	# ~ frames.append(cam.frame)

# ~ for x in range(400):
	# ~ cv2.imwrite(f'pictors/{100+x}.png', frames[x])


