import math
import numpy as np
import cv2

def show(img):
	cv2.imshow("image", img)
	cv2.waitKey(0)
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


class camera:
	def __init__(self, theta, phi, fov, width, height, x, y, z):
		rays = np.zeros((h, w, 5), np.float32)
		Y, X = np.mgrid[:height, :width]
		pAdj = (Y - ((Y.shape[0] - 1) / 2))
		tAdj = (X - ((X.shape[1] - 1) / 2))
		maxAdj = max(np.max(np.absolute(pAdj)), np.max(np.absolute(tAdj)))
		pAdj /= maxAdj
		pAdj *= (fov / 2)
		tAdj /= maxAdj
		tAdj *= (fov / 2)
		rays[:,:,0] = x
		rays[:,:,1] = y
		rays[:,:,2] = z
		rays[:,:,3], rays[:,:,4] = tAdj + theta, pAdj + phi
		print(f"{np.max(rays[:,:,3])}, {rays[rays.shape[0] // 2, rays.shape[1] // 2, 3]}, {np.min(rays[:,:,3])}")
		print(f"{np.max(rays[:,:,4])}, {rays[rays.shape[0] // 2, rays.shape[1] // 2, 4]}, {np.min(rays[:,:,4])}")
		
		self.rays = rays
		
		
		


class Op:
	def pM(a):
		return np.array([a, -a], np.float64)
	def mP(a):
		return np.array([-a, a], np.float64)

def findCircleIntersection(Circle, Ray):
	ly = (-Ray[0] / Ray[1])
	lyt = (-Ray[2] / Ray[1])
	a = 1 + ly ** 2
	b = -2 * Circle[0] + 2 * (lyt - Circle[1]) * ly
	c = Circle[0] ** 2 + (lyt - Circle[1]) ** 2 + (-Circle[2] * Circle[2])
	if(b**2 - 4 * a * c < 0):
		return False
	X = (-b + Op.pM(math.sqrt(b**2 - (4 * a * c)))) / (2 * a)
	Y = X * Ray[0] + Ray[2]
	inters = [[X[d], Y[d]] for d in range(0, X.shape[0])]

	return inters
	
precisionLevel = 0xFFFFFFFF
def truncateFloatPointError(a):
	global precisionLevel
	return np.int64(precisionLevel * a) / precisionLevel


def findSphereIntersection(Sphere, Ray):
	global precisionLevel
	r = Sphere[3]
	Sphere = np.array(Sphere, np.float64)[:3]
	vx = truncateFloatPointError(math.cos(Ray[3] / 180 * math.pi) * math.sin(Ray[4] / 180 * math.pi))
	vy = truncateFloatPointError(math.sin(Ray[3] / 180 * math.pi) * math.sin(Ray[4] / 180 * math.pi))
	vz = truncateFloatPointError(math.cos(Ray[4] / 180 * math.pi))
	# ~ print(Ray[3], Ray[4])
	vect = np.array([vx, vy, vz], np.float32)
	vect /= np.max(vect)
	
	maxindex = np.argmax(vect)
	
	lxt = (Ray[0] - Ray[maxindex]) / vect[maxindex]
	lyt = (Ray[1] - Ray[maxindex]) / vect[maxindex]
	lzt = (Ray[2] - Ray[maxindex]) / vect[maxindex]
	cord = np.array([lxt, lyt, lzt], np.float32)
	# ~ print(vect)
	a = np.sum(vect ** 2)
	b = np.sum(2 * (cord - Sphere) * vect)
	c = np.sum((cord - Sphere) ** 2) - r ** 2
	
	if(b**2 - 4 * a * c < 0):
		return [0, 0, 0]
	# ~ print("yippee!")
	ans = (-b + Op.pM(math.sqrt(b**2 - (4 * a * c)))) / (2 * a)
	X = ans * vect[0]
	Y = ans * vect[1]
	Z = ans * vect[2]
	inters = np.array([[X[d], Y[d], Z[d]] for d in range(0, X.shape[0])], np.float32)
	return inters[0] if np.sum(np.power(inters[0] - Ray[:3], 2)) ** 0.5 < np.sum(np.power(inters[1] - Ray[:3], 2)) ** 0.5 else inters[1]

def calculateNormal(Sphere, Intersection, Ray):
	if((Intersection[0] == 0 and Intersection[1] == 0 and Intersection[2] == 0)):
		return spherical2Cartesian(Ray[3], Ray[4], 0)
		# ~ return (Ray[3], Ray[4])
	# ~ print("aye carumba!")
	theta, phi, _ = cartesian2Spherical((Intersection[0] - Sphere[0]), (Intersection[1] - Sphere[1]), (Intersection[2] - Sphere[2]))
	return spherical2Cartesian(theta, phi, 1)
	# ~ return theta, phi
	
def Render(Cam, Sphere):
	intersections = np.array([[findSphereIntersection(Sphere, Cam.rays[a, b]) for b in range(0, Cam.rays.shape[1])] for a in range(0, Cam.rays.shape[0])], np.float32)
	print("intersected!")
	Normals = np.array([[calculateNormal(Sphere, intersections[a, b], Cam.rays[a, b]) for b in range(0, Cam.rays.shape[1])] for a in range(0, Cam.rays.shape[0])], np.float32)
	show(normalize(Normals))
	print(Normals.shape)
	print(Cam.rays.shape)
	print("normalized!")
	x, y, z = spherical2Cartesian(Cam.rays[:,:,3], Cam.rays[:, :, 4], 1)
	Rays = np.array([[[x[a, b], y[a, b], z[a, b]] for b in range(0, x.shape[1])] for a in range(0, x.shape[0])], np.float32)
	Frame = np.zeros((Cam.rays.shape[0], Cam.rays.shape[1]), np.float32)
	Frame[:,:] = truncateFloatPointError(Normals[:, :, 0] - Rays[:, :, 0] + Normals[:, :, 1] - Rays[:, :, 1] + Normals[:, :, 2] + Rays[:, :, 2])
	# ~ Frame[Normals[:,:,:] == [0, 0, 0]] = [0] 
	print(np.max(Frame))
	print(np.min(Frame))
	# ~ Frame[not intersections[:,:].any()] = np.min(Frame)
	
	return normalize(Frame)
	
	
	






XT = 20
YT = 10
ZT = -12
Sphere = np.array([20, 10, 1, 10], np.float32)
theta, phi, _ = cartesian2Spherical(Sphere[0] - XT, Sphere[1] - YT, Sphere[2] - ZT)
print(theta, phi)
fov = 90
fov = fov / 180 * math.pi
w = 500
h = 500
cam = camera(theta, phi, fov, w, h, XT, YT, ZT)
# ~ light = [0, 0]
show(Render(cam, Sphere))
