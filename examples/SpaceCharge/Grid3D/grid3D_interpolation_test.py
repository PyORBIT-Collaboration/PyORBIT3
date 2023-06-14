#-----------------------------------------------------
#Grid3D interpolation and gradient test
#-----------------------------------------------------

import sys
import math
import orbit.core

from spacecharge import Grid3D


print("Start.")

sizeX = 50
sizeY = 50
sizeZ = 50
xMin = -2.0
xMax = +2.0
yMin = -4.0
yMax = +4.0
zMin = -2.0
zMax = +2.0

grid3D = Grid3D(sizeX,sizeY,sizeZ)
grid3D.setGridX(xMin,xMax)
grid3D.setGridY(yMin,yMax)
grid3D.setGridZ(zMin,zMax)

grid3D.longWrapping(True)

def Func(x,y,z):
	return 1.0/(x*x+y*y+z*z)

def FuncGrad(x,y,z):
	return (-2*x/(x*x+y*y+z*z)**2,-2*y/(x*x+y*y+z*z)**2,-2*z/(x*x+y*y+z*z)**2)

for ix in range(sizeX):
	x = grid3D.getGridX(ix)
	for iy in range(sizeY):
		y = grid3D.getGridY(iy)
		for iz in range(sizeZ):
			z = grid3D.getGridZ(iz)
			val = Func(x,y,z)
			grid3D.setValue(val,ix,iy,iz)
	
#----------------------------------------------
# let's create a straight line along gradient
#----------------------------------------------

nPoints = 10000
gradX = 0.0
gradY = 0.0
gradZ = 0.5

x0 = 1.
y0 = 1.
z0 = 0.

t_min = -20.
t_max = +20.
step_t = (t_max - t_min)/nPoints

traj_arr = []
t_arr = []

for ip in range(nPoints):
	t = t_min + step_t*ip
	x = x0 + gradX*t
	y = y0 + gradY*t
	z = z0 + gradZ*t
	if(x > xMin and x < xMax and y > yMin and y < yMax and z > zMin and z < zMax):
		traj_arr.append([x,y,z])
		t_arr.append(t)

#------- trajectory is ready

grid_val_arr = []
func_val_arr = []
diff_val_arr = []

grid_gradX_val_arr = []
func_gradX_val_arr = []
diff_gradX_val_arr = []

grid_gradY_val_arr = []
func_gradY_val_arr = []
diff_gradY_val_arr = []

grid_gradZ_val_arr = []
func_gradZ_val_arr = []
diff_gradZ_val_arr = []

for t_ind in range(len(t_arr)):
	t = t_arr[t_ind]
	[x,y,z] = traj_arr[t_ind]
	grid_val = grid3D.getValue(x,y,z)
	func_val = Func(x,y,z)
	grid_val_arr.append(grid_val)
	func_val_arr.append(func_val)
	diff_val_arr.append(abs(func_val - grid_val))
	#-------------------------------
	(func_gradX,func_gradY,func_gradZ) = FuncGrad(x,y,z)
	(grid_gradX,grid_gradY,grid_gradZ) = grid3D.calcGradient(x,y,z)
	#---------
	func_gradX_val_arr.append(func_gradX)
	grid_gradX_val_arr.append(grid_gradX)
	diff_gradX_val_arr.append(abs(func_gradX - grid_gradX))
	#---------
	func_gradY_val_arr.append(func_gradY)
	grid_gradY_val_arr.append(grid_gradY)
	diff_gradY_val_arr.append(abs(func_gradY - grid_gradY))	
	#---------
	func_gradZ_val_arr.append(func_gradZ)
	grid_gradZ_val_arr.append(grid_gradZ)
	diff_gradZ_val_arr.append(abs(func_gradZ - grid_gradZ))	
	
import matplotlib.pyplot as plt

plt.plot(t_arr,grid_val_arr)
plt.plot(t_arr,func_val_arr)
plt.ylabel('Values, arb. units')
plt.xlabel('t')
#plt.axis([-1., 1., 0, 0.001])

plt.savefig('interpolation_in_grid3d.png')
plt.show()

plt.plot(t_arr,diff_val_arr)
plt.ylabel('Diff. in Values, arb. units')
plt.xlabel('t')
#plt.axis([-1., 1., 0, 0.001])

plt.savefig('interpolation_in_grid3d_diff.png')
plt.show()

plt.plot(t_arr,func_gradX_val_arr)
plt.plot(t_arr,grid_gradX_val_arr)
plt.ylabel('Gradient X Values, arb. units')
plt.xlabel('t')
#plt.axis([-1., 1., 0, 0.001])

plt.savefig('interpolation_in_grid3d_gradientX.png')
plt.show()

plt.plot(t_arr,func_gradY_val_arr)
plt.plot(t_arr,grid_gradY_val_arr)
plt.ylabel('Gradient Y Values, arb. units')
plt.xlabel('t')
#plt.axis([-1., 1., 0, 0.001])

plt.savefig('interpolation_in_grid3d_gradientY.png')
plt.show()

plt.plot(t_arr,func_gradZ_val_arr)
plt.plot(t_arr,grid_gradZ_val_arr)
plt.ylabel('Gradient Z Values, arb. units')
plt.xlabel('t')
#plt.axis([-1., 1., 0, 0.001])

plt.savefig('interpolation_in_grid3d_gradientZ.png')
plt.show()

plt.plot(t_arr,diff_gradZ_val_arr)
plt.ylabel('Diff. Gradient Z Values, arb. units')
plt.xlabel('t')
#plt.axis([-1., 1., 0, 0.001])

plt.savefig('interpolation_in_grid3d_diff_gradientZ.png')
plt.show()

sys.exit(0)
