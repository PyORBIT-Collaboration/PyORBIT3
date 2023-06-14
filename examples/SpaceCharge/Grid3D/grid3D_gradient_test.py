#-----------------------------------------------------
#Grid3D gradient test
#-----------------------------------------------------

import sys
import math
import orbit.core

from spacecharge import Grid3D

print("Start.")

sizeX = 50
sizeY = 50
sizeZ = 50
xMin =  2.0
xMax = +4.0
yMin =  1.0
yMax = +4.0
zMin =  1.0
zMax = +4.0

grid3D = Grid3D(sizeX,sizeY,sizeZ)
grid3D.setGridX(xMin,xMax)
grid3D.setGridY(yMin,yMax)
grid3D.setGridZ(zMin,zMax)

def Func(x,y,z):
	return 1.0/(x*x+y*y+z*z)
	#return (x**2+y**2+z**2)

def FuncGrad(x,y,z):
	return (-2*x/(x*x+y*y+z*z)**2,-2*y/(x*x+y*y+z*z)**2,-2*z/(x*x+y*y+z*z)**2)
	#return (2*x,2*y,2*z)


for ix in range(sizeX):
	x = grid3D.getGridX(ix)
	for iy in range(sizeY):
		y = grid3D.getGridY(iy)
		for iz in range(sizeZ):
			z = grid3D.getGridZ(iz)
			val = Func(x,y,z)
			grid3D.setValue(val,ix,iy,iz)

diff_max = 0.
max_point = (0.,0.)
for ix in range(sizeX):
	x = grid3D.getGridX(ix) + 0.5*(grid3D.getMaxX() - grid3D.getMinX())/grid3D.getSizeX()
	for iy in range(sizeY):
		y = grid3D.getGridY(iy) + 0.5*(grid3D.getMaxY() - grid3D.getMinY())/grid3D.getSizeY()
		for iz in range(sizeZ):
			z = grid3D.getGridZ(iz) + 0.5*(grid3D.getMaxZ() - grid3D.getMinZ())/grid3D.getSizeZ()
			#print "(x,y,z)=",(x,y,z)," val=",grid3D.getValue(x,y,z)," grad=",grid3D.calcGradient(x,y,z)
			(gradX0,gradY0,gradZ0) = FuncGrad(x,y,z)
			(gradX,gradY,gradZ) = grid3D.calcGradient(x,y,z)
			grad_val = math.sqrt((gradX)**2 + (gradY)**2 + (gradZ)**2)
			diff = 100.*math.sqrt((gradX - gradX0)**2 + (gradY - gradY0)**2 + (gradZ - gradZ0)**2)/grad_val
			if(diff > diff_max):
				diff_max = diff
				max_point = (x,y,z)
		
print("max diff % =",diff_max)
print("max point =",max_point)
(x,y,z) = max_point
print("Func at max point         =",Func(x,y,z))
print("Grid3D value at max point =",grid3D.getValue(x,y,z))
print("Grad at max point =",FuncGrad(x,y,z))
print("Grid3D Grad at max point =",grid3D.calcGradient(x,y,z))

print("Stop.")

