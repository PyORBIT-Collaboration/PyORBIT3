#-----------------------------------------------------
#Grid1D gradient test
# There are grid1D.calcGradient(z) or grid1D.calcGradientSmoothed(z) functions.
# The smoothed version should be more precise.
#-----------------------------------------------------

import sys
import math

from orbit.space_charge import Grid1D


print("Start.")

sizeZ = 50
sizeZp = 51
zMin = 1.0
zMax = 2.0
sizeTest = 100
sizeTestp = 101
dz = (zMax - zMin) / sizeTest

grid1D = Grid1D(sizeZ, zMin, zMax)

def Func(z):
	return 1.0 /(z * z)

def FuncGrad(z):
	return -2.0 / (z * z * z)


for iz in range(sizeZ):
	z = grid1D.getGridZ(iz)
	val = Func(z)
	grid1D.setValue(val,iz)

diff_max  = 0.
diffS_max = 0.

max_point  = z
max_pointS = z

for iz in range(sizeTestp):
	z = zMin + iz * dz
	func0 = Func(z)
	func  = grid1D.getValue(z)
	funcS = grid1D.getValueSmoothed(z)
	grad0 = FuncGrad(z)
	grad  = grid1D.calcGradient(z)
	gradS = grid1D.calcGradientSmoothed(z)
	
	diff  = math.sqrt((1.0 - grad  / grad0)**2)
	diffS = math.sqrt((1.0 - gradS / grad0)**2)

	print(z, "   ", func0, "   ", func, "   ", funcS, "   ", grad0, "   ", grad, "   ", gradS, "   ", diff, "   ", diffS)

	if(diff  > diff_max):
		diff_max   = diff
		max_point  = z

	if(diffS > diffS_max):
		diffS_max  = diffS
		max_pointS = z


print("max diff  =",diff_max)
print("max point =",max_point)
z = max_point
print("Func at max point  =",Func(z))
print("Grad at max point  =",FuncGrad(z))

print("max diffS  =",diffS_max)
print("max pointS =",max_pointS)
z = max_pointS
print("Func at max pointS =",Func(z))
print("Grad at max pointS =",FuncGrad(z))

print("Stop.")
