#-----------------------------------------------------
# Grid1D interpolation test
# User can change number of grid points (sizeZ) and the type of interpolation 
# (grid1D.getValue or grid1D.getValueSmoothed)
# The test should be performed for sizeZ = 1,2,3,  20
# You expect to see good agreement between theoretical and grid1D values
#-----------------------------------------------------

import sys
import math
import orbit.core

from spacecharge import Grid1D


print("Start.")

sizeZ = 20
zMin =  1.0
zMax = +4.0

extent = 2.0

grid1D = Grid1D(sizeZ,zMin,zMax)


def FuncTest(z):
	y = 2.0 + 0.5*(z-1.) + z**2
	return y
	
def FuncTestGrad(z):
	y = 0.5 + 2*z
	return y	

#---------------------------------
# Setup values on the grid
#---------------------------------
z_arr = []
val_arr = []
for ind in range(sizeZ):
	z = grid1D.getGridZ(ind)
	val = FuncTest(z)
	grid1D.setValue(val,ind)
	z_arr.append(z)
	val_arr.append(val)
	
#-------------------------------------------
# Interpolation data preparation
#-------------------------------------------
nPlotPoints = 1000
x_polt_min = zMin - extent
x_polt_max = zMax + extent
step = (x_polt_max - x_polt_min)/(nPlotPoints - 1)

x_arr = []
y_arr = []
yp_th_arr = []
yp_arr = []
for ind_p in range(nPlotPoints):
	x = x_polt_min + ind_p*step
	y = grid1D.getValue(x)
	#y = grid1D.getValueSmoothed(x)
	yp_th = FuncTestGrad(x)
	yp = grid1D.calcGradient(x)
	#yp = grid1D.calcGradientSmoothed(x)
	x_arr.append(x)
	y_arr.append(y)
	yp_th_arr.append(yp_th)
	yp_arr.append(yp)
	
	
#-------------------------------------
# Plots
#-------------------------------------
	
import matplotlib.pyplot as plt	
	
plt.plot(x_arr,y_arr)
plt.plot(z_arr,val_arr,"ro")
plt.ylabel('y(x)')
plt.xlabel('x')

plt.show()

plt.plot(x_arr,yp_arr)
plt.plot(x_arr,yp_th_arr)
plt.ylabel("y'(x)")
plt.xlabel('x')

plt.show()
