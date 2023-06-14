#----------------------------------------------------------------------
#Grid3D wrapped bunch binning test
#It is a test when we want to bin bunch with a certain length L
#into the Grid3D with z-length smaller than L. It is used 
#in the case of neighboring bunches space charge effect calculations.
#----------------------------------------------------------------------

import sys
import math
import random
import orbit.core

from spacecharge import Grid3D
from bunch import Bunch

print("==== Start ====")

sizeX = 64
sizeY = 64
sizeZ = 64

xMin = -1.0
xMax = +1.0
yMin = -1.0
yMax = +1.0
zMin = -1.0
zMax = +1.0

grid3D = Grid3D(sizeX,sizeY,sizeZ)
grid3D.setGridX(xMin,xMax)
grid3D.setGridY(yMin,yMax)
grid3D.setGridZ(zMin,zMax)

grid3D.longWrapping(True)

def getBunch(nParts,sigma):
	b = Bunch()
	b.macroSize(1.0e-6)
	z_min = +1.0e+64
	z_max = -1.0e+64
	for i in range(nParts):
		x  = 0.
		xp = 0.
		y  = 0.
		yp = 0.
		zp = 0.
		z = random.gauss(0,sigma)
		if(z < z_min): z_min = z
		if(z > z_max): z_max = z
		b.addParticle(x,xp,y,yp,z,zp)
	b.compress()
	return (b,z_min,z_max)


labmda = (zMax - zMin)


def getRhoArr(grid3D,sigma,labmda,nParts = 1000000):
	(bunch,z_min,z_max) = getBunch(nParts,sigma)
	grid3D.setZero()
	grid3D.binBunch(bunch,labmda)
	sizeX = grid3D.getSizeX()
	sizeY = grid3D.getSizeY()
	sizeZ = grid3D.getSizeZ()
	longRho_arr = []
	z_arr = []
	max_val = 0.
	sum_total = 0.
	for ind_z in range(sizeZ):
		sum_rho = 0.
		z = grid3D.getGridZ(ind_z)
		for ind_x in range(sizeX):
			for ind_y in range(sizeY):
				val = grid3D.getValueOnGrid(ind_x,ind_y,ind_z)
				#print "debug (ind_x,ind_y,ind_z)=",(ind_x,ind_y,ind_z)," val = ",val
				sum_rho += val
		longRho_arr.append(sum_rho)
		z_arr.append(z)
		sum_total += sum_rho
		if(max_val < sum_rho): max_val = sum_rho
	print("debug sigma =",sigma," sum_total (should be the same)=",sum_total)
	return [longRho_arr,z_arr,max_val]

#---- longitudinal sigma values for the bunch 
sigma_arr = [0.1,0.3,0.5,0.8,1.0,2.0,3.0]

results_arr = []
for sigma in sigma_arr:
	[longRho_arr,z_arr,max_val] = getRhoArr(grid3D,sigma,labmda)
	results_arr.append([longRho_arr,z_arr,max_val])
	
#---- now let's print the results
fl_out = open("wrapped_z_distr.dat","w")

st = "z "
for sigma in sigma_arr:
	st += " S"+str(sigma)+" "
	
fl_out.write(st + "\n")

z_arr = results_arr[0][1]

for ind_z in range(len(z_arr)):
	z = z_arr[ind_z]
	st = " %8.6f "%z
	for ind_s in range(len(sigma_arr)):
			st += " %12.5g "%results_arr[ind_s][0][ind_z]
	fl_out.write(st + "\n")
	
fl_out.close()
	

import matplotlib.pyplot as plt

for ind_s in range(len(sigma_arr)):
	plt.plot(z_arr,results_arr[ind_s][0])

plt.ylabel('RhoZ')
plt.xlabel('z')
#plt.axis([-1., 1., 0, 0.00026])

plt.savefig('rho_wrapped_plot_unnormalized.png')
plt.show()

#------- perform normalizations of all results to 1
for [longRho_arr,z_arr,max_val] in results_arr:
	for ind in range(sizeZ):
		longRho_arr[ind] /= max_val	
#--------------------------------------------------

for ind_s in range(len(sigma_arr)):
	plt.plot(z_arr,results_arr[ind_s][0])

plt.ylabel('RhoZ')
plt.xlabel('z')
#plt.axis([-1., 1., 0, 0.00026])

plt.savefig('rho_wrapped_plot.png')
plt.show()

print("Stop.")

