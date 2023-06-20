#-----------------------------------------------------
# Creates Grid3D for charge density and potential. 
# Uses Poisson Solver to find potentials for 3D Gaussian
# distribution. The longitudinal binning is performed 
# assuming the longitudinally periodic distribution.
# It means we will use wrapped version of binning and gradient
# calculations.
#-----------------------------------------------------
# 

import sys
import math
import random
import orbit.core

import orbit_mpi

from spacecharge import Grid3D
from spacecharge import PoissonSolverFFT3D

from bunch import Bunch

print("==== Start ====")

#---- defines how to treat wrapping of the bunches
use_wrapping = True
#use_wrapping = False

sizeX = 64
sizeY = 64
sizeZ = 64

xMin = -1.0
xMax = +1.0
yMin = -1.0
yMax = +1.0
zMin = -1.0
zMax = +1.0

rho3D = Grid3D(sizeX,sizeY,sizeZ)
rho3D.setGridX(xMin,xMax)
rho3D.setGridY(yMin,yMax)
rho3D.setGridZ(zMin,zMax)

phi3D = Grid3D(sizeX,sizeY,sizeZ)
phi3D.setGridX(xMin,xMax)
phi3D.setGridY(yMin,yMax)
phi3D.setGridZ(zMin,zMax)

if(use_wrapping):
	print("Longitudinal wrapping for rho3D = ",rho3D.longWrapping())
	print("Longitudinal wrapping for phi3D = ",phi3D.longWrapping())
	
	rho3D.longWrapping(True)
	phi3D.longWrapping(True)
	
	print("Longitudinal wrapping for rho3D = ",rho3D.longWrapping())
	print("Longitudinal wrapping for phi3D = ",phi3D.longWrapping())

def getBunch(nParts,sigmaX,sigmaY,sigmaZ):
	b = Bunch()
	b.macroSize(1.0e-6)
	z_min = +1.0e+64
	z_max = -1.0e+64
	for i in range(nParts):
		x  = random.gauss(0,sigmaX)
		xp = 0.
		y  = random.gauss(0,sigmaY)
		yp = 0.
		zp = 0.
		z = random.gauss(0,sigmaZ)
		if(z < z_min): z_min = z
		if(z > z_max): z_max = z
		b.addParticle(x,xp,y,yp,z,zp)
	b.compress()
	return (b,z_min,z_max)


lambdaZ = (zMax - zMin)

nParts = 1000000
sigmaZ = 0.8

if(not use_wrapping):
	sigmaZ = 0.1

 
(sigmaX,sigmaY) = (0.1,0.1)
(bunch,z_min,z_max) = getBunch(nParts,sigmaX,sigmaY,sigmaZ)

rho3D.setZero()
phi3D.setZero()

if(use_wrapping):
	rho3D.binBunch(bunch,lambdaZ)
else:
	rho3D.binBunch(bunch)

#------------------------------------------
#  SC 3D Poisson solver
#------------------------------------------

print("Solver 3D grid (x,y,z) : ",(sizeX,sizeY,sizeZ))
solver = PoissonSolverFFT3D(sizeX,sizeY,sizeZ,xMin,xMax,yMin,yMax,zMin,zMax)

#---- this should be even number
nBunches = 2*10

if(use_wrapping):
	solver.numExtBunches(nBunches)
	solver.distBetweenBunches(lambdaZ)

#---- we have to update Green function FFT after we changed nBunches or lambda
solver.updateGeenFunction()

print("debug nBunches=",solver.numExtBunches())
print("debug distance=",solver.distBetweenBunches())
print("debug sigma_z =",sigmaZ)
print("debug ===================== Solver is ready")

solver.findPotential(rho3D,phi3D)

z_arr = []
rho_arr = []
phi_arr = []
Ez_arr = []
ind_x = sizeX/2
ind_y = sizeY/2
for ind_z in range(sizeZ):
	z = rho3D.getGridZ(ind_z)
	rho = rho3D.getValueOnGrid(ind_x,ind_y,ind_z)
	phi = phi3D.getValueOnGrid(ind_x,ind_y,ind_z)
	#rho = rho3D.getValue(0.,0.,z)
	#phi = phi3D.getValue(0.,0.,z)
	Ez = phi3D.calcGradient(0.,0.,z)[2]
	#print "debug ind_z=",ind_z," z=",z," rho =",rho," phi=",phi
	z_arr.append(z)
	rho_arr.append(rho)
	phi_arr.append(phi)
	Ez_arr.append(Ez)
	
import matplotlib.pyplot as plt

plt.plot(z_arr,rho_arr,"bo",z_arr,rho_arr)
plt.ylabel('RhoZ')
plt.xlabel('z')
#plt.axis([-1., 1., 0, 0.001])

plt.savefig('rho_plot.png')
plt.show()

plt.plot(z_arr,phi_arr,"bo",z_arr,phi_arr)
plt.ylabel('PhiZ')
plt.xlabel('z')
#plt.axis([-1., 1., 0., 9.0])

plt.savefig('phi_plot.png')
plt.show()

plt.plot(z_arr,Ez_arr,"bo",z_arr,Ez_arr)
plt.ylabel('Ez')
plt.xlabel('z')
#plt.axis([-1., 1., -10., 10.])

plt.savefig('ez_plot.png')
plt.show()

sys.exit(0)

