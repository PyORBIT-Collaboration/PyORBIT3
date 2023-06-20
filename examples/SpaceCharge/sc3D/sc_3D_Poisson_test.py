#-----------------------------------------------------
#Creates Grid3D for charge density and Poisson Solver
#for a pointlike charge (3D case) 
#Copmares with the exact result.
#-----------------------------------------------------
# We can scale the grids with respect to the solver by
# a constant coefficient. It is result of scaling of the Green 
# function which is 1/r

import sys
import math
import orbit.core

import orbit_mpi

from spacecharge import Grid3D
from spacecharge import PoissonSolverFFT3D

print("Start.")

sizeX = 64
sizeY = 64
sizeZ = 64
xMin = -5.0
xMax = +5.0
yMin = -5.5
yMax = +5.5
zMin = -6.0
zMax = +6.0
print(" x,y,z sizes: ",sizeX," ",sizeY," ",sizeZ)
solver = PoissonSolverFFT3D(sizeX,sizeY,sizeZ,xMin,xMax,yMin,yMax,zMin,zMax)

scale_coeff = 1.2

gridRho = Grid3D(sizeX,sizeY,sizeZ)
gridRho.setGridX(scale_coeff*xMin,scale_coeff*xMax)
gridRho.setGridY(scale_coeff*yMin,scale_coeff*yMax)
gridRho.setGridZ(scale_coeff*zMin,scale_coeff*zMax)

gridPhi = Grid3D(sizeX,sizeY,sizeZ)
gridPhi.setGridX(scale_coeff*xMin,scale_coeff*xMax)
gridPhi.setGridY(scale_coeff*yMin,scale_coeff*yMax)
gridPhi.setGridZ(scale_coeff*zMin,scale_coeff*zMax)

chrage_pos_x = 1.0
chrage_pos_y = 1.0
chrage_pos_z = 1.0
charge = 1.0
gridRho.binValue(charge,chrage_pos_x,chrage_pos_y,chrage_pos_z)

print("Start solver.")

solver.findPotential(gridRho,gridPhi)

r_test = 4.0
n_angle_steps = 10
angle_step = 360./(n_angle_steps - 1)
print("  i    x       y        z         r          phi         phi_theory    ratio phi/theory  ")
for i in range(n_angle_steps):
	angle = math.pi*i*angle_step/180.
	x = r_test*math.cos(angle)
	y = r_test*math.sin(angle)
	z = 1.
	phi = gridPhi.getValue(x,y,z)
	dist = math.sqrt((chrage_pos_x - x)*(chrage_pos_x - x) + (chrage_pos_y - y)*(chrage_pos_y - y) + (chrage_pos_z - z)*(chrage_pos_z - z))
	phi_th = 1./dist
	print("",i," %7.4f  %7.4f  %7.4f  %7.4f  %12.5g  %12.5g  %12.7g  "%(x,y,z,dist,phi,phi_th,(phi/phi_th))) 

print("Stop.")

count = 0
while(1 < 2):
	solver.findPotential(gridRho,gridPhi)
	count += 1
	if(count % 10 == 0): print("solved n=",count)

