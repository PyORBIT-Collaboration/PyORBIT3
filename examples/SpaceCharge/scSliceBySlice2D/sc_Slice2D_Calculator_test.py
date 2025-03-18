#-----------------------------------------------------
# This example calculates the field gradients for bunch
# with particles uniformly distributed inside the circle 
# in X-Y plane for each z and with triangle density 
# distribution along z-axis. 
# The space charge is calculated by using Slice2D SC node.
# It means field potential for each X-Y slice solved as
# 2D problem. The longitudinal field can be calculated 
# from the potential of each X-Y slice as z-dependent function
# 
# The fields are compared with the analytical formula.
#
# Keep in mind that 3D density in Grid3D treated as 
# long charge threads at particular X,Y position when
# the 2D Poisson problem solved for each slice. So, when
# you interpret potential you have to account for z_step
# accordingly.
#-----------------------------------------------------

import sys
import math
import random

from orbit.core.bunch import Bunch

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.teapot import teapot

from orbit.core.spacecharge import Boundary2D

from orbit.space_charge.sc2dslicebyslice import scAccNodes, scLatticeModifications
from orbit.core.spacecharge import SpaceChargeCalcSliceBySlice2D

#------------------------------
# Let's make everything random
#------------------------------
random.seed(10)

#-------------------------------------
# Longitudinal distribution functions
#-------------------------------------

def getLongZ_Uniform(bunch_length):
	""" Uniform distribution alond z-axis """
	z = bunch_length*0.5*(1.0 - 2*random.random())
	return z
	
def getLongZ_Triangle(bunch_length):
	""" Triangle distribution alond z-axis """
	z_is_good = False
	while(not z_is_good):
		z = bunch_length*0.5*(1.0 - 2*random.random())
		threshold = 1. - 2.*abs(z)/bunch_length
		if(random.random() > threshold):
			continue
		z_is_good = True
	return z

print ("Start.")
#---------------------------------
#make an initial bunch
#---------------------------------
b_init = Bunch()
bunch_radius = 0.005        # m
bunch_length = 10.0       # m
nParts = 3000000
total_macroSize = 1.0e+14
energy = 1.4               # GeV

#---- In PoissonSolverFFT2D.cc we use our own definition of the field and potential
#---- for charged cylinder Er = lambda/r instead of CGS standard 
#---- Er = 2*lambda/r outside cylinder.
#---- Inside Er=lambda*r/(r.cylinder)**2 instead of Er=2*lambda*r/(r.cylinder)**2

gradXY_theory_CGS_standart = 2*(total_macroSize/bunch_length)/bunch_radius

#---- Our definition in PyORBIT
gradXY_theory = gradXY_theory_CGS_standart/2

#------------------------------------------------
#generate an uniform cylinder distribution
#no energy or momentum spread
#------------------------------------------------
for ip in range(nParts):
	r = bunch_radius*math.sqrt(random.random())
	phi = 2*math.pi*random.random()
	x = r*math.sin(phi)
	y = r*math.cos(phi)
	#z = getLongZ_Uniform(bunch_length)
	z = getLongZ_Triangle(bunch_length)
	b_init.addParticle(x,0.,y,0.,z,0.)

#-------------------------------
# set bunch parameters
#-------------------------------
nParticlesGlobal = b_init.getSizeGlobal()
b_init.macroSize(total_macroSize/nParticlesGlobal)
syncPart = b_init.getSyncParticle()
syncPart.kinEnergy(energy)

#-------------------------------
# Let's make boundray
#-------------------------------

nboundarypoints = 128
n_freespacemodes = 32
#---- boundary diameter
d_boundary = 0.020
boundary = Boundary2D(nboundarypoints, n_freespacemodes, "Circle", d_boundary, d_boundary)

#----------------------------------------------
#    make Slice2D space charge calculator
#----------------------------------------------
sizeX = 64   #number of grid points in horizontal direction
sizeY = 64   #number of grid points in vertical direction
sizeZ = 15   #number of longitudinal slices in the 2.5D space charge solver
calcsliced = SpaceChargeCalcSliceBySlice2D(sizeX, sizeY, sizeZ)


#-------------------------------------------------------
# Track bunch through Slice2D Space Charge solver node
# We are not going to analyse the bunch particles in
# this script. We just need the potential phiGrid3d
# that was generated during the tracking.
#-------------------------------------------------------
sc_distance = 1.
#---- Longitudinal field is correct only when the boundary is present
#---- Transverse is always working 
#calcsliced.trackBunch(b_init,sc_distance)
calcsliced.trackBunch(b_init,sc_distance,boundary)

#---- N.B. after calcsliced.trackBunch(...) the sizes and values of
#---- rhoGrid3d changed to get correct physical values for potential
#---- in phiGrid3d.
#---- The bin longitudinal size of rhoGrid3d (and accordingly in phiGrid3d)
#---- increased by gamma for Lorentz transformation from lab. system to bunch
#---- rest system of coordinates.
#---- After 3D bunch binning into rhoGrid3d we have to multiply each rho value
#---- by 2/grid_step_size_z because we have to transform the charge at all
#---- point of rhoGrid3d from point-like charge to thin thread charge per 
#---- meter. We have to do this because we will solve 2D problem for each X-Y 
#---- slice as a set of longitudinal charged threads. Also we have to account 
#---- for 2D solver Green function been lambda/r (it was used for simplicity :( ) 
#---- instead of 2*lambda/ as it should be in the CGS system of units.

rhoGrid3d = calcsliced.getRhoGrid()
phiGrid3d = calcsliced.getPhiGrid()

print ("==================================")
print ("rhoGrid3d min-max X= %+8.5f  %+8.5f"%(rhoGrid3d.getMinX(), rhoGrid3d.getMaxX()))
print ("rhoGrid3d min-max Y= %+8.5f  %+8.5f"%(rhoGrid3d.getMinY(), rhoGrid3d.getMaxY()))
print ("rhoGrid3d min-max Z= %+8.5f  %+8.5f"%(rhoGrid3d.getMinZ(), rhoGrid3d.getMaxZ()))
print ("==================================")

#---- array of macrozise of each slice
rho_slice_arr = []
for iz in range(rhoGrid3d.getSizeZ()):
	rho_slice = 0.
	for ix in range(rhoGrid3d.getSizeX()):
		for iy in range(rhoGrid3d.getSizeY()):
			rho_slice +=  rhoGrid3d.getValueOnGrid(ix,iy,iz)
	rho_slice_arr.append(rho_slice)

fl_out_name = "slice2d_sc_calculator_test.dat"
fl_out = open(fl_out_name,"w")

st = "z_ind   z[m]  rho_slice r/2  gradXY  gradXY_th   gradZ  gradZ_th   phi0"
print (st)
fl_out.write(st + "\n")

z_step = rhoGrid3d.getGridZ(1) -  rhoGrid3d.getGridZ(0)

rho_total = 0.

#---- we use it to calculate d(lambda)/dz = (rho_slice(z+z_step) - rho_slice(z))/z_step
rho_slice_old = 0.
for iz in range(rhoGrid3d.getSizeZ()):
	z_grid = rhoGrid3d.getGridZ(iz)
	rho_slice = rho_slice_arr[iz]
	#---------------------------------------------------------------------
	#---- gradZ_theory should not depend on number of slices and z_step
	#---- The 0.5 coefficient should be used because rho_slice was corrected
	#---- during the bunch tracking to accout for SC 2D solver with lambda/r
	#---- Green function instead of 2*lambda/r in CSG system of units.
	gradZ_theory = 0.5*((rho_slice - rho_slice_old)/(z_step))*(1+2.0*math.log(d_boundary/(2*bunch_radius)))
	#---- gradXY_theory should not depend on number of slices and z_step
	gradXY_theory = 0.5*(2*rho_slice/bunch_radius)
	#---------------------------------------------------------------------
	rho_slice_old = rho_slice
	#----------------------------
	x = 0.
	y = 0.
	ez0 = phiGrid3d.calcGradient(x,y,z_grid)[2]
	gradZ = ez0
	phi0 = phiGrid3d.getValue(x,y,z_grid)
	theta = math.pi/4
	x = (bunch_radius/2)*math.cos(theta)
	y = (bunch_radius/2)*math.sin(theta)
	r = math.sqrt(x**2 + y**2)
	(ex,ey,ez) = phiGrid3d.calcGradient(x,y,z_grid)
	gradXY = math.sqrt(ex**2 + ey**2)
	st  = "%2d "%iz + " %+8.4f "%z_grid + " %12.5g "%rho_slice
	st += "   %+8.4f  %12.5g  %12.5g  "%(r,gradXY, gradXY_theory*(r/bunch_radius))
	st += " %+12.5g "%gradZ + " %+12.5g "%gradZ_theory + " %+12.5g "%phi0
	print (st)
	fl_out.write(st + "\n")
	rho_total += rho_slice

fl_out.close()

print ("=============================================================")
print ("This charge is different form total_macroSize = %12.5g "%total_macroSize)
print ("rho_total = %12.5g "%rho_total)
print ("Reasons for that are explaned in the text of the script.")
print ("Corrected with (z_step/2) rho_total = %12.5g "%(rho_total*z_step/2))
print ("=============================================================")

print ("Stop.")
sys.exit(0)


