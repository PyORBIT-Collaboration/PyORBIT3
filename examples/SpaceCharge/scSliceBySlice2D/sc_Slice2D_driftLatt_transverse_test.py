#-----------------------------------------------------
# This example tracks the particles distributed 
# uniformly in a cylinder through a drift.
# The space charge is calculated by using Slice2D SC node.
# The change in the particles transverse momentum compared
# with the analytical formula.
# User can see that results is not sensitive to the presence of
# the beam pipe boundary.
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

print ("Start.")
#---------------------------------
#make an initial bunch
#---------------------------------
b_init = Bunch()
bunch_radius = 0.005        # m
bunch_length = 10.0       # m
nParts = 3000000
total_macroSize = 1.0e+14
energy = 0.4               # GeV

#------------------------------------------------
#generate an uniform cylinder distribution
#no energy or momentum spread
#------------------------------------------------
for ip in range(nParts):
	r = bunch_radius*math.sqrt(random.random())
	phi = 2*math.pi*random.random()
	x = r*math.sin(phi)
	y = r*math.cos(phi)
	z = bunch_length*0.5*(1.0 - 2*random.random())
	b_init.addParticle(x,0.,y,0.,z,0.)

#-------------------------------
# set bunch parameters
#-------------------------------
nParticlesGlobal = b_init.getSizeGlobal()
b_init.macroSize(total_macroSize/nParticlesGlobal)
syncPart = b_init.getSyncParticle()
syncPart.kinEnergy(energy)


#--------------------------------------
#---- Dump ORBIT_MPI Bunch
#--------------------------------------
"""
from orbit.utils.orbit_mpi_utils.bunch_pyorbit_to_orbit import bunch_pyorbit_to_orbit
fl_bunch_orbit_name = "orbit_100k_cilinder.dat" 
ringLength = 348.333/2.0
bunch_pyorbit_to_orbit(ringLength,b_init,fl_bunch_orbit_name)
"""

#------------------------------------------
#          Initial Bunch is ready
#
#          Now let's make a lattice
#------------------------------------------

nboundarypoints = 128
n_freespacemodes = 32
# r_boundary seems here diameter!? In ORBIT it's 0.145m
r_boundary = 0.020
boundary = Boundary2D(nboundarypoints, n_freespacemodes, "Circle", r_boundary, r_boundary)

def getLattice(lattice_length,n_parts,calcsliced):
	elem = teapot.DriftTEAPOT("my_drift")
	elem.setLength(lattice_length)
	elem.setnParts(n_parts)	
	teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
	teapot_lattice.addNode(elem)
	teapot_lattice.initialize()
	# we will put SC nodes as frequently as possible.
	# In this case one for each part of the Drift node
	sc_path_length_min = 0.00000001
	scLatticeModifications.setSC2DSliceBySliceAccNodes(teapot_lattice, sc_path_length_min,calcsliced)
	#scLatticeModifications.setSC2DSliceBySliceAccNodes(teapot_lattice, sc_path_length_min,calcsliced,boundary)
	return teapot_lattice

#----------------------------------------------
#make Slice2D space charge calculator
#----------------------------------------------
sizeX = 64   #number of grid points in horizontal direction
sizeY = 64   #number of grid points in vertical direction
sizeZ = 3  #number of longitudinal slices in the 2.5D space charge solver
calcsliced = SpaceChargeCalcSliceBySlice2D(sizeX, sizeY, sizeZ)

lattice_length = 1.0    # the length of the drift
n_parts = 1  # number of parts on what the drift will be chopped, or the number of SC nodes
lattice = getLattice(lattice_length,n_parts,calcsliced)

#-------------------------------
#  Lattice is ready
#-------------------------------

nodes = lattice.getNodes()
for node in nodes:
	print ("node=", node.getName()," s start,stop = %4.3f %4.3f "%lattice.getNodePositionsDict()[node])
	print ("There are ", node.getNumberOfBodyChildren()," child nodes.")


#--------------------------------------------------------------------
# Create a new bunch, copy everything from b_init to this new bunch
# Later we will compare momentum kicks between two bunches
#--------------------------------------------------------------------
b = Bunch()
b_init.copyBunchTo(b)

#---------------------------------------
#track the bunch through the lattice
#---------------------------------------
lattice.trackBunch(b)

#---------------------------------------
#  Now compare bunches and the theory
#  The theory dp = r * 2 * r0 * total_macrosize * charge^2 * drift_length / (gamma^3*beta^2*bunch_radius^2*bunch_length)
#  r0 is a classical radius of the particle
#---------------------------------------
n_graph_points = 10000
n_step = int(b.getSize()/n_graph_points) + 1

st = " z[m]  r[m]  dp[rad] dp_theory[rad] "
print (st)

#--- r angle angle_theory
r_dp_arr = [[],[],[],[]]
for ip in range(0,b.getSize(),n_step):
	z = b.z(ip)
	x = b.x(ip)
	y = b.y(ip)
	xp = b.xp(ip)
	yp = b.yp(ip)
	x0 = b_init.x(ip)
	y0 = b_init.y(ip)
	xp0 = b_init.xp(ip)
	yp0 = b_init.yp(ip)	
	r = math.sqrt(x0**2+y0**2)
	dp = math.sqrt((xp - xp0)**2 + (yp - yp0)**2)
	dp_th = 2.0*r*b.classicalRadius()*total_macroSize*(b.charge())**2*lattice_length / (b.getSyncParticle().gamma()**3*b.getSyncParticle().beta()**2*bunch_radius**2*bunch_length)	
	print (" %+12.5e  %12.5e  %12.5e  %12.5e "%(z,r,dp,dp_th))
	r_dp_arr[0].append(z)
	r_dp_arr[1].append(r)
	r_dp_arr[2].append(dp)
	r_dp_arr[3].append(dp_th)
	
fl_out = open("test_slice2d_uniform_cylinder_result_slices_"+str(sizeZ) +".dat","w")

st = " z[m]  r[m]  dp[rad] dp_theory[rad] "
fl_out.write(st + "\n")

for ind in range(len(r_dp_arr[0])):
	z = r_dp_arr[0][ind]
	r = r_dp_arr[1][ind]
	dp = r_dp_arr[2][ind]
	dp_th = r_dp_arr[3][ind]
	st = "%+12.5e %12.5e %12.5e %12.5e"%(z,r,dp,dp_th)
	fl_out.write(st + "\n")
fl_out.close()



print("Finish.")
sys.exit()

