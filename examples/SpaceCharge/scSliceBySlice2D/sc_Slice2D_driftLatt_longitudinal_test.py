#-----------------------------------------------------
# This example tracks the particles distributed 
# with triangle shaped longitudinal density and 
# uniformly distributed inside the circle in transverse 
# direction.
# The space charge is calculated by using Slice2D SC node.
# The change in the particles transverse and longitudinal 
# momentum compared with the analytical formulas.
# User can see that results in transverse plane is not sensitive 
# to the presence of the beam pipe boundary.
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
nParts = 1000000
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
beta = syncPart.beta()
gamma = syncPart.gamma()

#------------------------------------------
#          Initial Bunch is ready
#
#          Now let's make a lattice
#------------------------------------------

nboundarypoints = 128
n_freespacemodes = 32
# d_boundary diameter
d_boundary = 0.040
boundary = Boundary2D(nboundarypoints, n_freespacemodes, "Circle", d_boundary, d_boundary)

def getLattice(lattice_length,n_parts,calcsliced):
	elem = teapot.DriftTEAPOT("my_drift")
	elem.setLength(lattice_length)
	elem.setnParts(n_parts)	
	teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
	teapot_lattice.addNode(elem)
	teapot_lattice.initialize()
	# we will put SC nodes as frequently as possible.
	# In this case one for each part of the Drift node
	sc_path_length_min = 0.1
	scNodes_arr = scLatticeModifications.setSC2DSliceBySliceAccNodes(teapot_lattice, sc_path_length_min,calcsliced,boundary)
	#print ("debug n sc nodes=",len(scNodes_arr))
	return teapot_lattice

#----------------------------------------------
#make Slice2D space charge calculator
#----------------------------------------------
sizeX = 64   #number of grid points in horizontal direction
sizeY = 64   #number of grid points in vertical direction
sizeZ = 20  #number of longitudinal slices in the 2.5D space charge solver
calcsliced = SpaceChargeCalcSliceBySlice2D(sizeX, sizeY, sizeZ)

#---- We are going to account for longitudinal electric SC field 
res = calcsliced.longTracking(True)

lattice_length = 1.0    # the length of the drift
n_parts = 1  # number of parts on what the drift will be chopped, or the number of SC nodes
lattice = getLattice(lattice_length,n_parts,calcsliced)

#-------------------------------
#  Lattice is ready
#-------------------------------
"""
nodes = lattice.getNodes()
for node in nodes:
	print ("node=", node.getName()," s start,stop = %4.3f %4.3f "%lattice.getNodePositionsDict()[node])
	print ("There are ", node.getNumberOfBodyChildren()," child nodes.")
"""
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
#  The theory dE = (Ldrift/gamma^2)*(1+2*ln(R/r0))*r_classic*(4*Nmacro/l_bunch^2)*mass
#  r0 is a classical radius of the particle
#---------------------------------------
r_classic = b.classicalRadius()
mass = b.mass()
dE_th = (lattice_length/gamma**2)*(1. + 2*math.log(d_boundary/(2*bunch_radius)))*r_classic*(4*total_macroSize/bunch_length**2)*mass
dE_th *= 1000*1000 # in keV

n_graph_points = 100

#---- place to collect dE vs z z_dE_arr[[z,dE]]
z_dE_arr = []

#---- we will use the particles near the z-axis
#---- If number of particles is too small, you can increase the r_max
r_max = bunch_radius/100

print ("=============================================")

st = "index   z[m]   dE[keV]    dE_th[keV]  "
print (st)

n_points_count = 0
for ip in range(b.getSize()):
	z = b.z(ip)
	x = b.x(ip)
	y = b.y(ip)
	r = math.sqrt(x**2 + y**2)
	if(r > r_max): continue
	z = b.z(ip)
	dE = (b.dE(ip) - b_init.dE(ip))*1000.*1000.
	n_points_count += 1
	z_dE_arr.append([z,dE])
	#print ("debug i= %3d "%n_points_count, " z= %+8.4f  dE[keV] = %+13.6g "%(z,dE))
	if(n_points_count > n_graph_points): break

z_dE_arr = sorted(z_dE_arr, key= lambda val_arr: val_arr[0])

for ind,[z,dE] in enumerate(z_dE_arr):
	dE_th = math.copysign(abs(dE_th),z)
	st  = "%3d "%ind
	st += "%+8.4f  %+10.3g  %+10.3g "%(z,dE,dE_th)
	print (st)

print("Finish.")
sys.exit(0)
