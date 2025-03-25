#---------------------------------------------------------
# This example tracks the particles uniformly distributed 
# inside the sphere in the bunch rest coordinate system.
# After tracking through the drift one or several times 
# the ratios of x,y,z rms sizes should be the same.
# Here we use 3D Uniform Ellipse Poisson Solver.
#---------------------------------------------------------

import sys
import math
import random

from orbit.core.bunch import Bunch, BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.teapot import teapot

from orbit.core.spacecharge import SpaceChargeCalc3D
from orbit.space_charge.sc3d import setSC3DAccNodes

from orbit.space_charge.sc3d import setUniformEllipsesSCAccNodes
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse

from orbit.space_charge.sc2p5d import setSC2p5DrbAccNodes
from orbit.core.spacecharge import SpaceChargeCalc2p5Drb

#------------------------------
# Let's make everything random
#------------------------------
random.seed(10)

#-------------------------------------
# Longitudinal distribution functions
#-------------------------------------

def getUniformXYZ(radius,gamma):
	""" 
	Uniform distribution inside the sphere. 
	Z-direction is contracted by  relativistic gamma
	"""
	(x,y,z) = (radius,radius,radius)
	r = math.sqrt(x**2 + y**2 + z**2)
	while(r > radius):
		x = radius*0.5*(1.0 - 2*random.random())
		y = radius*0.5*(1.0 - 2*random.random())
		z = radius*0.5*(1.0 - 2*random.random())
		r = math.sqrt(x**2 + y**2 + z**2)
	return (x,y,z/gamma)

print ("Start.")
#---------------------------------
#make an initial bunch
#---------------------------------
b_init = Bunch()
bunch_radius = 0.005        # m
bunch_length = 10.0       # m
nParts = 1000000
total_macroSize = 1.0e+11
energy = 1.0               # GeV

syncPart = b_init.getSyncParticle()
syncPart.kinEnergy(energy)
beta = syncPart.beta()
gamma = syncPart.gamma()
print ("==================================================")
print ("Mass[MeV]           = %8.3f "%(b_init.mass()*1000.))
print ("Kinetic energy[MeV] = %8.3f "%(energy*1000.))
print ("gamma               = %8.5f "%gamma)
print ("Total macrosize     = %12.5e "%total_macroSize)
print ("==================================================")

#------------------------------------------------
#---- generate an uniform sphere distribution
#---- without energy or momentum spread
#------------------------------------------------
for ip in range(nParts):
	(x,y,z) = getUniformXYZ(bunch_radius,gamma)
	b_init.addParticle(x,0.,y,0.,z,0.)

#-------------------------------
# set bunch parameters
#-------------------------------
nParticlesGlobal = b_init.getSizeGlobal()
b_init.macroSize(total_macroSize/nParticlesGlobal)
print ("Total numbers of macro-particles = ",nParticlesGlobal)

#----------------------------------------------
#---- Analysis of the initial bunch
#----------------------------------------------
twiss_analysis = BunchTwissAnalysis()
twiss_analysis.analyzeBunch(b_init)

x_rms = math.sqrt(twiss_analysis.getCorrelation(0,0)) * 1000.0
y_rms = math.sqrt(twiss_analysis.getCorrelation(2,2)) * 1000.0
z_rms = math.sqrt(twiss_analysis.getCorrelation(4,4)) * 1000.0

print ("==========Before Traking=====================")
print ("Initial x_rms[mm] = %6.5f  x_rms/y_rms         = %6.5f "%(x_rms,x_rms/y_rms))
print ("Initial y_rms[mm] = %6.5f  x_rms/(z_rms*gamma) = %6.5f "%(y_rms,x_rms/(z_rms*gamma)))
print ("Initial z_rms[mm] = %6.5f  x_rms/(z_rms*gamma) = %6.5f "%(z_rms,x_rms/(z_rms*gamma)))
print ("=============================================")
#------------------------------------------
#          Initial Bunch is ready
#          Now let's make a lattice
#------------------------------------------

def getLattice(lattice_length,n_parts,calc3d):
	elem = teapot.DriftTEAPOT("my_drift")
	elem.setLength(lattice_length)
	elem.setnParts(n_parts)	
	teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
	teapot_lattice.addNode(elem)
	teapot_lattice.initialize()
	# we will put SC nodes as frequently as possible.
	# In this case one for each part of the Drift node
	sc_path_length_min = 0.1
	scNodes_arr = setUniformEllipsesSCAccNodes(teapot_lattice, sc_path_length_min, calc3d)
	#print ("debug n sc nodes=",len(scNodes_arr))
	return teapot_lattice

# -------------------------------------------------------------
# set of uniformly charged ellipses Space Charge
# -------------------------------------------------------------
nEllipses = 1
calc3d = SpaceChargeCalcUnifEllipse(nEllipses)

lattice_length = 10.0    # the length of the drift
n_parts = 30  # number of parts on what the drift will be chopped, or the number of SC nodes
lattice = getLattice(lattice_length,n_parts,calc3d)

print ("Lattice length[m] = %6.3f "%lattice_length)
print ("Number of parts in the lattice = ",n_parts)

#--------------------------------------------------------------------
# Create a new bunch, copy everything from b_init to this new bunch
#--------------------------------------------------------------------
bunch = Bunch()
b_init.copyBunchTo(bunch)

#---------------------------------------
#track the bunch through the lattice
#---------------------------------------
lattice.trackBunch(bunch)

#----------------------------------------------
#---- Analysis of the final bunch
#----------------------------------------------
twiss_analysis.analyzeBunch(bunch)

x_rms = math.sqrt(twiss_analysis.getCorrelation(0,0)) * 1000.0
y_rms = math.sqrt(twiss_analysis.getCorrelation(2,2)) * 1000.0
z_rms = math.sqrt(twiss_analysis.getCorrelation(4,4)) * 1000.0

print ("==========After Traking=======================")
print ("Initial x_rms[mm] = %6.5f  x_rms/y_rms         = %6.5f "%(x_rms,x_rms/y_rms))
print ("Initial y_rms[mm] = %6.5f  x_rms/(z_rms*gamma) = %6.5f "%(y_rms,x_rms/(z_rms*gamma)))
print ("Initial z_rms[mm] = %6.5f  x_rms/(z_rms*gamma) = %6.5f "%(z_rms,x_rms/(z_rms*gamma)))
print ("=============================================")

print("Finish.")
sys.exit(0)
