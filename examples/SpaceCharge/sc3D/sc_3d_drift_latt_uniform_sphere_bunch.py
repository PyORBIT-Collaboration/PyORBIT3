#---------------------------------------------------------
# This example tracks the particles uniformly distributed 
# inside the sphere in the bunch rest coordinate system.
# After tracking through the drift one or several times 
# the ratios of x,y,z rms sizes should be the same.
# Here we use 3D Uniform Ellipse Poisson Solver.
#---------------------------------------------------------
import argparse
import sys
import math
import random
from time import perf_counter

from orbit.core.orbit_mpi import mpi_comm, MPI_Comm_rank, MPI_Comm_size, MPI_Bcast
from orbit.core.orbit_mpi import mpi_datatype, mpi_op
from orbit.core.orbit_mpi import MPI_Bcast
from orbit.core import orbit_mpi

from orbit.core.bunch import Bunch, BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.teapot import teapot

from orbit.core.spacecharge import SpaceChargeCalc3D
from orbit.space_charge.sc3d import setSC3DAccNodes

from orbit.space_charge.sc3d import setUniformEllipsesSCAccNodes
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse

from orbit.space_charge.sc2p5d import setSC2p5DrbAccNodes
from orbit.core.spacecharge import SpaceChargeCalc2p5Drb



def print0(*args, **kwargs):
	if orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD) == 0:
		print(*args, **kwargs)

#------------------------------
# Let's make everything random
#------------------------------
random.seed(10)

#------------------------------
# Auxilary functions
#------------------------------

def func_theory(x):
	"""
	The function used in the charged sphere expansion dynamics calculation.
	See the presentation.
	"""
	val = math.sqrt(abs(x**2-1))
	val = math.log(x+val)/2 + x*val/2
	return val
	
def getR_Theory(s_path,r_init,R,Ntot,bunch):
	"""
	The root finding in the charged sphere expansion dynamics calculation.
	See the presentation.
	s_path - path length / drift length []
	r_init - initial distance from the sphere center [m]
	R - bunch sphere radius [m]
	Ntot - total macro-size of the bunch
	beta - relativistic beta
	r_classic - classical radius of particle
	"""
	r_classic = bunch.classicalRadius()
	gamma = bunch.getSyncParticle().gamma()
	beta = bunch.getSyncParticle().beta()
	val = math.sqrt((r_classic*Ntot/R**3)/2)*s_path/(gamma*beta)
	#---- Now we solve equation val = func_theory(x) for x 
	x_min = 1.000000001
	x_max = 10.
	x = (x_min + x_max)/2.
	#print ("debug s_path=",s_path," val=",val)
	n_iter = 15
	count = 0
	f_min = func_theory(x_min)
	f_max = func_theory(x_max)
	f_val = func_theory((x_min + x_max)/2.)
	while(count < n_iter):
		x = (x_min + x_max)/2.
		f_val = func_theory(x)
		if(val <= f_val):
			f_max = f_val
			x_max = x
			count += 1
			continue
		if(val >= f_val):
			f_min = f_val
			x_min = x
			count += 1
			continue
	r_res = 1000.*r_init*x**2
	#print ("debug s_path=",s_path," val=",val," f_val=",f_val," x=",x," r_res=",r_res)
	return r_res

#----------------------------------------
# Uniform sphere distribution functions
#----------------------------------------

def getUniformXYZ(radius,gamma):
	""" 
	Uniform distribution inside the sphere. 
	Z-direction is contracted by  relativistic gamma
	"""
	(x,y,z) = (radius,radius,radius)
	r = math.sqrt(x**2 + y**2 + z**2)
	while(r > radius):
		x = radius*(1.0 - 2*random.random())
		y = radius*(1.0 - 2*random.random())
		z = radius*(1.0 - 2*random.random())
		r = math.sqrt(x**2 + y**2 + z**2)
	return (x,y,z/gamma)



def getLattice(lattice_length,n_parts,sc_type):
	elem = teapot.DriftTEAPOT("my_drift")
	elem.setLength(lattice_length)
	elem.setnParts(n_parts)
	teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
	teapot_lattice.addNode(elem)
	teapot_lattice.initialize()
	# we will put SC nodes as frequently as possible.
	# In this case one for each part of the Drift node
	sc_path_length_min = 0.1
	if sc_type == 'ellipsoid':
		nEllipses = 1
		calc3d = SpaceChargeCalcUnifEllipse(nEllipses)
		scNodes_arr = setUniformEllipsesSCAccNodes(teapot_lattice, sc_path_length_min, calc3d)
	elif sc_type == 'fft3d':
		sizeX = 64
		sizeY = 64
		sizeZ = 64
		calc3d = SpaceChargeCalc3D(sizeX, sizeY, sizeZ)
		scNodes_arr = setSC3DAccNodes(teapot_lattice, sc_path_length_min, calc3d)
	#print ("debug n sc nodes=",len(scNodes_arr))
	return teapot_lattice


def main(nParts: int, sc_type: str):
	print0("Start.")
	#---------------------------------
	#make an initial bunch
	#---------------------------------
	b_init = Bunch()
	bunch_radius = 0.005        # m
	bunch_length = 10.0       # m
	# nParts = 1000000
	total_macroSize = 1.0e+11
	energy = 1.4               # GeV

	syncPart = b_init.getSyncParticle()
	syncPart.kinEnergy(energy)
	beta = syncPart.beta()
	gamma = syncPart.gamma()
	print0("==================================================")
	print0("Mass[MeV]           = %8.3f "%(b_init.mass()*1000.))
	print0("Kinetic energy[MeV] = %8.3f "%(energy*1000.))
	print0("gamma               = %8.5f "%gamma)
	print0("Total macrosize     = %12.5e "%total_macroSize)
	print0("==================================================")


	rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
	mpi_size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)

	#------------------------------------------------
	#---- generate a uniform sphere distribution
	#---- without energy or momentum spread
	#------------------------------------------------

	r_in = 0.5 * bunch_radius
	if rank == 0:
		# -------------------------------------------------------------------
		# ---- Here we add three particles - along x,y,z axis at
		# ---- r_bunch/2 distance from the center. We will use them to
		# ---- compare analytical solution with the PyORBIT tracking,
		# ---- N.B. : adding 3 particles did not change total charge much.
		# ------------------------------------------------------------------
		b_init.addParticle(-r_in, 0., 0., 0., 0., 0.)
		b_init.addParticle(0., 0., r_in, 0., 0., 0.)
		b_init.addParticle(0., 0., 0., 0., r_in / gamma, 0.)

	r_in *= 1000  # now in [mm]


	cr = 0
	for ip in range(nParts):
		(x,y,z) = getUniformXYZ(bunch_radius,gamma)
		if rank == cr % mpi_size:
			b_init.addParticle(x,0.,y,0.,z,0.)
		cr += 1

	#-------------------------------
	# set bunch parameters
	#-------------------------------
	nParticlesGlobal = b_init.getSizeGlobal()
	b_init.macroSize(total_macroSize/nParticlesGlobal)
	print0("Total numbers of macro-particles = ",nParticlesGlobal)

	#----------------------------------------------
	#---- Analysis of the initial bunch
	#----------------------------------------------
	twiss_analysis = BunchTwissAnalysis()
	twiss_analysis.analyzeBunch(b_init)

	x_rms = math.sqrt(twiss_analysis.getCorrelation(0,0)) * 1000.0
	y_rms = math.sqrt(twiss_analysis.getCorrelation(2,2)) * 1000.0
	z_rms = math.sqrt(twiss_analysis.getCorrelation(4,4)) * 1000.0

	print0("==========Before Traking=====================")
	print0("Initial x_rms[mm] = %6.5f  x_rms/y_rms         = %6.5f "%(x_rms,x_rms/y_rms))
	print0("Initial y_rms[mm] = %6.5f  x_rms/(z_rms*gamma) = %6.5f "%(y_rms,x_rms/(z_rms*gamma)))
	print0("Initial z_rms[mm] = %6.5f  x_rms/(z_rms*gamma) = %6.5f "%(z_rms,x_rms/(z_rms*gamma)))
	print0("=============================================")
	#------------------------------------------
	#          Initial Bunch is ready
	#          Now let's make a lattice
	#------------------------------------------



	# -------------------------------------------------------------
	# set of uniformly charged ellipses Space Charge
	# -------------------------------------------------------------


	lattice_length = 0.5    # the length of the drift
	n_parts = 30  # number of parts on what the drift will be chopped, or the number of SC nodes

	lattice = getLattice(lattice_length, n_parts, sc_type)

	print0("Lattice length[m] = %6.3f "%lattice_length)
	print0("Number of parts in the lattice = ",n_parts)

	#--------------------------------------------------------------------
	# Create a new bunch, copy everything from b_init to this new bunch
	#--------------------------------------------------------------------
	bunch = Bunch()
	b_init.copyBunchTo(bunch)


	#-------------------------------------------------------
	#---- Let's track bunch through the drift several times
	#-------------------------------------------------------
	s_path = 0.

	print0("s[m] r_rms[mm] r1[mm]  r2[mm] r3[mm] r_theory[mm] ")
	start = perf_counter()
	for ind in range(20):
		#---------------------------------------
		#track the bunch through the lattice
		#---------------------------------------
		lattice.trackBunch(bunch)

		s_path += lattice_length

		#----------------------------------------------
		#---- Analysis of the final bunch
		#----------------------------------------------
		twiss_analysis.analyzeBunch(bunch)

		x_rms = math.sqrt(twiss_analysis.getCorrelation(0,0)) * 1000.0
		y_rms = math.sqrt(twiss_analysis.getCorrelation(2,2)) * 1000.0
		z_rms = math.sqrt(twiss_analysis.getCorrelation(4,4)) * 1000.0
		r_rms = math.sqrt(x_rms**2 + y_rms**2 + (z_rms*gamma)**2)

		st = ""
		r_out = 0.
		for ip in range(3):
			z = bunch.z(ip)*gamma
			x = bunch.x(ip)
			y = bunch.y(ip)
			xp = bunch.xp(ip)*1000
			r_out = 1000.*math.sqrt(x**2 + y**2 + z**2)
			st += " %6.3f "%(r_out)


		r_out_th = getR_Theory(s_path,0.5*bunch_radius,bunch_radius,total_macroSize,bunch)

		print0(" %7.3f  %6.3f  "%(s_path,r_rms),st," %6.3f "%r_out_th)



	finish = perf_counter()
	print0("==========After Tracking=======================")
	print0("x_rms[mm] = %6.5f  x_rms/y_rms         = %6.5f "%(x_rms,x_rms/y_rms))
	print0("y_rms[mm] = %6.5f  y_rms/(z_rms*gamma) = %6.5f "%(y_rms,y_rms/(z_rms*gamma)))
	print0("z_rms[mm] = %6.5f  x_rms/(z_rms*gamma) = %6.5f "%(z_rms,x_rms/(z_rms*gamma)))
	print0("=============================================")

	# print0(f"Finish run on {mpi_size} CPUs.\n"
	# 	   f"Tracking Wall Time {(finish - start):3f} seconds.\n"
	# 	   f"Tracking rate {nParticlesGlobal/(finish - start)/mpi_size:1f} particles/second/node")
	print0("Final: %7.3f  %6.3f  " % (s_path, r_rms), st, " %6.3f " % r_out_th)
	sys.exit(0)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Track uniformly distributed charged sphere.')

	parser.add_argument('--N', type=int, default=1000000, help='Number of particles (default: 1000000)')
	parser.add_argument('--SC', type=str, default='ellipsoid', choices=['ellipsoid', 'fft3d', 'none'],
						help='Type of space charge (default: ellipsoid)')
	args = parser.parse_args()


	main(args.N, args.SC)
