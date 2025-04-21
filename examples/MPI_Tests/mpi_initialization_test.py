#-------------------------------------------------------------
# This is a test that MPI is initialized and Bunch class
# instance can feel this MPI
#-------------------------------------------------------------
"""
>mpirun -np 2 python mpi_initialization_test.py
"""

import sys
import math
import random

from orbit.core.orbit_mpi import mpi_comm, MPI_Comm_rank, MPI_Comm_size, MPI_Bcast
from orbit.core.orbit_mpi import mpi_datatype, mpi_op
from orbit.core.orbit_mpi import MPI_Bcast
from orbit.core import orbit_mpi

from orbit.core.bunch import Bunch


mpi_init = orbit_mpi.MPI_Initialized()

rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
size = orbit_mpi.MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)

if(rank == 0):
	print("debug mpi is initialized =", mpi_init, " should be not zero.")
	print ("debug there should only one line like this rank=",rank," N CPUs = ",size)

n_particles = 10

bunch = Bunch()
for ind in range(n_particles):
	bunch.addParticle(0.,0.,0.,0.,0.,0.)

n_particles_global = bunch.getSizeGlobal()
if(rank == 0):
	print ("debug should be only one: rank=",rank," n_particles_global =",n_particles_global," it should be =",n_particles*size)

if(rank == 0):
	print ("debug should be only one Stop.")

sys.exit(mpi_init)
