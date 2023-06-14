#!/usr/bin/env python

#--------------------------------------------------------
# The classes will generates bunches for pyORBIT SNS linac 
# at the entrance of SNS MEBT accelerator line (by default)
# It is parallel, but it is not efficient.
#--------------------------------------------------------

import math
import sys
import os
import random

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import KVDist2D, KVDist3D
from orbit.bunch_generators import GaussDist2D, GaussDist3D
from orbit.bunch_generators import WaterBagDist2D, WaterBagDist3D
from orbit.bunch_generators import TwissAnalysis

from bunch import Bunch

class SNS_Linac_BunchGenerator:
	"""
	Generates the pyORBIT SNS Linac Bunches.
	Twiss parameters has the fol following units: x in [m], xp in [rad]
	and the X and Y emittances are un-normalized. The longitudinal emittance 
	is in [GeV*m].
	"""
	def __init__(self,twissX, twissY, twissZ, frequency = 402.5e+6):
		self.twiss = (twissX, twissY, twissZ)
		self.bunch_frequency = frequency
		self.bunch = Bunch()
		syncPart = self.bunch.getSyncParticle()
		#set H- mass
		#self.bunch.mass(0.9382723 + 2*0.000511)
		self.bunch.mass(0.939294)
		self.bunch.charge(-1.0)
		syncPart.kinEnergy(0.0025)
		self.c = 2.99792458e+8    # speed of light in m/sec
		self.beam_current = 38.0  # beam current in mA , design = 38 mA
		self.rf_wave_lenght = self.c/self.bunch_frequency
		self.si_e_charge = 1.6021773e-19
		
	def getKinEnergy(self):
		"""
		Returns the kinetic energy in GeV
		"""
		return self.bunch.getSyncParticle().kinEnergy()
		
	def setKinEnergy(self, e_kin = 0.0025):
		"""
		Sets the kinetic energy in GeV
		"""
		self.bunch.getSyncParticle().kinEnergy(e_kin)
		
	def getZtoPhaseCoeff(self,bunch):
		"""
		Returns the coefficient to calculate phase in degrees from the z-coordinate.
		"""
		bunch_lambda = bunch.getSyncParticle().beta()*self.rf_wave_lenght 
		phase_coeff = 360./bunch_lambda
		return phase_coeff
		
	def getBeamCurrent(self):
		"""
		Returns the beam currect in mA
		"""
		return self.beam_current
		
	def setBeamCurrent(self, current):
		"""
		Sets  the beam currect in mA
		"""
		self.beam_current = current
	
	def getBunch(self, nParticles = 0, distributorClass = WaterBagDist3D, cut_off = -1.):
		"""
		Returns the pyORBIT bunch with particular number of particles.
		"""
		comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
		rank = orbit_mpi.MPI_Comm_rank(comm)
		size = orbit_mpi.MPI_Comm_size(comm)		
		data_type = mpi_datatype.MPI_DOUBLE		
		main_rank = 0		
		bunch = Bunch()
		self.bunch.copyEmptyBunchTo(bunch)		
		macrosize = (self.beam_current*1.0e-3/self.bunch_frequency)
		macrosize /= (math.fabs(bunch.charge())*self.si_e_charge)
		distributor = None
		if(distributorClass == WaterBagDist3D):
			distributor = distributorClass(self.twiss[0],self.twiss[1],self.twiss[2])
		else:
			distributor = distributorClass(self.twiss[0],self.twiss[1],self.twiss[2], cut_off)
		bunch.getSyncParticle().time(0.)	
		for i in range(nParticles):
			(x,xp,y,yp,z,dE) = distributor.getCoordinates()
			(x,xp,y,yp,z,dE) = orbit_mpi.MPI_Bcast((x,xp,y,yp,z,dE),data_type,main_rank,comm)
			if(i%size == rank):
				bunch.addParticle(x,xp,y,yp,z,dE)
		nParticlesGlobal = bunch.getSizeGlobal()
		bunch.macroSize(macrosize/nParticlesGlobal)
		return bunch
	

