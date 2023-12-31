#!/usr/bin/env python

"""
This is not a parallel version!
"""

# for mpi operations
from orbit.core import orbit_mpi
from orbit.core.orbit_mpi import mpi_comm, mpi_datatype, mpi_op

import math
import random
import sys
from orbit.core import bunch

BunchTwissAnalysis = bunch.BunchTwissAnalysis
BunchTuneAnalysis = bunch.BunchTuneAnalysis
from ..utils.consts import speed_of_light


class StatLats:
    """
    This class gathers delivers the statistical twiss parameters
    """

    def __init__(self, filename):
        self.file_out = open(filename, "a")
        self.bunchtwissanalysis = BunchTwissAnalysis()

    def writeStatLats(self, s, bunch, lattlength=0):
        self.bunchtwissanalysis.analyzeBunch(bunch)
        emitx = self.bunchtwissanalysis.getEmittance(0)
        betax = self.bunchtwissanalysis.getBeta(0)
        alphax = self.bunchtwissanalysis.getAlpha(0)
        betay = self.bunchtwissanalysis.getBeta(1)
        alphay = self.bunchtwissanalysis.getAlpha(1)
        emity = self.bunchtwissanalysis.getEmittance(1)
        dispersionx = self.bunchtwissanalysis.getDispersion(0)
        ddispersionx = self.bunchtwissanalysis.getDispersionDerivative(0)
        dispersiony = self.bunchtwissanalysis.getDispersion(1)
        ddispersiony = self.bunchtwissanalysis.getDispersionDerivative(1)

        sp = bunch.getSyncParticle()
        time = sp.time()
        if lattlength > 0:
            time = sp.time() / (lattlength / (sp.beta() * speed_of_light))

        # if mpi operations are enabled, this section of code will
        # determine the rank of the present node
        rank = 0  # default is primary node
        mpi_init = orbit_mpi.MPI_Initialized()
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        if mpi_init:
            rank = orbit_mpi.MPI_Comm_rank(comm)

        # only the primary node needs to output the calculated information
        if rank == 0:
            self.file_out.write(
                str(s)
                + "\t"
                + str(time)
                + "\t"
                + str(emitx)
                + "\t"
                + str(emity)
                + "\t"
                + str(betax)
                + "\t"
                + str(betay)
                + "\t"
                + str(alphax)
                + "\t"
                + str(alphay)
                + "\t"
                + str(dispersionx)
                + "\t"
                + str(ddispersionx)
                + "\n"
            )

    def closeStatLats(self):
        self.file_out.close()


class StatLatsSetMember:
    """
    This class delivers the statistical twiss parameters
    """

    def __init__(self, file):
        self.file_out = file
        self.bunchtwissanalysis = BunchTwissAnalysis()

    def writeStatLats(self, s, bunch, lattlength=0):
        self.bunchtwissanalysis.analyzeBunch(bunch)
        emitx = self.bunchtwissanalysis.getEmittance(0)
        betax = self.bunchtwissanalysis.getBeta(0)
        alphax = self.bunchtwissanalysis.getAlpha(0)
        betay = self.bunchtwissanalysis.getBeta(1)
        alphay = self.bunchtwissanalysis.getAlpha(1)
        emity = self.bunchtwissanalysis.getEmittance(1)
        dispersionx = self.bunchtwissanalysis.getDispersion(0)
        ddispersionx = self.bunchtwissanalysis.getDispersionDerivative(0)
        # dispersiony = self.bunchtwissanalysis.getDispersion(1, bunch)
        # ddispersiony = self.bunchtwissanalysis.getDispersionDerivative(1, bunch)

        sp = bunch.getSyncParticle()
        time = sp.time()

        if lattlength > 0:
            time = sp.time() / (lattlength / (sp.beta() * speed_of_light))

        # if mpi operations are enabled, this section of code will
        # determine the rank of the present node
        rank = 0  # default is primary node
        mpi_init = orbit_mpi.MPI_Initialized()
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        if mpi_init:
            rank = orbit_mpi.MPI_Comm_rank(comm)

        # only the primary node needs to output the calculated information
        if rank == 0:
            self.file_out.write(
                str(s)
                + "\t"
                + str(time)
                + "\t"
                + str(emitx)
                + "\t"
                + str(emity)
                + "\t"
                + str(betax)
                + "\t"
                + str(betay)
                + "\t"
                + str(alphax)
                + "\t"
                + str(alphay)
                + "\t"
                + str(dispersionx)
                + "\t"
                + str(ddispersionx)
                + "\n"
            )

    def closeStatLats(self):
        self.file_out.close()

    def resetFile(self, file):
        self.file_out = file


class Moments:
    """
    This class delivers the beam moments
    """

    def __init__(self, filename, order, nodispersion, emitnorm):
        self.file_out = open(filename, "a")
        self.bunchtwissanalysis = BunchTwissAnalysis()
        self.order = order
        if nodispersion == False:
            self.dispterm = -1
        else:
            self.dispterm = 1

        if emitnorm == True:
            self.emitnormterm = 1
        else:
            self.emitnormterm = -1

    def writeMoments(self, s, bunch, lattlength=0):
        sp = bunch.getSyncParticle()
        time = sp.time()
        if lattlength > 0:
            time = sp.time() / (lattlength / (sp.beta() * speed_of_light))

        self.bunchtwissanalysis.computeBunchMoments(bunch, self.order, self.dispterm, self.emitnormterm)

        # if mpi operations are enabled, this section of code will
        # determine the rank of the present node
        rank = 0  # default is primary node
        mpi_init = orbit_mpi.MPI_Initialized()
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        if mpi_init:
            rank = orbit_mpi.MPI_Comm_rank(comm)

        # only the primary node needs to output the calculated information
        if rank == 0:
            self.file_out.write(str(s) + "\t" + str(time) + "\t")
            for i in range(0, self.order + 1):
                for j in range(0, i + 1):
                    self.file_out.write(str(self.bunchtwissanalysis.getBunchMoment(i - j, j)) + "\t")
            self.file_out.write("\n")

    def closeMoments(self):
        self.file_out.close()


class MomentsSetMember:
    """
    This class delivers the beam moments
    """

    def __init__(self, file, order, nodispersion, emitnorm):
        self.file_out = file
        self.order = order
        self.bunchtwissanalysis = BunchTwissAnalysis()
        if nodispersion == False:
            self.dispterm = -1
        else:
            self.dispterm = 1

        if emitnorm == True:
            self.emitnormterm = 1
        else:
            self.emitnormterm = -1

    def writeMoments(self, s, bunch, lattlength=0):
        sp = bunch.getSyncParticle()
        time = sp.time()

        if lattlength > 0:
            time = sp.time() / (lattlength / (sp.beta() * speed_of_light))

        self.bunchtwissanalysis.computeBunchMoments(bunch, self.order, self.dispterm, self.emitnormterm)

        # if mpi operations are enabled, this section of code will
        # determine the rank of the present node
        rank = 0  # default is primary node
        mpi_init = orbit_mpi.MPI_Initialized()
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        if mpi_init:
            rank = orbit_mpi.MPI_Comm_rank(comm)

        # only the primary node needs to output the calculated information
        if rank == 0:
            self.file_out.write(str(s) + "\t" + str(time) + "\t")
            for i in range(0, self.order + 1):
                for j in range(0, i + 1):
                    self.file_out.write(str(self.bunchtwissanalysis.getBunchMoment(i - j, j)) + "\t")
            self.file_out.write("\n")

    def resetFile(self, file):
        self.file_out = file


class BPMSignal:
    """
    This class delivers the average value for coordinate x and y
    """

    def __init__(self):
        self.bunchtwissanalysis = BunchTwissAnalysis()
        self.xAvg = 0.0
        self.yAvg = 0.0
        self.xpAvg = 0.0
        self.ypAvg = 0.0

    def analyzeSignal(self, bunch):
        self.bunchtwissanalysis.analyzeBunch(bunch)

        # if mpi operations are enabled, this section of code will
        # determine the rank of the present node
        rank = 0  # default is primary node
        mpi_init = orbit_mpi.MPI_Initialized()
        comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
        if mpi_init:
            rank = orbit_mpi.MPI_Comm_rank(comm)

        # only the primary node needs to output the calculated information
        if rank == 0:
            self.xAvg = self.bunchtwissanalysis.getAverage(0)
            self.xpAvg = self.bunchtwissanalysis.getAverage(1)
            self.yAvg = self.bunchtwissanalysis.getAverage(2)
            self.ypAvg = self.bunchtwissanalysis.getAverage(3)

    def getSignalX(self):
        return self.xAvg

    def getSignalXP(self):
        return self.xpAvg

    def getSignalY(self):
        return self.yAvg

    def getSignalYP(self):
        return self.ypAvg
