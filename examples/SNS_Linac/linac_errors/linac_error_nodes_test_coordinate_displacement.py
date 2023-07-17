#!/usr/bin/env python

# ------------------------------------------------------------------------
# Test of the linac error node - Coordinate Displacement Error Node
#
# Script creates a lattice with two quads, one dipole and three drifts.
#
# The generated bunch tracked through this lattice and we can see trajectory
# and correlations at the end
#
# The lattice is for the test only. It is nonphysical.
# ------------------------------------------------------------------------

import math
import sys
import os
import random

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer

from orbit.py_linac.lattice import LinacAccLattice
from orbit.py_linac.lattice import Sequence
from orbit.py_linac.lattice import Drift
from orbit.py_linac.lattice import Quad
from orbit.py_linac.lattice import Bend


from orbit.py_linac.errors import ErrorCntrlCoordDisplacement

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist3D
from orbit.bunch_generators import TwissAnalysis

from orbit.core.bunch import Bunch, BunchTwissAnalysis

random.seed(100)


def bunchCentering(bunch):
    """
    Bunch center after generating can have small deviation from the (0,0,0,0,0,0)
    This function will force centering the bunch.
    """
    twiss_analysis = BunchTwissAnalysis()
    twiss_analysis.analyzeBunch(bunch)
    # -----------------------------------------------
    # let's center the beam
    (x_avg, y_avg) = (twiss_analysis.getAverage(0), twiss_analysis.getAverage(2))
    (xp_avg, yp_avg) = (twiss_analysis.getAverage(1), twiss_analysis.getAverage(3))
    (z_avg, dE_avg) = (twiss_analysis.getAverage(4), twiss_analysis.getAverage(5))
    for part_id in range(bunch.getSize()):
        bunch.x(part_id, bunch.x(part_id) - x_avg)
        bunch.y(part_id, bunch.y(part_id) - y_avg)
        bunch.xp(part_id, bunch.xp(part_id) - xp_avg)
        bunch.yp(part_id, bunch.yp(part_id) - yp_avg)
        bunch.z(part_id, bunch.z(part_id) - z_avg)
        bunch.dE(part_id, bunch.dE(part_id) - dE_avg)
    # -----------------------------------------------
    return (x_avg, y_avg, dE_avg)


class SNS_Linac_BunchGenerator:
    """
    Generates the pyORBIT SNS Linac Bunches using the Gauss distribution.
    Twiss parameters have the following units: x in [m], xp in [rad]
    and the X and Y emittances are un-normalized. The longitudinal emittance
    is in [GeV*m].
    """

    def __init__(self, frequency=402.5e6):
        self.bunch_frequency = frequency
        # set H- mass
        # self.bunch.mass(0.9382723 + 2*0.000511)
        self.bunch = Bunch()
        self.init_coords = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.bunch.mass(0.939294)
        self.bunch.charge(-1.0)
        self.bunch.getSyncParticle().kinEnergy(0.0025)
        self.c = 2.99792458e8  # speed of light in m/sec
        self.beam_current = 38.0  # beam current in mA , design = 38 mA
        self.rf_wave_lenght = self.c / self.bunch_frequency
        self.si_e_charge = 1.6021773e-19
        # ----------------------------------------
        self.twiss_analysis = BunchTwissAnalysis()

    def setInitialCorrdsCenter(self, x0, xp0, y0, yp0, z0, dE0):
        self.init_coords = (x0, xp0, y0, yp0, z0, dE0)

    def getInitialCorrdsCenter(self):
        return self.init_coords

    def setParticleCharge(self, charge):
        """
        Sets the particle charge H- => -1.0   and proton => +1.0
        """
        self.bunch.charge(charge)

    def getKinEnergy(self):
        """
        Returns the kinetic energy in GeV
        """
        return self.bunch.getSyncParticle().kinEnergy()

    def setKinEnergy(self, e_kin=0.0025):
        """
        Sets the kinetic energy in GeV
        """
        self.bunch.getSyncParticle().kinEnergy(e_kin)

    def getZtoPhaseCoeff(self, bunch):
        """
        Returns the coefficient to calculate phase in degrees from the z-coordinate.
        """
        bunch_lambda = bunch.getSyncParticle().beta() * self.rf_wave_lenght
        phase_coeff = 360.0 / bunch_lambda
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

    def getBunch(self, nParticles, twissX, twissY, twissZ, cut_off=-1.0):
        """
        Returns the pyORBIT bunch with particular number of particles.
        """
        (x0, xp0, y0, yp0, z0, dE0) = self.init_coords
        bunch = Bunch()
        self.bunch.copyEmptyBunchTo(bunch)
        macrosize = self.beam_current * 1.0e-3 / self.bunch_frequency
        macrosize /= math.fabs(bunch.charge()) * self.si_e_charge
        distributor = GaussDist3D(twissX, twissY, twissZ, cut_off)
        bunch.getSyncParticle().time(0.0)
        for i in range(nParticles):
            (x, xp, y, yp, z, dE) = distributor.getCoordinates()
            bunch.addParticle(x + x0, xp + xp0, y + y0, yp + yp0, z + z0, dE + dE0)
        nParticlesGlobal = bunch.getSizeGlobal()
        bunch.macroSize(macrosize / nParticlesGlobal)
        return bunch


# ---------------------------------------------------------------
# ---- Let's make a linac lattice
# ---------------------------------------------------------------

accSeq = Sequence("Test_Errors_AccSeq")


drift_1 = Drift("drift_1")
drift_2 = Drift("drift_2")
drift_3 = Drift("drift_3")
drift_4 = Drift("drift_4")

for drift in [drift_1, drift_2, drift_3, drift_4]:
    drift.setLength(0.5)
    drift.setnParts(10)


quad_1 = Quad("quad_1")
quad_2 = Quad("quad_2")

quad_1.setLength(0.05)
quad_2.setLength(0.05)

quad_1.setParam("dB/dr", +15.0)
quad_2.setParam("dB/dr", -15.0)

quad_1.setnParts(6)
quad_2.setnParts(6)


bend = Bend("bend")
bend.setLength(0.5)
bend.setParam("theta", 90.0 * (math.pi / 180.0))
bend.setnParts(10)

# ---- This is for fringe field control of the bend magnet
# bend.setUsageFringeFieldIN(False)
# bend.setUsageFringeFieldOUT(False)

accSeq.addNode(drift_1)
accSeq.addNode(quad_1)
accSeq.addNode(drift_2)
accSeq.addNode(quad_2)
accSeq.addNode(drift_3)
accSeq.addNode(bend)
accSeq.addNode(drift_4)

accLattice = LinacAccLattice("Error_Test_Lattice")

for node in accSeq.getNodes():
    accLattice.addNode(node)
accLattice.initialize()

# ---- This will force the usage of linac type quad tracker
# ---- This tracker is defined in /scr/linac/tracking/linac_tracking.cc
# ---- If it is not True, the usual Teapot like quad tracking functions
# ---- will be used from /scr/teapot/teapotbase.cc
accLattice.setLinacTracker(True)

print("==============Lattice======================")
nodes = accLattice.getNodes()
for node in nodes:
    print("node =", node.getName())
print("===========================================")


# ---------------------------------------------------------------------------
# Here the example how to assign the error functions to a particular node
# in the lattice.
# ---------------------------------------------------------------------------

node_with_error = accLattice.getNodeForName("quad_1")

errorCntrl_1 = ErrorCntrlCoordDisplacement(node_with_error.getName() + "_errorNode")
errorCntrl_1.setOneNodeParent(node_with_error)
errorCntrl_1.setDisplacementParameter("dx", 0.01)  # in meters
errorCntrl_1.setDisplacementParameter("dy", 0.01)  # in meters
errorCntrl_1.setDisplacementParameter("dE", 0.0)  # in GeV


# -------------------------------------------------------
#    Now let's generate bunch - nothing new
# -------------------------------------------------------

bunch_generator = SNS_Linac_BunchGenerator()

peak_current = 0.0  # mA
bunch_generator.setBeamCurrent(peak_current)

# ------ PyORBIT emittances
(alphaX, betaX, emittX) = (-1.0, 4.0, 1.0 * 1.0e-6)
(alphaY, betaY, emittY) = (+1.0, 2.0, 1.0 * 1.0e-6)
(alphaZ, betaZ, emittZ) = (-0.02, 100.0, 0.016 * 1.0e-6)

print(" ========= PyORBIT Twiss ===========")
print(" aplha beta emitt[mm*mrad] X= %+6.4f %6.4f %6.4f " % (alphaX, betaX, emittX * 1.0e6))
print(" aplha beta emitt[mm*mrad] Y= %+6.4f %6.4f %6.4f " % (alphaY, betaY, emittY * 1.0e6))
print(" aplha beta emitt[m*GeV]   Z= %+6.4f %6.2f %6.4f " % (alphaZ, betaZ, emittZ * 1.0e6))

twissX = TwissContainer(alphaX, betaX, emittX)
twissY = TwissContainer(alphaY, betaY, emittY)
twissZ = TwissContainer(alphaZ, betaZ, emittZ)

nParticles = 50000

bunch = bunch_generator.getBunch(nParticles, twissX, twissY, twissZ)

(x_avg, y_avg, dE_avg) = bunchCentering(bunch)
print("debug (x_avg,y_avg,dE_avg) =", (x_avg, y_avg, dE_avg))

accLattice.trackDesignBunch(bunch)


# track through the lattice
paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.001}
actionContainer = AccActionsContainer("Bunch Tracking")

twiss_analysis = BunchTwissAnalysis()

pos_start = 0.0

print("===============================================================")
print(
    "     node              pos[m]        avg_x[mm] avg_y[mm] avg_dE[MeV]  x_rms[mm] y_rms[mm] dE_rms[MeV]   x_z_corr[mm*mm] x_dE_corr[mm*MeV]  y_z_corr[mm*mm] y_dE_corr[mm*MeV]     eKin[MeV]  nParts "
)


def actions_all(paramsDict):
    node = paramsDict["node"]
    bunch = paramsDict["bunch"]
    pos = paramsDict["path_length"]
    if paramsDict["old_pos"] == pos:
        return
    if paramsDict["old_pos"] + paramsDict["pos_step"] > pos:
        return
    paramsDict["old_pos"] = pos
    paramsDict["count"] += 1
    gamma = bunch.getSyncParticle().gamma()
    beta = bunch.getSyncParticle().beta()
    twiss_analysis.analyzeBunch(bunch)
    x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1] * twiss_analysis.getTwiss(0)[3]) * 1000.0
    y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1] * twiss_analysis.getTwiss(1)[3]) * 1000.0
    z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1] * twiss_analysis.getTwiss(2)[3]) * 1000.0
    dE_rms = math.sqrt(twiss_analysis.getTwiss(2)[2] * twiss_analysis.getTwiss(2)[3]) * 1000.0
    x_z_corr = twiss_analysis.getCorrelation(0, 4) * 1000.0 * 1000.0
    x_dE_corr = twiss_analysis.getCorrelation(0, 5) * 1000.0 * 1000.0
    y_z_corr = twiss_analysis.getCorrelation(2, 4) * 1000.0 * 1000.0
    y_dE_corr = twiss_analysis.getCorrelation(2, 5) * 1000.0 * 1000.0
    z_to_phase_coeff = bunch_generator.getZtoPhaseCoeff(bunch)
    z_rms_deg = z_to_phase_coeff * z_rms / 1000.0
    nParts = bunch.getSizeGlobal()
    x_avg = twiss_analysis.getAverage(0) * 1000.0
    y_avg = twiss_analysis.getAverage(2) * 1000.0
    z_deg_avg = twiss_analysis.getAverage(4) * z_to_phase_coeff
    dE_avg = twiss_analysis.getAverage(5) * 1.0e3
    eKin = bunch.getSyncParticle().kinEnergy() * 1.0e3
    s = " %20s  %4.5f " % (node.getName(), pos + pos_start)
    s += "   %8.2f  %8.2f  %8.5f   " % (x_avg, y_avg, dE_avg)
    s += "   %8.2f  %8.2f  %8.6f " % (x_rms, y_rms, dE_rms)
    s += "   %8.2f  %8.6f    " % (x_z_corr, x_dE_corr)
    s += "   %8.2f  %8.6f    " % (y_z_corr, y_dE_corr)
    s += "     %6.3f   %8d " % (eKin, nParts)
    print(s)


actionContainer.addAction(actions_all, AccActionsContainer.ENTRANCE)
actionContainer.addAction(actions_all, AccActionsContainer.BODY)
actionContainer.addAction(actions_all, AccActionsContainer.EXIT)

accLattice.trackBunch(bunch, paramsDict=paramsDict, actionContainer=actionContainer)

# --------------------------------------------------------------------------
# You can see the shift of average coordinates in the quad_1
# --------------------------------------------------------------------------

print("Done.")
