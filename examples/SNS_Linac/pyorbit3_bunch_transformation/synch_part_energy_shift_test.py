#! /usr/bin/env python

"""
This script will track one particle in the bunch through the
drift + RF Gap + drift lattice. The absolute position and the energy
of the partcle should be the same after the synchronous particle energy
shifted by some amount.
"""

import sys
import math
import random
import time

from orbit.py_linac.lattice import LinacAccLattice
from orbit.py_linac.lattice import LinacAccNodes

from orbit.py_linac.lattice import BaseLinacNode, LinacNode, LinacMagnetNode, MarkerLinacNode, Drift, Quad, AbstractRF_Gap, Bend
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice import RF_Cavity, Sequence
from orbit.py_linac.lattice import BaseRF_Gap

# from linac import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF
from linac import BaseRfGap_slow, RfGapTTF_slow

from bunch import Bunch
from bunch import SynchPartRedefinitionZdE

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.lattice import LinacPhaseApertureNode


# ---- auxiliary function for the node array sorting according to the positions
def positionComp(node1_da, node2_da):
    if node1_da.getParam("pos") > node2_da.getParam("pos"):
        return 1
    else:
        if node1_da.getParam("pos") == node2_da.getParam("pos"):
            return 0
    return -1


# ---- auxiliary function to print bunch
def printBunch(bunch):
    n_parts = bunch.getSize()
    eKin = bunch.getSyncParticle().kinEnergy()
    print("synch part Ekin [MeV] =  %16.9g " % (eKin * 1000.0))
    print("synch part time [us]  =  %16.9g " % (bunch.getSyncParticle().time() * 1.0e6))
    for ind in range(n_parts):
        (x, xp, y, yp, z, dE) = (bunch.x(ind), bunch.xp(ind), bunch.y(ind), bunch.yp(ind), bunch.z(ind), bunch.dE(ind))
        print("i= %3d       (z,dE) [mm,MeV] = ( %16.9g , %16.9g )" % (ind, z * 1000.0, (eKin + dE) * 1000.0))
    print("     --- distances ---")
    for ind in range(n_parts - 1):
        (x, xp, y, yp, z, dE) = (bunch.x(ind), bunch.xp(ind), bunch.y(ind), bunch.yp(ind), bunch.z(ind), bunch.dE(ind))
        (x1, xp1, y1, yp1, z1, dE1) = (
            bunch.x(ind + 1),
            bunch.xp(ind + 1),
            bunch.y(ind + 1),
            bunch.yp(ind + 1),
            bunch.z(ind + 1),
            bunch.dE(ind + 1),
        )
        print("i= %3d  delta(z,dE) [mm,MeV] = ( %16.9g , %16.9g )" % (ind, (z - z1) * 1000.0, (dE - dE1) * 1000.0))


# ------------------------------------------------------
# ------------------- Start of Script ------------------
# ------------------------------------------------------
L_drift_1 = 1.0
L_drift_2 = 1.0
seq_start_pos = 0.0

accLattice = LinacAccLattice("LinacDriftRFGapDrift")

# ---- lattice has one accelerator sequence
accSeq = Sequence("TEST")
accSeq.setLinacAccLattice(accLattice)
accSeq.setLength(L_drift_1 + L_drift_2)
accSeq.setPosition(seq_start_pos)

bpmFrequency = 805.0e6
accSeq.addParam("bpmFrequency", bpmFrequency)

# ---- lattice has one RF cavity with one RF Gap
# ---- Cavity is not an accelerator node in the lattice
cav_amp = 1.0
frequency = 402.5e6
cav_pos = L_drift_1
cav = RF_Cavity("RF_Cavity")
cav.setAmp(cav_amp)
cav.setFrequency(frequency)
cav.setPosition(cav_pos)
accSeq.addRF_Cavity(cav)

# ---- make RF Gap node
E0TL = 0.000100  # GeV
gap_phase = -90.0
rfGap_node = BaseRF_Gap("RF_Cavity:RF_Gap")
rfGap_node.setLength(0.0)
rfGap_node.setParam("E0TL", E0TL)
rfGap_node.setParam("E0L", E0TL)
rfGap_node.setParam("mode", 0.0)
rfGap_node.setParam("gap_phase", gap_phase * math.pi / 180.0)
rfGap_node.setParam("EzFile", "")
rfGap_node.setParam("beta_min", 0.01)
rfGap_node.setParam("beta_max", 0.99)
rfGap_node.setParam("pos", L_drift_1)
(polyT, polyS, polyTp, polySp) = rfGap_node.getTTF_Polynimials()
polyT.order(0)
polyS.order(0)
polyTp.order(0)
polySp.order(0)
polyT.coefficient(0, 1.0)
polyS.coefficient(0, 0.0)
polyTp.coefficient(0, 0.0)
polySp.coefficient(0, 0.0)


# rfGap_node.setCppGapModel(BaseRfGap_slow())
rfGap_node.setCppGapModel(BaseRfGap())
# rfGap_node.setCppGapModel(RfGapTTF())
# rfGap_node.setCppGapModel(RfGapTTF_slow())

# ---- set the RF Gap as a part of the cavity
cav.addRF_GapNode(rfGap_node)
cav.setPhase(gap_phase * math.pi / 180.0)

# ---- drift nodes
drift1 = Drift("drifr1")
drift2 = Drift("drifr2")
drift1.setLength(L_drift_1)
drift2.setLength(L_drift_2)
drift1.setParam("pos", L_drift_1 / 2)
drift2.setParam("pos", L_drift_1 + L_drift_2 / 2.0)

# ---- add nodes to the sequence: the order of nodes will be defined by positions
accSeq.addNode(drift2)
accSeq.addNode(rfGap_node)
accSeq.addNode(drift1)
accSeq.getNodes().sort(positionComp)

# add all AccNodes to the linac lattice
for accNode in accSeq.getNodes():
    accLattice.addNode(accNode)
accLattice.initialize()


# ---- setup the the linac specific trackers
accLattice.setLinacTracker(True)

# for accNode in accLattice.getNodes():
# 	print "debug node=",accNode.getName()," L=",accNode.getLength()," pos=",accNode.getParam("pos")

# ---- let's make bunch with one particle
bunch = Bunch()
# set H- mass
# self.bunch.mass(0.9382723 + 2*0.000511)
bunch.mass(0.939294)
bunch.charge(-1.0)
bunch.getSyncParticle().kinEnergy(0.0025)

bunch.addParticle(0.000, 0.000, 0.000, 0.000, 0.001, 0.00005)
bunch.addParticle(0.000, 0.000, 0.000, 0.000, -0.001, -0.00005)

# set up design
accLattice.trackDesignBunch(bunch)


for delta_E in (0.0, 0.00001):
    # delta_E = 0.000
    bunch_tmp = Bunch()
    bunch.copyBunchTo(bunch_tmp)
    print("======================================================== CASE")
    print("delta_E [MeV] = ", delta_E * 1000.0)
    redefCalc = SynchPartRedefinitionZdE()
    redefCalc.shift_dE(bunch_tmp, delta_E)

    print("===== bunch at start =====")
    printBunch(bunch_tmp)

    accLattice.trackBunch(bunch_tmp)

    print("===== bunch at end =====")
    printBunch(bunch_tmp)

print("==============================================")
print("==============================================")

for delta_Z in (0.0, 0.001):
    bunch_tmp = Bunch()
    bunch.copyBunchTo(bunch_tmp)
    print("======================================================== CASE")
    print("delta_Z [mm] = ", delta_Z * 1000.0)
    redefCalc = SynchPartRedefinitionZdE()
    redefCalc.shift_Z(bunch_tmp, delta_Z)

    print("===== bunch at start =====")
    printBunch(bunch_tmp)

    accLattice.trackBunch(bunch_tmp)

    print("===== bunch at end =====")
    printBunch(bunch_tmp)
