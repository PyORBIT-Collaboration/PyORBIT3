#! /usr/bin/env python

"""

This is an example of TrajectoryCorrection usage.
This script will correct orbit in HEBT1-HEBT2
by using the existing DCorrectors.

"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from orbit.core.bunch import Bunch

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.orbit_correction import TrajectoryCorrection

names = ["HEBT1", "HEBT2"]
# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.5)

# ---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print("Linac lattice is ready. L=", accLattice.getLength())

# ----------------------------------------------------------
# Set Linac style quads and drifts instead of TEAPOT style
# That can be useful when energy spread is huge and TEAPOT
# accuracy is not enough for tracking.
# This will slow down tracking and it is not symplectic.
# ----------------------------------------------------------
# accLattice.setLinacTracker(True)


# -----Let's create bunc
eKin = 1.0

bunch = Bunch()
syncPart = bunch.getSyncParticle()
# set H- mass
# self.bunch.mass(0.9382723 + 2*0.000511)
bunch.mass(0.939294)
bunch.charge(-1.0)
syncPart.kinEnergy(eKin)

# ---- add one particle to the bunch
x0 = 1.0 / 1000.0
xp0 = 1.0 / 1000.0
y0 = -1.0 / 1000.0
yp0 = -1.0 / 1000.0

z0 = 0.0
dE = 0.0

bunch.addParticle(x0, xp0, y0, yp0, z0, dE)

bunch_init = Bunch()
bunch.copyBunchTo(bunch_init)

# set up design
accLattice.trackDesignBunch(bunch_init)

print("Design tracking completed.")

# track through the lattice
paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.1}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")

pos_start = 0.0

file_out = open("pyorbit_hebt_distorted_trajectory.dat", "w")

st = " Node   position   x[mm]   xp[mrad]   y[mm]   yp[mrad]  "
file_out.write(st + "\n")
print(st)


def action_entrance(paramsDict):
    node = paramsDict["node"]
    bunch = paramsDict["bunch"]
    pos = paramsDict["path_length"]
    if paramsDict["old_pos"] == pos:
        return
    if paramsDict["old_pos"] + paramsDict["pos_step"] > pos:
        return
    paramsDict["old_pos"] = pos
    paramsDict["count"] += 1
    x = bunch.x(0) * 1000.0
    y = bunch.y(0) * 1000.0
    xp = bunch.xp(0) * 1000.0
    yp = bunch.yp(0) * 1000.0
    st = " %35s  %4.5f " % (node.getName(), pos + pos_start)
    st += "   %6.4f  %6.4f  %6.4f  %6.4f" % (x, xp, y, yp)
    file_out.write(st + "\n")
    file_out.flush()
    print(st)


def action_exit(paramsDict):
    action_entrance(paramsDict)


actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

time_start = time.process_time()

accLattice.trackBunch(bunch_init, paramsDict=paramsDict, actionContainer=actionContainer)

time_exec = time.process_time() - time_start
print("time[sec]=", time_exec)

file_out.close()


trajCorrection = TrajectoryCorrection(accLattice)

"""
bpms = trajCorrection.getBPMs()
for bpm in bpms:
	print "bpm=",bpm.getName()

print "==============================="

dchs = trajCorrection.getDCHs()
for dch in dchs:
	print "dch=",dch.getName()

print "==============================="

dcvs = trajCorrection.getDCVs()
for dcv in dcvs:
	print "dcv=",dcv.getName()

print "==============================="

quads = trajCorrection.getQuads()
for quad in quads:
	print "quad=",quad.getName()

print "==============================="

trnsvBPMs = trajCorrection.getTransverseBPMs()
for trnsvBPM in trnsvBPMs:
	print "TransverseBPM = ",trnsvBPM.getName()

print "==============================="

trnsvBPMs = trajCorrection.getQuadTransverseBPMs()
for trnsvBPM in trnsvBPMs:
	print "TransverseBPM = ",trnsvBPM.getName()
"""

print("=====================================")
dch_arr = trajCorrection.getDCHs()
dcv_arr = trajCorrection.getDCVs()

print("==============Initial Fields in Correctors ")
for dch in dch_arr:
    print("dch=", dch.getName(), " field [T] = %+8.7f" % dch.getParam("B"))
print("========")
for dcv in dcv_arr:
    print("dcv=", dcv.getName(), " field [T] = %+8.7f" % dcv.getParam("B"))

print("===============Trajectory=============")
# ----- print initial trajectory
bunch_init = Bunch()
bunch.copyBunchTo(bunch_init)
trajCorrection.calculateTrajectory(bunch_init, print_info=True)

# ---- correct trajectory
bunch_init = Bunch()
bunch.copyBunchTo(bunch_init)
trajCorrection.correctTrajectory(bunch_init)

print("==============New Fields in Correctors ")
for dch in dch_arr:
    print("dch=", dch.getName(), " field [T] = %+8.7f" % dch.getParam("B"))
print("========")
for dcv in dcv_arr:
    print("dcv=", dcv.getName(), " field [T] = %+8.7f" % dcv.getParam("B"))

print("===============Trajectory=============")
# ----- print fixed trajectory
bunch_init = Bunch()
bunch.copyBunchTo(bunch_init)
trajCorrection.calculateTrajectory(bunch_init, print_info=True)
print("=====================================")

print("Done.")
