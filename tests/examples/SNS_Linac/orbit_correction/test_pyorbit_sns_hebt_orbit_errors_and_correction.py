#! /usr/bin/env python

"""

This is an example of TrajectoryCorrection usage.
This script will correct orbit in HEBT1-HEBT2
by using the existing DCorrectors.

The orbit distortion is created by the offsets and rotations
of quads.

"""

import random
import pytest
import os

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from bunch import Bunch

from orbit.py_linac.orbit_correction import TrajectoryCorrection

from orbit.py_linac.lattice_modifications import CoordinateDisplacementNodesModification
from orbit.py_linac.lattice_modifications import StraightRotationZ_NodesModification
from orbit.py_linac.lattice_modifications import StraightRotationX_NodesModification
from orbit.py_linac.lattice_modifications import StraightRotationY_NodesModification
from orbit.py_linac.lattice_modifications import QuadFieldsErrorsDeployment

random.seed(100)
script_dir = os.path.dirname(__file__)

names = ["HEBT1", "HEBT2"]
# ---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.5)

# ---- the XML file name with the structure
xml_file_name = os.path.join(script_dir, "../sns_linac_xml/sns_linac.xml")

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


quads = accLattice.getQuads()


# -----Let's create bunc
eKin = 1.0

bunch_init = Bunch()
syncPart = bunch_init.getSyncParticle()
# set H- mass
# self.bunch_init.mass(0.9382723 + 2*0.000511)
bunch_init.mass(0.939294)
bunch_init.charge(-1.0)
syncPart.kinEnergy(eKin)

# ---- add one particle to the bunch
x0 = 0.0
xp0 = 0.0
y0 = 0.0
yp0 = 0.0
z0 = 0.0
dE = 0.0

bunch_init.addParticle(x0, xp0, y0, yp0, z0, dE)

bunch = Bunch()
bunch_init.copyBunchTo(bunch)

# ---- cat off levels for Gaussian distribution in sigma values
cut_off_level = 3.0

# ---- setup random offsets for quads
offset_error = 0.0001
coordDisplModification = CoordinateDisplacementNodesModification()
coordDisplModification.addLatticeNodes(quads)
coordDisplModification.setGaussDistributedDisplacementParameter("dx", offset_error, cut_off_level)
coordDisplModification.setGaussDistributedDisplacementParameter("dy", offset_error, cut_off_level)

# ---- setup random roration angles
angle_rms = 0.0001
# ---- setup random X-axis roration angles
angleModificationX = StraightRotationX_NodesModification()
angleModificationX.addLatticeNodes(quads)
angleModificationX.setGaussDistributedAngle(angle_rms, cut_off_level)
# ---- setup random Y-axis roration angles
angleModificationY = StraightRotationY_NodesModification()
angleModificationY.addLatticeNodes(quads)
angleModificationY.setGaussDistributedAngle(angle_rms, cut_off_level)
# ---- setup random Z-axis roration angles
angleModificationZ = StraightRotationZ_NodesModification()
angleModificationZ.addLatticeNodes(quads)
angleModificationZ.setGaussDistributedAngle(angle_rms, cut_off_level)

# ---- setup random quads fields
relative_error = 0.001
quadFieldModification = QuadFieldsErrorsDeployment()
quadFieldModification.addQuads(quads)
quadFieldModification.setGaussDistributedRealtiveErrors(relative_error, cut_off_level)

# set up design
accLattice.trackDesignBunch(bunch)
accLattice.trackBunch(bunch)

trajCorrection = TrajectoryCorrection(accLattice)

print("===============Initial Trajectory  =============")
bunch = Bunch()
bunch_init.copyBunchTo(bunch)
trajCorrection.calculateTrajectory(bunch, print_info=True)
print("======================================")


# ---- correct trajectory
bunch = Bunch()
bunch_init.copyBunchTo(bunch)
trajCorrection.correctTrajectory(bunch)

# ---- print the DCH and DCV fields
print("Maximal achivable field in DCH(V) in RTBT1 is 0.017 [T] ")
dcvs = trajCorrection.getDCVs()
dcvs_string = ""
for dc in dcvs:
    dcvs_string += f"DCV = {dc.getName()} B [T] = {dc.getParam('B'): 6.5f}\n"
print(dcvs_string)
print("==============")
dchs = trajCorrection.getDCHs()
dchs_string = ""
for dc in dchs:
    dchs_string += f"DCH = {dc.getName()} B [T] = {dc.getParam('B'): 6.5f}\n"
print(dchs_string)
print("===============Corrected Trajectory=============")
bunch = Bunch()
bunch_init.copyBunchTo(bunch)
trajCorrection.calculateTrajectory(bunch, print_info=True)
print("======================================")

print("Done.")


def test_dcvs():
    expected_dcvs_string = """DCV = HEBT_Mag:DCV05 B [T] =  0.00067
DCV = HEBT_Mag:DCV07 B [T] =  0.00026
DCV = HEBT_Mag:DCV15 B [T] = -0.00043
DCV = HEBT_Mag:DCV17 B [T] =  0.00034
DCV = HEBT_Mag:DCV21 B [T] =  0.00069
DCV = HEBT_Mag:DCV23 B [T] = -0.00047
DCV = HEBT_Mag:DCV29 B [T] = -0.00142
DCV = HEBT_Mag:DCV31 B [T] =  0.00104\n"""
    assert dcvs_string == expected_dcvs_string


def test_dchs():
    expected_dcvs_string = """DCH = HEBT_Mag:DCH06 B [T] =  0.00123
DCH = HEBT_Mag:DCH08 B [T] = -0.00010
DCH = HEBT_Mag:DCH14 B [T] =  0.00001
DCH = HEBT_Mag:DCH16 B [T] = -0.00033
DCH = HEBT_Mag:DCH22 B [T] =  0.00041
DCH = HEBT_Mag:DCH24 B [T] = -0.00119
DCH = HEBT_Mag:DCH28 B [T] =  0.00041
DCH = HEBT_Mag:DCH30 B [T] =  0.00027\n"""
    assert dchs_string == expected_dcvs_string
