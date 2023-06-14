#! /usr/bin/env python

"""

This is an example of TrajectoryCorrection usage.
This script will correct orbit in HEBT1-HEBT2
by using the existing DCorrectors.

The orbit distortion is created by the offsets and rotations
of quads.

"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from bunch import Bunch

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.orbit_correction import TrajectoryCorrection

from orbit.py_linac.lattice_modifications import CoordinateDisplacementNodesModification
from orbit.py_linac.lattice_modifications import StraightRotationZ_NodesModification
from orbit.py_linac.lattice_modifications import StraightRotationX_NodesModification
from orbit.py_linac.lattice_modifications import StraightRotationY_NodesModification
from orbit.py_linac.lattice_modifications import QuadFieldsErrorsDeployment

random.seed(100)


names = ["HEBT1","HEBT2"]
#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.5)

#---- the XML file name with the structure
xml_file_name = "../sns_linac_xml/sns_linac.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

print("Linac lattice is ready. L=",accLattice.getLength())

#----------------------------------------------------------
# Set Linac style quads and drifts instead of TEAPOT style 
# That can be useful when energy spread is huge and TEAPOT
# accuracy is not enough for tracking.
# This will slow down tracking and it is not symplectic.
#----------------------------------------------------------
#accLattice.setLinacTracker(True)


quads = accLattice.getQuads()


#-----Let's create bunc
eKin = 1.0

bunch_init = Bunch()
syncPart = bunch_init.getSyncParticle()
#set H- mass
#self.bunch_init.mass(0.9382723 + 2*0.000511)
bunch_init.mass(0.939294)
bunch_init.charge(-1.0)
syncPart.kinEnergy(eKin)

#---- add one particle to the bunch
x0  = 0.
xp0 = 0.
y0  = 0.
yp0 = 0.
z0  = 0.
dE  = 0.

bunch_init.addParticle(x0,xp0,y0,yp0,z0,dE)

bunch = Bunch()
bunch_init.copyBunchTo(bunch)

#---- cat off levels for Gaussian distribution in sigma values
cut_off_level = 3.0

#---- setup random offsets for quads 
offset_error = 0.0001 
coordDisplModification = CoordinateDisplacementNodesModification()
coordDisplModification.addLatticeNodes(quads)
coordDisplModification.setGaussDistributedDisplacementParameter("dx",offset_error,cut_off_level)
coordDisplModification.setGaussDistributedDisplacementParameter("dy",offset_error,cut_off_level)

#---- setup random roration angles
angle_rms = 0.0001
#---- setup random X-axis roration angles
angleModificationX = StraightRotationX_NodesModification()
angleModificationX.addLatticeNodes(quads)
angleModificationX.setGaussDistributedAngle(angle_rms,cut_off_level)
#---- setup random Y-axis roration angles
angleModificationY = StraightRotationY_NodesModification()
angleModificationY.addLatticeNodes(quads)
angleModificationY.setGaussDistributedAngle(angle_rms,cut_off_level)
#---- setup random Z-axis roration angles
angleModificationZ = StraightRotationZ_NodesModification()
angleModificationZ.addLatticeNodes(quads)
angleModificationZ.setGaussDistributedAngle(angle_rms,cut_off_level)

#---- setup random quads fields
relative_error = 0.001
quadFieldModification = QuadFieldsErrorsDeployment()
quadFieldModification.addQuads(quads)
quadFieldModification.setGaussDistributedRealtiveErrors(relative_error,cut_off_level)
	
#set up design
accLattice.trackDesignBunch(bunch)
accLattice.trackBunch(bunch)

trajCorrection = TrajectoryCorrection(accLattice)

print("===============Initial Trajectory  =============")
bunch = Bunch()
bunch_init.copyBunchTo(bunch)
trajCorrection.calculateTrajectory(bunch, print_info = True)
print("======================================")


#---- correct trajectory
bunch = Bunch()
bunch_init.copyBunchTo(bunch)
trajCorrection.correctTrajectory(bunch)

#---- print the DCH and DCV fields
print("Maximal achivable field in DCH(V) in RTBT1 is 0.017 [T] ") 
dcvs = trajCorrection.getDCVs()
for dc in dcvs:
	print("DCV = ",dc.getName()," B [T] = %+6.5f "%dc.getParam("B"))
print("==============")
dchs = trajCorrection.getDCHs()
for dc in dchs:
	print("DCH = ",dc.getName()," B [T] = %+6.5f "%dc.getParam("B"))

print("===============Corrected Trajectory=============")
bunch = Bunch()
bunch_init.copyBunchTo(bunch)
trajCorrection.calculateTrajectory(bunch, print_info = True)
print("======================================")

print("Done.")

