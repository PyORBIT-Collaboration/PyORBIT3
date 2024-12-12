#! /usr/bin/env python

"""
This script is the linac solenoid test.
We will compare results with TEAPOT solenoid and will test the lattice parser. 
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory

from orbit.core.bunch import Bunch

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.lattice import Quad,Solenoid

from orbit.teapot import SolenoidTEAPOT

#------------------------------------------------
#    Start of Script
#------------------------------------------------

def printBunch(text,bunch):
	"""
	Print particles coordinates in the bunch.
	"""
	print (text)
	b = bunch
	for ind in range(bunch.getSize()):
		(x, xp, y, yp, z, dE) = (b.x(ind),b.xp(ind),b.y(ind),b.yp(ind),b.z(ind),b.dE(ind))
		st  = " %2d "%ind 
		st += " %12.5g , %12.5g "%(x*1000., xp*1000.)
		st += " %12.5g , %12.5g "%(y*1000., yp*1000.)
		st += " %12.5g , %12.5g "%(z, dE*1000.)
		print (st)
	print ("==================================")

names = ["TEST_Seq",]

sns_linac_factory = SNS_LinacLatticeFactory()
sns_linac_factory.setMaxDriftLength(0.5)

# ---- the XML file name with the structure
xml_file_name = "test_solenoid_lattice.xml"

# ---- make lattice from XML file
accLattice = sns_linac_factory.getLinacAccLattice(names, xml_file_name)

print ("============ Test Lattice =============")
for ind,node in enumerate(accLattice.getNodes()):
	st  = "ind= %2d "%ind
	st += " node= %35s "%node.getName()
	st += " length[m] = %6.4f "%node.getLength()
	st += " nParts= %2d "%node.getnParts()
	st += " type= %10s "%node.getType()
	print (st)
print ("========================================")


node_soln = accLattice.getNodesOfClass(Solenoid)[0]
node_soln_ind = accLattice.getNodeIndex(node_soln)
print ("Solenoid field B=",node_soln.getParam("B"))

#---- Results should be the same for any parts in the solenoid node
nSolParts = 10
node_soln.setnParts(nSolParts)
print (" nParts= %2d "%node_soln.getnParts())

bunch_init = Bunch()
syncPart = bunch_init.getSyncParticle()
bunch_init.mass(0.9382723)
bunch_init.charge(+1.0)
syncPart.kinEnergy(0.003)

(x, xp, y, yp, z, dE) = (0.,0.,0.,0.,0.,0.)
bunch_init.addParticle(x, xp, y, yp, z, dE)

(x, xp, y, yp, z, dE) = (0.001,0.001,-0.001,-0.001,0.1,0.00002)
bunch_init.addParticle(x, xp, y, yp, z, dE)
bunch_init.addParticle(x/2, xp/2, y/2, yp/2, z/2, dE/2)

(x, xp, y, yp, z, dE) = (-0.001,-0.001,+0.001,+0.001,-0.1,-0.00002)
bunch_init.addParticle(x, xp, y, yp, z, dE)
bunch_init.addParticle(x/2, xp/2, y/2, yp/2, z/2, dE/2)

#---- set up design
accLattice.trackDesignBunch(bunch_init)

#---- track bunch to the solenoid start
accLattice.trackBunch(bunch_init,None,None,0,node_soln_ind-1)

printBunch("========= Initial Bunch: before Solenoid ",bunch_init)

bunch = Bunch()
bunch_init.copyBunchTo(bunch)

#---- track bunch through the solenoid
accLattice.trackBunch(bunch,None,None,node_soln_ind,node_soln_ind)

printBunch("========= Bunch After: after Solenoid ",bunch)

#----------------------------------------------
# Now let's create TEAPOT solenoid and track
# the same bunch through it.
#----------------------------------------------

soln_teapot = SolenoidTEAPOT("solenoid_teapot")
soln_teapot.setParam("B",node_soln.getParam("B"))
soln_teapot.setLength(node_soln.getLength())

bunch = Bunch()
bunch_init.copyBunchTo(bunch)
soln_teapot.trackBunch(bunch)

printBunch("========= Bunch After: after TEAPOT Solenoid - should be the same",bunch)

print ("Stop.")


