from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.core.bunch import Bunch
from orbit.rf_cavities import RFNode, RFLatticeModifications

import math
import os
import pytest

print("Start.")


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


print("Start.")
temp_bunch = "temp_bunch.txt"
b = Bunch()
b.addParticle(1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0)
b.addParticle(0.0, 1.0e-3, 0.0, 0.0, 0.0, 0.0)
b.addParticle(0.0, 0.0, 1.0e-3, 0.0, 0.0, 0.0)
b.addParticle(0.0, 0.0, 0.0, 1.0e-3, 0.0, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3)
b.compress()

syncPart = b.getSyncParticle()
energy = 1.0  # energy in GeV
# p = syncPart.energyToMomentum(energy)
# syncPart.pz(p)
syncPart.kinEnergy(energy)

lattice = AccLattice("test_lattice")

elem0 = teapot.DriftTEAPOT("drift0")

lattice.addNode(elem0)

# -----------------------------
# Set TEAPOT nodes parameters
# -----------------------------

elem0.setLength(4.0)

lattice.initialize()

# ///////////////////////////////////////////////////////////
ZtoPhi = 2.0 * math.pi / lattice.getLength()
dESync = 0.0
RFHNum = 1.0
RFVoltage = 0.1
RFPhase = 0.0
length = 0.0
name = "harmonic_rfnode"
rf_node = RFNode.Harmonic_RFNode(ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase, length, name)
position = 1.0
RFLatticeModifications.addRFNode(lattice, position, rf_node)

print("Lattice length = ", lattice.getLength())
print("ZtoPhi = ", ZtoPhi)
print("dESync = ", dESync)
print("RFHNum = ", RFHNum)
print("RFVoltage = ", RFVoltage)
print("RFPhase = ", RFPhase)
# ///////////////////////////////////////////////////////////

print("==============BEFORE============================")
b.dumpBunch()
print("==========================================")


# =====track action ============
def bodyAction(paramsDict):
    node = paramsDict["node"]
    node.track(paramsDict)


accContainer = AccActionsContainer()
accContainer.addAction(bodyAction, AccActionsContainer.BODY)

paramsDict = {}
paramsDict["bunch"] = b

lattice.trackActions(accContainer, paramsDict)
print("=============AFTER=============================")
b.dumpBunch(temp_bunch)
bunch_after = read_lines(temp_bunch)
print(bunch_after)
print("==========================================")

print("lattice length=", lattice.getLength())
print("beta=", b.getSyncParticle().beta())
print("TEAPOT time[sec]=", b.getSyncParticle().time())
print("SIMPLE time[sec]=", lattice.getLength() / (b.getSyncParticle().beta() * 2.99792458e8))
print("Stop.")


def test_simple_rf_cavity():
    expected_bunch_after = """0.001 0 0 0 0 0
0.0039999998 0.001 0 0 -1.9627964e-06 7.8539816e-08
0 0 0.001 0 0 0
0 0 0.0039999998 0.001 -1.9627964e-06 7.8539816e-08
0 0 0 0 0.94737387 -0.1
0 0 0 0 0.00061923385 0.00097522276"""
    assert bunch_after == expected_bunch_after


os.remove(temp_bunch)
