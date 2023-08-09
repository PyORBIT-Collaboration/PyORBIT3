from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.core.bunch import Bunch
from orbit.rf_cavities import RFNode, RFLatticeModifications

import math
import os
import pytest


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


temp_bunch = "temp_bunch.txt"

print("Start.")

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
accelDict = {}
accelDict["gammaTrans"] = 1.0e10
accelDict["RFHNum"] = 1
accelDict["n_tuple"] = 8
accelDict["time"] = (0, 5.0e-09, 10.0e-09, 15.0e-09, 20.0e-09, 25.0e-09, 30.0e-09, 35.0e-09, 40.0e-09)
accelDict["SyncPhase"] = (0.0, 10.0, 20.0, 30.0, 40.0, 30.0, 20.0, 10.0, 0.0)
accelDict["RFVoltage"] = (0.1, 0.11, 0.12, 0.13, 0.14, 0.13, 0.12, 0.11, 0.1)
accelDict["RFPhase"] = (0.0, 30.0, 60.0, 90.0, 90.0, 90.0, 60.0, 30.0, 0.0)
length = 0.0
name = "harmonic_rfnode"
rf_node = RFNode.SyncPhaseDep_Harmonic_RFNode(ZtoPhi, accelDict, b, length, name)
position = 1.0
RFLatticeModifications.addRFNode(lattice, position, rf_node)

print("Lattice length = ", lattice.getLength())

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
lattice.trackActions(accContainer, paramsDict)
lattice.trackActions(accContainer, paramsDict)
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


def test_syncphasedep_rf_cavity():
    expected_bunch_after = """0.001 0 0 0 -0.00078986123 -0.0013118323
0.01447853 0.00085864263 0 0 -0.00079565342 -0.0013104034
0 0 0.001 0 -0.00078986123 -0.0013118323
0 0 0.01447853 0.00085864263 -0.00079565342 -0.0013104034
0 0 0 0 0.30176901 -0.45642588
0 0 0 0 0.00090121999 -0.00073868826"""
    assert bunch_after == expected_bunch_after


os.remove(temp_bunch)
