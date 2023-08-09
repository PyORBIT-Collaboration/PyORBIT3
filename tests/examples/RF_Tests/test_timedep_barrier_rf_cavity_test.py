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


print("Start.")

temp_bunch = "temp_bunch.txt"

b = Bunch()
b.addParticle(0.0, 0.0, 0.0, 0.0, -1.8, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, -1.5, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, -1.2, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, -0.9, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, -0.6, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, -0.3, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 0.3, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 0.6, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 0.9, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 1.2, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 1.5, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 1.8, 0.0)
b.compress()

syncPart = b.getSyncParticle()
energy = 100.0  # energy in GeV
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
accelDict["n_tuple"] = 8
accelDict["time"] = (0, 5.0e-09, 10.0e-09, 15.0e-09, 20.0e-09, 25.0e-09, 30.0e-09, 35.0e-09, 40.0e-09)
accelDict["RFVoltage"] = (0.1, 0.11, 0.12, 0.13, 0.14, 0.13, 0.12, 0.11, 0.1)
accelDict["RFPhasep"] = (150.0, 140.0, 130.0, 120.0, 110.0, 120.0, 130.0, 140.0, 150.0)
accelDict["RFPhasem"] = (-150.0, -140.0, -130.0, -120.0, -110.0, -120.0, -130.0, -140.0, -150.0)
accelDict["dRFPhasep"] = (30.0, 40.0, 50.0, 60.0, 70.0, 60.0, 50.0, 40.0, 30.0)
accelDict["dRFPhasem"] = (30.0, 40.0, 50.0, 60.0, 70.0, 60.0, 50.0, 40.0, 30.0)

length = 0.0
name = "barrier_rfnode"
rf_node = RFNode.TimeDep_Barrier_RFNode(ZtoPhi, accelDict, b, length, name)
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


def test_timedep_rf_cavity():
    expected_bunch_after = """0 0 0 0 -1.7999979 0.2771211
0 0 0 0 -1.4999967 0.40895008
0 0 0 0 -1.1999981 0.22871667
0 0 0 0 -0.89999918 0.088031113
0 0 0 0 -0.59999998 0.0023661764
0 0 0 0 -0.3 0
0 0 0 0 0 0
0 0 0 0 0.3 0
0 0 0 0 0.59999998 -0.0023661764
0 0 0 0 0.89999918 -0.088031113
0 0 0 0 1.1999981 -0.22871667
0 0 0 0 1.4999967 -0.40895007
0 0 0 0 1.7999979 -0.2771211"""
    assert bunch_after == expected_bunch_after


os.remove(temp_bunch)
