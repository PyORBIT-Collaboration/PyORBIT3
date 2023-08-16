from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.core.bunch import Bunch
from orbit.rf_cavities import RFNode, RFLatticeModifications

import orbit.core
from rfcavities import Harmonic_Cav
from rfcavities import Dual_Harmonic_Cav

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
RFHNum = 1.0
RatioRFHNum = 2.0
RFVoltage = 0.1
RatioVoltage = 0.5
RFPhase = 1.0
RFPhase2 = 1.0
length = 0.0
name = "dual_harmonic_rfnode"
position = 0.0

rf_node = RFNode.Dual_Harmonic_RFNode(ZtoPhi, RFHNum, RatioRFHNum, RFVoltage, RatioVoltage, RFPhase, RFPhase2, length, name)
RFLatticeModifications.addRFNode(lattice, position, rf_node)


print("Lattice length = ", lattice.getLength())
print("ZtoPhi = ", ZtoPhi)
print("RFHNum = ", RFHNum)
print("RFVoltage = ", RFVoltage)
print("RatioVoltage = ", RatioVoltage)
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

for i in range(10):
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


def test_dual_rf_cavity():
    expected_bunch_after = """0.001 0 0 0 0 0
0.04 0.001 0 0 -1.9995023e-05 2.1480161e-09
0 0 0.001 0 0 0
0 0 0.04 0.001 -1.9995023e-05 2.1480161e-09
0 0 0 0 -0.98893646 -0.15147618
0 0 0 0 0.0063068599 0.000998734"""
    assert bunch_after == expected_bunch_after


os.remove(temp_bunch)
