from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch

# ///////////////////////////////////////////////////////////
from orbit.rf_cavities import RFNode, RFLatticeModifications

# ///////////////////////////////////////////////////////////

import math
import orbit.core

print("Start.")


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


myfile = "out3.txt"

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
RFVoltage = 0.1
RFPhasep = 150.0
RFPhasem = -150.0
dRFPhasep = 30.0
dRFPhasem = 30.0
length = 0.0
name = "barrier_rfnode"
rf_node = RFNode.Barrier_RFNode(ZtoPhi, RFVoltage, RFPhasep, RFPhasem, dRFPhasep, dRFPhasem, length, name)
position = 1.0
RFLatticeModifications.addRFNode(lattice, position, rf_node)

print("Lattice length = ", lattice.getLength())
print("ZtoPhi = ", ZtoPhi)
print("RFVoltage = ", RFVoltage)
print("RFPhasep  = ", RFPhasep)
print("RFPhasem  = ", RFPhasem)
print("dRFPhasep = ", dRFPhasep)
print("dRFPhasem = ", dRFPhasem)

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
b.dumpBunch(myfile)
mystring = read_lines(myfile)
print(mystring)
print("==========================================")

print("lattice length=", lattice.getLength())
print("beta=", b.getSyncParticle().beta())
print("TEAPOT time[sec]=", b.getSyncParticle().time())
print("SIMPLE time[sec]=", lattice.getLength() / (b.getSyncParticle().beta() * 2.99792458e8))
print("Stop.")


def test_barrier_rf_cavity_bunch_after():
    expected = """0 0 0 0 -1.7645982 0.080901699
0 0 0 0 -1.4687555 0.070710678
0 0 0 0 -1.2 0
0 0 0 0 -0.9 0
0 0 0 0 -0.6 0
0 0 0 0 -0.3 0
0 0 0 0 0 0
0 0 0 0 0.3 0
0 0 0 0 0.6 0
0 0 0 0 0.9 0
0 0 0 0 1.2 0
0 0 0 0 1.4639496 -0.070710678
0 0 0 0 1.7582998 -0.080901699"""
    assert mystring == expected
