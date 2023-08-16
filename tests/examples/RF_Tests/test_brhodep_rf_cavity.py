import math
import pytest
import os
from orbit.teapot import teapot
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.core.bunch import Bunch
from orbit.rf_cavities import RFNode, RFLatticeModifications


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
accelDict["BRho"] = (
    5.6573729517,
    5.7573729517,
    5.8573729517,
    5.9573729517,
    6.0573729517,
    5.9573729517,
    5.8573729517,
    5.7573729517,
    5.6573729517,
)
accelDict["RFVoltage"] = (0.1, 0.11, 0.12, 0.13, 0.14, 0.13, 0.12, 0.11, 0.1)
accelDict["RFPhase"] = (0.0, 30.0, 60.0, 90.0, 90.0, 90.0, 60.0, 30.0, 0.0)
length = 0.0
name = "harmonic_rfnode"
rf_node = RFNode.BRhoDep_Harmonic_RFNode(ZtoPhi, accelDict, b, length, name)
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

# ==============BEFORE============================
# ==========================================
# =============AFTER=============================
# % PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES
# % BUNCH_ATTRIBUTE_DOUBLE charge   1
# % BUNCH_ATTRIBUTE_DOUBLE classical_radius   1.5347e-18
# % BUNCH_ATTRIBUTE_DOUBLE macro_size   0
# % BUNCH_ATTRIBUTE_DOUBLE mass   0.938272
# %  SYNC_PART_COORDS 0 0 0  x, y, z positions in [m]
# %  SYNC_PART_MOMENTUM 0 0 1.696037912  px, py, pz momentum component in GeV/c
# %  info only: energy of the synchronous particle [GeV] = 1
# %  info only: momentum of the synchronous particle [GeV/c] = 1.696037912
# %  info only: beta=v/c of the synchronous particle = 0.8750256155
# %  info only: gamma=1/sqrt(1-(v/c)**2) of the synchronous particle = 2.065788684
# %  SYNC_PART_TIME 1.029253363e-08  time in [sec]
# % x[m] px[rad] y[m] py[rad] z[m]  (pz or dE [GeV])
# 0.00082842613 -0.00015325645 0 0 0.00029901071 -1.647656e-16
# 0.0023631085 0.00076945561 0 0 0.00043014195 -3.2953121e-16
# 2.3859816e-08 1.3257047e-07 0.0010236879 4.1172715e-05 4.4922855e-08 0
# -4.7883174e-07 6.8646524e-07 0.0028841539 0.0010928603 1.5869592e-06 -3.2953121e-16
# 0 0 0 0 1 0
# 0.00023639406 0.00019942888 0 0 -0.00041525426 0.001
# ==========================================
# lattice length= 2.7
# beta= 0.875025615499
# TEAPOT time[sec]= 1.02925336251e-08
# SIMPLE time[sec]= 1.02925336251e-08
# Stop.


def test_brhodep_rf_cavity():
    expected_bunch_after = """0.001 0 0 0 0.23665159 0.2136866
0.014534863 0.001 0 0 0.23664579 0.21368806
0 0 0.001 0 0.23665159 0.2136866
0 0 0.014534863 0.001 0.23664579 0.21368806
0 0 0 0 0.52876637 -0.2351167
0 0 0 0 0.23832426 0.21425499"""
    assert bunch_after == expected_bunch_after


os.remove(temp_bunch)
