##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting
# several collimation nodes
##############################################################

import pytest
import os
import random

from orbit.core.orbit_utils import random as orbit_random
from orbit.teapot import teapot
from orbit.teapot import DriftTEAPOT
from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit

from orbit.collimation import TeapotCollimatorNode
from orbit.collimation import addTeapotCollimatorNode


def read_values_from_file(file_path):
    values = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line:
                values.extend(map(float, line.split()))

    return values


orbit_random.seed(0)
random.seed(100)
script_dir = os.path.dirname(__file__)

print("Start.")

teapot_latt = teapot.TEAPOT_Lattice()
print("Read MAD.")
teapot_latt.readMAD(os.path.join(script_dir, "../Apertures/MAD_Lattice/LATTICE"), "RING")
print("Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

n_drifts = 0
drift_length = 0.0
for node in teapot_latt.getNodes():
    if isinstance(node, DriftTEAPOT):
        n_drifts += 1
        drift_length += node.getLength()

print("number of drifts =", n_drifts)
print("total drift length =", drift_length)
print("============= nodes inside the region ===========")
pos_start = 16.5
pos_stop = 20.5
# print all nodes around the specified position
for node in teapot_latt.getNodes():
    (node_pos_start, node_pos_stop) = teapot_latt.getNodePositionsDict()[node]
    if node_pos_start > pos_start and node_pos_stop < pos_stop:
        print("node=", node.getName(), " type=", node.getType(), "  pos=", node_pos_start, " L=", node.getLength())


length = 0.5
ma = 3
density_fac = 1.0
shape = 1
a = 0.0
b = 0
c = 0
d = 0
angle = 0
pos = 18.5

collimator = TeapotCollimatorNode(length, ma, density_fac, shape, a, b, c, d, angle)
# collimator = TeapotCollimatorNode(length, ma, density_fac, shape, a, b, c, d, angle)
addTeapotCollimatorNode(teapot_latt, pos, collimator)

print("===========Lattice modified =======================================")
print("New Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

n_drifts = 0
drift_length = 0.0
for node in teapot_latt.getNodes():
    if isinstance(node, DriftTEAPOT):
        n_drifts += 1
        drift_length += node.getLength()

print("number of drifts =", n_drifts)
print("total drift length =", drift_length)

print("============= nodes inside the region ===========")
# print all nodes around the specified position
for node in teapot_latt.getNodes():
    (node_pos_start, node_pos_stop) = teapot_latt.getNodePositionsDict()[node]
    if node_pos_start > pos_start and node_pos_stop < pos_stop:
        print("node=", node.getName(), " type=", node.getType(), "  pos=", node_pos_start, " L=", node.getLength())


# ------------------------------
# Main Bunch init
# ------------------------------
b = Bunch()
print("Read Bunch.")
runName = "Benchmark_Collimator"

b.mass(0.93827231)
b.macroSize(1.0e1)
energy = 1.0  # Gev
# get initial bunch from ORBIT_MPI input
bunch_orbit_to_pyorbit(teapot_latt.getLength(), energy, os.path.join(script_dir, "Bm_KV_Uniform_100"), b)
# b.readBunch("parts_in.dat")

b.getSyncParticle().kinEnergy(energy)

# =====track bunch through Collimator Node============
paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"] = lostbunch
paramsDict["bunch"] = b

lostbunch.addPartAttr("LostParticleAttributes")

# collimator.trackBunch(b,paramsDict)
collimator.track(paramsDict)

# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), lostbunch, "lostbunch.dat")
print("Stop.")

mainbunch = read_values_from_file("mainbunch.dat")


def test_collimator_lattice_mainbuch():
    expected_mainbunch = os.path.join(script_dir, "expectedmainbunch.dat")
    expected_mainbunch = read_values_from_file(expected_mainbunch)

    assert len(mainbunch) == len(expected_mainbunch)

    for e, a in zip(expected_mainbunch, mainbunch):
        assert e == pytest.approx(a, abs=0.000000001)


os.remove("mainbunch.dat")
os.remove("lostbunch.dat")
