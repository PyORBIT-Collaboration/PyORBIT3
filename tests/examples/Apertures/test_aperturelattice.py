##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting
# diagnotics nodes
##############################################################

import math
import sys
import orbit.core
import pytest
import os

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit


from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject
from orbit.aperture import addTeapotApertureNode
from orbit.aperture import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from orbit.aperture import addCircleApertureSet, addEllipseApertureSet, addRectangleApertureSet

print("Start.")


def read_values_from_file(file_path):
    values = []

    with open(file_path) as f:
        for line in f:
            if not line.startswith("%"):
                line = line.strip()
                if line:
                    values.extend(map(float, line.split()))

    return values


script_dir = os.path.dirname(__file__)

teapot_latt = teapot.TEAPOT_Lattice()
print("Read MAD.")
teapot_latt.readMAD(os.path.join(script_dir, "MAD_Lattice/LATTICE"), "RING")
print("Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

Aperturenode = CircleApertureNode(0.01)
# Aperturenode = CircleApertureNode(.01, .002, .0025)
# Aperturenode = RectangleApertureNode(.02, .015, -.002, -.0025)
addTeapotApertureNode(teapot_latt, 240, Aperturenode)
# addEllipseApertureSet(.02, .01,  teapot_latt, 150, 195)

print("===========Lattice modified =======================================")
print("New Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

print("============= nodes inside the region ===========")
# print all nodes around the specified position
for node in teapot_latt.getNodes():
    print("node=", node.getName(), " type=", node.getType(), " L=", node.getLength())

# ------------------------------
# Main Bunch init
# ------------------------------

b = Bunch()
print("Read Bunch.")
runName = "Benchmark_Diagnostics"
b.mass(0.93827231)
b.macroSize(1.0e1)
energy = 1.0  # Gev
bunch_orbit_to_pyorbit(teapot_latt.getLength(), energy, os.path.join(script_dir, "Bm_KV_Uniform_1000"), b)
b.getSyncParticle().kinEnergy(energy)
paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"] = lostbunch
paramsDict["bunch"] = b
lostbunch.addPartAttr("LostParticleAttributes")

# =====track bunch ============

print("Tracking...")
for turn in range(1):
    teapot_latt.trackBunch(b, paramsDict)
b.dumpBunch("bunch_final.dat")
lostbunch.dumpBunch("lostbunch_final.dat")

bunch_final = read_values_from_file("bunch_final.dat")

print("Stop.")


def test_final_bunch():
    expected_bunch_final = os.path.join(script_dir, "expected_bunch_final.dat")
    expected_bunch_final = read_values_from_file(expected_bunch_final)

    assert len(bunch_final) == len(expected_bunch_final)

    for e, a in zip(expected_bunch_final, bunch_final):
        assert e == pytest.approx(a, abs=0.00000000001)


os.remove("bunch_final.dat")
os.remove("lostbunch_final.dat")
