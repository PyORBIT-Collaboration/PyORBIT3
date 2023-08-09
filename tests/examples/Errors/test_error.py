# ------------------------------------------------
# pyORBIT error module benchmark
# ------------------------------------------------

import sys
import math
import pytest
import orbit.core
import os

import orbit_mpi

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch
from orbit.errors import AddErrorNode
from orbit.errors import AddErrorSet

from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

print("Start.")
# ------------------------------
# Bunch init
# ------------------------------
script_dir = os.path.dirname(__file__)
bunch_input_file = os.path.join(script_dir, "pyorbit_bunch_input.dat")
test_lattice = os.path.join(script_dir, "./LATTICES/Test.LAT")

b = Bunch()
print("Read Bunch.")
b.readBunch(bunch_input_file)
b.mass(0.93827231)
b.macroSize(1.0)

energy = 1.0  # Gev
b.getSyncParticle().kinEnergy(energy)

# ------------------------------
# Make a Teapot Lattice
# ------------------------------

print("Generate Lattice.")
lattice = TEAPOT_Lattice("no_sc_lattice")
lattice.readMAD(test_lattice, "TEST")

setDict = {}
paramsDict = {}

positioni = 2.8
positionf = 3.2
paramsDict["errtype"] = "RotationError"
paramsDict["elementtype"] = "straight"
paramsDict["subtype"] = "YS"
paramsDict["sample"] = "Fixed"
paramsDict["angle"] = 0.1

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict)

z = 0.0
string = "Start lattice, z =  {}\n".format(z)

# Set the number of sections in quads to the same as for ORBIT_MPI
for acc_elem in lattice.getNodes():
    z += acc_elem.getLength()
    string += "Node =  {}  type =  {}  L =  {}  N child nodes =  {}  z =  {}\n".format(
        acc_elem.getName(), acc_elem.getType(), acc_elem.getLength(), acc_elem.getNumberOfChildren(), z
    )
    if acc_elem.getType() == "quad teapot":
        acc_elem.setnParts(5)

print(string)
# dump initial bunch for ORBIT_MPI input
# bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch_input.dat")

# =====track bunch ============
ACC_TURNS = 1
print("Tracking.")
for i in range(ACC_TURNS):
    lattice.trackBunch(b)
    print("Turn ", i)

# dump ORBIT_MPI bunch to compare results
# bunch_pyorbit_to_orbit(lattice.getLength(), b, "bunch_output.dat")

print("STOP.")


def test_string():
    expected_string = """Start lattice, z =  0.0
Node =  QFH  type =  quad teapot  L =  0.25  N child nodes =  4  z =  0.25
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  0.5
Node =  S  type =  solenoid teapot  L =  0.25  N child nodes =  4  z =  0.75
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  1.0
Node =  B  type =  bend teapot  L =  1.0  N child nodes =  4  z =  2.0
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  2.25
Node =  K  type =  kick teapot  L =  0.25  N child nodes =  4  z =  2.5
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  2.75
Node =  QDH  type =  quad teapot  L =  0.25  N child nodes =  5  z =  3.0
Node =  QDH  type =  quad teapot  L =  0.25  N child nodes =  5  z =  3.25
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  3.5
Node =  M  type =  multipole teapot  L =  0.25  N child nodes =  4  z =  3.75
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  4.0
Node =  B  type =  bend teapot  L =  1.0  N child nodes =  4  z =  5.0
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  5.25
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  5.5
Node =  D  type =  drift teapot  L =  0.25  N child nodes =  4  z =  5.75
Node =  QFH  type =  quad teapot  L =  0.25  N child nodes =  4  z =  6.0\n"""
    assert string == expected_string
