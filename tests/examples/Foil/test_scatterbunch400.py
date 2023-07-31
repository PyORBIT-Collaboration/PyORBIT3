import pytest

import os
import math
import sys
from orbit.core.bunch import Bunch
from orbit.core.foil import Foil
from orbit.injection import InjectParts
import random
from orbit.core.orbit_utils import random as orbit_random


# Setting the random seed for the Python random module and for the orbit_utils random number generator
random.seed(100)
orbit_random.seed(0)


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


print("Start.")
myfile = "out10.txt"
myfile2 = "out11.txt"

xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050
# Below is 1000 times the width of normal foil but will do only one turn.
thick = 400

foil = Foil(xmin, xmax, ymin, ymax, thick)
script_dir = os.path.dirname(__file__)
parts_file = os.path.join(script_dir, "parts.dat")

# ------------------------------
# Main Bunch init
# ------------------------------
b = Bunch()
print("Read Bunch.")
# runName = "Benchmark_Collimator"

b.mass(0.93827231)
b.macroSize(1.0e1)
energy = 1.0  # Gev
b.readBunch(parts_file)
b.getSyncParticle().kinEnergy(energy)

# =====track bunch through Foil============

lostbunch = Bunch()
lostbunch.addPartAttr("LostParticleAttributes")

foil.traverseFoilFullScatter(b, lostbunch)
b.dumpBunch(myfile)
main_bunch_400 = read_lines(myfile)
print(main_bunch_400)
print("=========LOST BUNCH==========")
lostbunch.dumpBunch(myfile2)
lost_bunch_400 = read_lines(myfile2)
print(lost_bunch_400)


def test_main_bunch_400_thickness():
    # expected bunch is kept in a seperate file since it is too large
    expected_main_bunch_400 = os.path.join(script_dir, "scatteredbunch400.dat")
    expected_main_bunch_400 = read_lines(expected_main_bunch_400)

    assert main_bunch_400 == expected_main_bunch_400


def test_lost_bunch_400_thickness():
    assert lost_bunch_400 == ""


os.remove(myfile)
os.remove(myfile2)
