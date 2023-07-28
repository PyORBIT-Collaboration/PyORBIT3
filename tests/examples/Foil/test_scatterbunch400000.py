import orbit.core
import pytest

import os
import math
import sys
from bunch import Bunch
from foil import Foil
from orbit.injection import InjectParts
import random
from orbit.core.orbit_utils import random as orbit_random

random.seed(100)
orbit_random.seed(0)


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


print("Start.")
myfile = "out12.txt"
myfile2 = "out13.txt"

xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050
# Below is 1000 times the width of normal foil but will do only one turn.
thick = 400000

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
main_bunch_400000 = read_lines(myfile)
print(main_bunch_400000)
print("=========LOST BUNCH==========")
lostbunch.dumpBunch(myfile2)
lost_bunch_400000 = read_lines(myfile2)
print(lost_bunch_400000)


def test_main_bunch_400000_thickness():
    # expected bunch is kept in a seperate file since it is too large
    expected_main_bunch_400000 = os.path.join(script_dir, "scatteredbunch400000.dat")
    expected_main_bunch_400000 = read_lines(expected_main_bunch_400000)

    assert main_bunch_400000 == expected_main_bunch_400000


def test_lost_bunch_400000_thickness():
    assert (
        lost_bunch_400000
        == """-0.010324893 0.0028276738 0.0052921438 -0.0017178431 107.75573 -0.00064100469 0
-0.00082089906 -0.0031856719 -0.0022896579 -0.00061952887 68.840279 -0.00071778829 0
-0.0087950427 -0.0036073099 -0.0017151565 -0.00017626514 -96.623375 -0.0006486367 0
-0.010176136 -2.6248166e-06 -0.0026141972 -0.0038619162 57.901529 -0.00078837163 0
-0.011389086 -0.0026861197 0.010329863 -5.7315716e-06 -102.19636 -0.00013807739 0
-0.017664141 0.00040303335 -0.0087506216 0.0042985987 -33.425301 -0.00068363343 0"""
    )


os.remove(myfile)
os.remove(myfile2)
