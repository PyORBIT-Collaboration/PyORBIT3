##############################################################
# This script will test the functionality of the EnergyAperture
# class of the "aperture" module
##############################################################

import math
import sys
import pytest
import os

from orbit.core.bunch import Bunch
from orbit.core.aperture import EnergyAperture


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


myfile = "out6.txt"
myfile2 = "out7.txt"

print("Start.")

# ------------------------------
# Main Bunch init
# ------------------------------

b = Bunch()
b.mass(0.93827231)
b.macroSize(1.0e1)
energy = 1.0  # Gev
b.getSyncParticle().kinEnergy(energy)


dE_spread = 0.010

nParts = 10
for ind in range(nParts):
    x = 0.0
    xp = 0.0
    y = 0.0
    yp = 0.0
    z = 0.0
    dE = (2 * dE_spread / nParts) * (ind - nParts / 2)
    b.addParticle(x, xp, y, yp, z, dE)

lostbunch = Bunch()

# ==== make EnergyAperture class instance

energyAperture = EnergyAperture()

# ==== check set get parameters methods
energyAperture.setPosition(111.0)
print("energyAperture            pos =", energyAperture.getPosition())
energyAperture.setMinMaxEnergy(-0.005, +0.005)
print("energyAperture min max energy =", energyAperture.getMinMaxEnergy())

# =====track bunch through the EnergyAperture ============
print("Tracking...")

# ---- this  will collect lost particles in the lost bunch
energyAperture.checkBunch(b, lostbunch)

# ---- if you do not care about the lost particles you can do this:
# energyAperture.checkBunch(b)

print("==============init bunch          ==========")
b.dumpBunch(myfile)
init_bunch = read_lines(myfile)
print(init_bunch)
print("==============lost particles bunch==========")
lostbunch.dumpBunch(myfile2)
lost_particle_bunch = read_lines(myfile2)
print(lost_particle_bunch)

print("Stop.")


# test init bunch
def test_phase_aperture_smth():
    expected_init_bunch = "0 0 0 0 0 -0.004\n" "0 0 0 0 0 -0.002\n" "0 0 0 0 0 0\n" "0 0 0 0 0 0.002\n" "0 0 0 0 0 0.004"
    assert init_bunch == expected_init_bunch


def test_lost_particles_bunch():
    expected_lost_particles_bunch = (
        "0 0 0 0 0 -0.01 111\n" "0 0 0 0 0 -0.008 111\n" "0 0 0 0 0 -0.006 111\n" "0 0 0 0 0 0.006 111\n" "0 0 0 0 0 0.008 111"
    )
    assert lost_particle_bunch == expected_lost_particles_bunch


os.remove(myfile)
os.remove(myfile2)
