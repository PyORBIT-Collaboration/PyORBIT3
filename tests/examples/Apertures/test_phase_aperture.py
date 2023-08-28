##############################################################
# This script will test the functionality of the PhaseAperture
# class of the "aperture" module
##############################################################

import os

from orbit.core.bunch import Bunch
from orbit.core.aperture import PhaseAperture


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


myfile = "out4.txt"
myfile2 = "out5.txt"
print("Start.")

# ------------------------------
# Main Bunch init
# ------------------------------

rf_frequency = 402.5e6
c = 2.99792458e8

b = Bunch()
b.mass(0.93827231)
b.macroSize(1.0e1)
energy = 1.0  # Gev
b.getSyncParticle().kinEnergy(energy)
beta = b.getSyncParticle().beta()

lambda_rf = c / rf_frequency
phase_to_z_coeff = lambda_rf * beta / 360.0

nParts = 10
for ind in range(nParts):
    x = 0.0
    xp = 0.0
    y = 0.0
    yp = 0.0
    z = (360.0 / nParts) * (ind - nParts / 2) * phase_to_z_coeff
    dE = 0.0
    b.addParticle(x, xp, y, yp, z, dE)

lostbunch = Bunch()

# ==== make PhaseAperture class instance

phaseAperture = PhaseAperture(2 * rf_frequency)

# ==== check set get parameters methods
print("phaseAperture  1   frequency =", phaseAperture.getRfFrequency())
phaseAperture.setRfFrequency(rf_frequency)
print("phaseAperture  2   frequency =", phaseAperture.getRfFrequency())
phaseAperture.setPosition(111.0)
print("phaseAperture            pos =", phaseAperture.getPosition())
phaseAperture.setMinMaxPhase(-100.0, +100.0)
print("phaseAperture min max phases =", phaseAperture.getMinMaxPhase())

# =====track bunch through the PhaseAperture ============
print("Tracking...")

# ---- this  will collect lost particles in the lost bunch
phaseAperture.checkBunch(b, lostbunch)

# ---- if you do not care about the lost particles you can do this:
# phaseAperture.checkBunch(b)

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
    expected_init_bunch = (
        "0 0 0 0 -0.13034836 0\n" "0 0 0 0 -0.065174181 0\n" "0 0 0 0 0 0\n" "0 0 0 0 0.065174181 0\n" "0 0 0 0 0.13034836 0"
    )
    assert init_bunch == expected_init_bunch


def test_lost_particles_bunch():
    expected_lost_particles_bunch = (
        "0 0 0 0 -0.32587091 0 111\n"
        "0 0 0 0 -0.26069673 0 111\n"
        "0 0 0 0 -0.19552254 0 111\n"
        "0 0 0 0 0.19552254 0 111\n"
        "0 0 0 0 0.26069673 0 111"
    )
    assert lost_particle_bunch == expected_lost_particles_bunch


os.remove(myfile)
os.remove(myfile2)
