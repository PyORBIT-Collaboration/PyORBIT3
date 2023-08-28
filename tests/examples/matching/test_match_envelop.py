import os
import math
import pytest

from orbit.teapot import TEAPOT_Lattice

from orbit.core.bunch import Bunch

from orbit.matching import Optics, EnvelopeSolver

script_dir = os.path.dirname(__file__)


# ------------------------------------------------------
# Read MADX twiss file, match and plot envelopes
# -------------------------------------------------------

bunch = Bunch()

energy = 11.4e-3  # beta must be on def. of disperion
bunch.getSyncParticle().kinEnergy(energy)

lattice = TEAPOT_Lattice("no_sc_lattice")
lattice.readMAD(os.path.join(script_dir, "fodo.lat"), "CELLA")  # lattice start at injection point
lattice.setUseRealCharge(useCharge=1)

beamline = Optics().readtwiss_teapot(lattice, bunch)

beamline.print_line()

emitx = 12.5e-6  # emittance_x
emity = 12.5e-6  # emittance_y
sigma_p = 1.0e-3  # rms momentum spread
Ksc = 1.0e-6 * 5  # space charge perveance
circum = 216.0  # 1080.0 # circumference
Ncells = 12.0  # number of cells

solve = EnvelopeSolver(beamline)


envx0, envxs0, envy0, envys0, Dx0, Dxs0, s = solve.match_root(emitx, emity, 0.0, 0.0)  # matched envelopes (no space charge)
phasex0, phasey0 = solve.phase_advance(envx0, envy0, Dx0, emitx, emity, 0.0, s)  # phase advance (no space charge)
print("zero current phase advance x [deg]:", phasex0 * 180.0 / math.pi)
print("zero current phase advance y [deg]:", phasey0 * 180.0 / math.pi)
envx, envxs, envy, envys, Dx, Dxs, s = solve.match_root(emitx, emity, sigma_p, Ksc)  # matched envelopes (with space charge)
phasex, phasey = solve.phase_advance(envx, envy, Dx, emitx, emity, sigma_p, s)  # phase advance (with space charge)

print("finite current phase advance x [deg]:", phasex * 180.0 / math.pi)
print("finite current phase advance y [deg]:", phasey * 180.0 / math.pi)


# values taken from original PyOrbit examples
expected_phasex = 0.4781374662241269
expected_phasey = 0.47845404552827536


def test_phaseX_phaseY():
    assert phasex == pytest.approx(expected_phasex, abs=0.0000001)
    assert phasey == pytest.approx(expected_phasey, abs=0.0000001)
