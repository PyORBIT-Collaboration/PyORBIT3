##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting
# injection nodes
##############################################################

import pytest
import os
import random

from orbit.teapot import teapot
from orbit.core.bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit
from orbit.injection import TeapotInjectionNode
from orbit.injection import addTeapotInjectionNode
from orbit.injection import JohoTransverse, UniformLongDist


def read_values_from_file(file_path):
    values = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line:
                values.extend(map(float, line.split()))

    return values


random.seed(100)
script_dir = os.path.dirname(__file__)

print("Start.")

teapot_latt = teapot.TEAPOT_Lattice()
print("Read MAD.")
teapot_latt.readMAD(os.path.join(script_dir, "../Apertures/MAD_Lattice/LATTICE"), "RING")
print("Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

# ====Injection aperature============
xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050

injectparams = (xmin, xmax, ymin, ymax)

# =====set up bunch stuff============
b = Bunch()

b.mass(0.93827231)
b.macroSize(1.0e1)
energy = 1.0  # Gev
b.getSyncParticle().kinEnergy(energy)

paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"] = lostbunch
paramsDict["bunch"] = b
lostbunch.addPartAttr("LostParticleAttributes")

# ------------------------------
# Initial Distribution Functions
# ------------------------------
sp = b.getSyncParticle()

order = 3.0
alphax = 0.063
betax = 10.209
alphay = 0.063
betay = 10.776
emitlim = 0.00152 * 2 * (order + 1) * 1e-6
xcenterpos = 0.0468
xcentermom = 0.001
ycenterpos = 0.0492
ycentermom = -0.00006
tailfrac = 0.1  # 10% in tails
tailfac = 1.5  # tail emittance is 50% greater than core
zlim = 120.0 * 248.0 / 360.0
zmin = -zlim
zmax = zlim
deltaEfrac = 0.001 / 2.0
eoffset = 0.1
sp = b.getSyncParticle()

xFunc = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom, tailfrac, tailfac)
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom, tailfrac, tailfac)
lFunc = UniformLongDist(zmin, zmax, sp, eoffset, deltaEfrac)

nparts = 10000.0

injectnode = TeapotInjectionNode(nparts, b, lostbunch, injectparams, xFunc, yFunc, lFunc)

addTeapotInjectionNode(teapot_latt, 0.0, injectnode)

print("===========Lattice modified =======================================")
print("New Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

injectnode.track(paramsDict)

# dump ORBIT_MPI bunch to compare results
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), lostbunch, "lostbunch.dat")
print("Stop.")

mainbunch = read_values_from_file("mainbunch.dat")


def test_injection_main_bunch():
    expected_mainbunch = os.path.join(script_dir, "expectedmainbunch.dat")
    expected_mainbunch = read_values_from_file(expected_mainbunch)

    assert len(mainbunch) == len(expected_mainbunch)

    for e, a in zip(expected_mainbunch, mainbunch):
        assert e == pytest.approx(a, abs=0.0000000001)


os.remove("mainbunch.dat")
os.remove("lostbunch.dat")
