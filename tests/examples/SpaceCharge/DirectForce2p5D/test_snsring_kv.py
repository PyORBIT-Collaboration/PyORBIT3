##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting
# injection nodes
##############################################################

import math
import os
import pytest

from orbit.teapot import teapot
from orbit.core.bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit
from orbit.space_charge.directforce2p5d import directforceLatticeModifications
from orbit.core.spacecharge import SpaceChargeForceCalc2p5D
from orbit.space_charge.sc1d import SC1D_AccNode
from orbit.rf_cavities import RFNode

print("Start.")


def read_values_from_file(file_path):
    values = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line:
                values.extend(map(float, line.split()))

    return values


# =====Main bunch parameters============
script_dir = os.path.dirname(__file__)
intensity = 1.0e14
turns = 1000.0

b = Bunch()
b.mass(0.93827231)

energy = 1.0  # Gev
b.readBunch(os.path.join(script_dir, "KV.dat"), 10)
b.getSyncParticle().kinEnergy(energy)
macrosize = intensity / b.getSize()
b.macroSize(macrosize)

print(macrosize)

paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"] = lostbunch
paramsDict["bunch"] = b
lostbunch.addPartAttr("LostParticleAttributes")

# =====Make a Teapot style lattice======

teapot_latt = teapot.TEAPOT_Ring()
print("Read MAD.")
teapot_latt.readMAD(os.path.join(script_dir, "SNSring_pyOrbitBenchmark.LAT"), "RING")
print("Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))
lattlength = teapot_latt.getLength()

# ----------------------------------------------
# Add one black absorber collimator to act like
# an aperture
# ----------------------------------------------
colllength = 0.00001
ma = 9
density_fac = 1.0
shape = 1
radius = 0.110

# -----------------------------
# Add RF Node
# -----------------------------

teapot_latt.initialize()
ZtoPhi = 2.0 * math.pi / lattlength
dESync = 0.0
RF1HNum = 1.0
RF1Voltage = 0.000016
RF1Phase = 0.0
RF2HNum = 2.0
RF2Voltage = -0.000003
RF2Phase = 0.0
length = 0.0

rf1_node = RFNode.Harmonic_RFNode(ZtoPhi, dESync, RF1HNum, RF1Voltage, RF1Phase, length, "RF1")
rf2_node = RFNode.Harmonic_RFNode(ZtoPhi, dESync, RF2HNum, RF2Voltage, RF2Phase, length, "RF2")
position1 = 196.0
position2 = 196.5
# RFLatticeModifications.addRFNode(teapot_latt, position1, rf1_node)
# RFLatticeModifications.addRFNode(teapot_latt, position2, rf2_node)

# ----------------------------------------------
# make 2.5D space charge calculator
# ----------------------------------------------

sizeX = 32  # number of grid points in horizontal direction
sizeY = 32  # number of grid points in vertical direction
sizeZ = 1  # number of longitudinal slices in the 2.5D space charge solver
forcecalc2p5d = SpaceChargeForceCalc2p5D(sizeX, sizeY, sizeZ)
sc_path_length_min = 0.00000001
directforceLatticeModifications.setDirectForce2p5DAccNodes(teapot_latt, sc_path_length_min, forcecalc2p5d)

# -----------------------------------------------
# Add longitudinal space charge node with Imped
# -----------------------------------------------

b_a = 10.0 / 3.0
length = lattlength
nMacrosMin = 1000
useSpaceCharge = 1
nBins = 128  # number of longitudinal slices in the 1D space charge solver
position = 64.0

# SNS Longitudinal Impedance tables. EKicker impedance from private communication
# with J.G. Wang. Seems to be for 7 of the 14 kickers (not sure why).
# Impedance in Ohms/n. Kicker and RF impedances are inductive with real part positive and imaginary is negative by Chao definition.

ZL_EKicker = [
    complex(42.0, -182),
    complex(35, -101.5),
    complex(30.3333, -74.6667),
    complex(31.5, -66.5),
    complex(32.2, -57.4),
    complex(31.5, -51.333),
    complex(31, -49),
    complex(31.5, -46.375),
    complex(31.8889, -43.556),
    complex(32.9, -40.6),
    complex(32.7273, -38.18),
    complex(32.25, -35.58),
    complex(34.46, -32.846),
    complex(35, -30.5),
    complex(35.4667, -28.0),
    complex(36.75, -25.81),
    complex(36.647, -23.88),
    complex(36.944, -21.1667),
    complex(36.474, -20.263),
    complex(36.4, -18.55),
    complex(35.333, -17),
    complex(35, -14.95),
    complex(33.478, -13.69),
    complex(32.375, -11.67),
    complex(30.8, -10.08),
    complex(29.615, -8.077),
    complex(28.519, -6.74),
    complex(27.5, -5),
    complex(26.552, -4.103),
    complex(25.433, -3.266),
    complex(24.3871, -2.7),
    complex(23.40625, -2.18),
]

ZL_RF = [
    complex(0.0, 0.0),
    complex(0.750, 0.0),
    complex(0.333, 0.0),
    complex(0.250, 0.0),
    complex(0.200, 0.0),
    complex(0.167, 0.0),
    complex(3.214, 0.0),
    complex(0.188, 0.0),
    complex(0.167, 0.0),
    complex(0.150, 0.0),
    complex(1.000, 0.0),
    complex(0.125, 0.0),
    complex(0.115, 0.0),
    complex(0.143, 0.0),
    complex(0.333, 0.0),
    complex(0.313, 0.0),
    complex(0.294, 0.0),
    complex(0.278, 0.0),
    complex(0.263, 0.0),
    complex(0.250, 0.0),
    complex(0.714, 0.0),
    complex(0.682, 0.0),
    complex(0.652, 0.0),
    complex(0.625, 0.0),
    complex(0.600, 0.0),
    complex(0.577, 0.0),
    complex(0.536, 0.0),
    complex(0.536, 0.0),
    complex(0.517, 0.0),
    complex(0.500, 0.0),
    complex(0.484, 0.0),
    complex(0.469, 0.0),
]


Z = []
for i in range(0, 32):
    zk = ZL_EKicker[i]
    zrf = ZL_RF[i]
    zreal = zk.real / 1.75 + zrf.real
    zimag = zk.imag / 1.75 + zrf.imag
    Z.append(complex(zreal, zimag))

sc1Dnode = SC1D_AccNode(b_a, length, nMacrosMin, useSpaceCharge, nBins)
# sc1Dnode.assignImpedance(Z);
# addLongitudinalSpaceChargeNode(teapot_latt, position, sc1Dnode)

# -------------------------------
#  Lattice is ready
# -------------------------------

nodes = teapot_latt.getNodes()
i = 0
for node in nodes:
    print(i, " node=", node.getName(), " s start,stop = %4.3f %4.3f " % teapot_latt.getNodePositionsDict()[node])
    print("There are ", node.getNumberOfBodyChildren(), " child nodes.")
    i = i + 1

# ================Do some turns===========================================

# newlatt = teapot_latt.getSubLattice(0, 693)
# newlatt.trackBunch(b, paramsDict)
teapot_latt.trackBunch(b, paramsDict)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "mainbunch.dat")

main_bunch = read_values_from_file("mainbunch.dat")


def test_ring_main_bunch():
    expected_main_bunch = read_values_from_file(os.path.join(script_dir, "expectedmainbunch.dat"))
    assert len(main_bunch) == len(expected_main_bunch)

    for e, a in zip(expected_main_bunch, main_bunch):
        assert e == pytest.approx(a, abs=0.000000001)


print("Stop.")

os.remove("mainbunch.dat")
