import math
import pytest
import os

from orbit.core.bunch import Bunch
from orbit.teapot import teapot
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.utils import consts
from orbit.impedances import addImpedanceNode
from orbit.impedances import BetFreqDep_LImpedance_Node


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

print("Start.")

# ------------------------------------------
# Make a lattice
# ------------------------------------------


def getLattice(lattice_length, n_parts):
    elem = teapot.DriftTEAPOT("a drift")
    elem.setLength(lattice_length)
    elem.setnParts(n_parts)
    teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
    teapot_lattice.addNode(elem)
    teapot_lattice.initialize()
    return teapot_lattice


lattice_length = 248.0  # the length of the lattice
# number of parts into which the drift will be chopped:
n_parts = 1
# also the number of SC nodes
lattice = getLattice(lattice_length, n_parts)

# ------------------------------
# Bunch initialization
# ------------------------------

b = Bunch()
print("Read Bunch.")
runName = "Benchmark_Collimator"

total_macroSize = 1.0e10
b.mass(consts.mass_proton)

ERef = 1.0  # Gev
energy = 1.5
print("energy is:", energy)
bunch_orbit_to_pyorbit(lattice.getLength(), energy, os.path.join(script_dir, "../SpaceCharge/sc1D/Bm_KV_Uniform_10000"), b)
b.getSyncParticle().kinEnergy(energy)
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(total_macroSize / nParticlesGlobal)

ERef = 1.0
GammaRef = 1.0 + ERef / b.mass()
BetaRef = math.sqrt(1.0 - 1.0 / (GammaRef * GammaRef))
FreqRef = (consts.speed_of_light * BetaRef) / lattice_length
print("BetaRef =", BetaRef)
print("FreqRef =", FreqRef)

Gamma = 1.0 + b.getSyncParticle().kinEnergy(energy) / b.mass()
Beta = math.sqrt(1.0 - 1.0 / (Gamma * Gamma))
Freq = (consts.speed_of_light * Beta) / lattice_length
print("Beta =", Beta)
print("Freq =", Freq)

# -------------------------------------------------------------------------
# Set up longitudinal impedance nodes:
# -------------------------------------------------------------------------

length = 248.0
nMacrosMin = 1
nBins = 64  # number of longitudinal slices
position = 0.0

# -------------------------------------------------------------------------
# Impedance data:
# SNS Longitudinal Impedance tables. EKicker impedance
# from private communication with J.G. Wang.
# Impedance in Ohms/n. Kicker and RF impedances are inductive with
# real part positive and imaginary part negative by Chao definition.
# -------------------------------------------------------------------------

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
for i in range(len(ZL_EKicker)):
    zk = ZL_EKicker[i]
    zrf = ZL_RF[i]
    # Multiply by 10 to make the effect bigger for benchmark purpose.
    zreal = 10.0 * (zk.real / 1.75 + zrf.real)
    zimag = 10.0 * (zk.imag / 1.75 + zrf.imag)
    Z.append(complex(zreal, zimag))

ZM = []
for i in range(len(ZL_EKicker) - 1):
    zk = ZL_EKicker[i]
    zrf = ZL_RF[i]
    # Multiply by 10 to make the effect bigger for benchmark purpose.
    zreal = 10.0 * (zk.real / 1.75 + zrf.real)
    zimag = 10.0 * (zk.imag / 1.75 + zrf.imag)
    ZM.append(complex(zreal, zimag))


# -------------------------------------------------------------------------
# Impedance dictionaries:
# f_impeDict is for frequency-dependent impedance, FreqDep_LImpedance_Node.
# compbf_impeDict is for comparison between BetFreqDep_LImpedance_Node
# and FreqDep_LImpedance_Node. Hence, impedance is independent of beta.
# In bf_impeDict, impedance does change with frequency for test of
# BetFreqDep_LImpedance_Node.
# -------------------------------------------------------------------------

f_list = []
for i in range(len(Z)):
    freq = i * FreqRef
    f_list.append(freq)

# -------------------------------------------------------------------------
# Additional data for beta dependent tests.
# -------------------------------------------------------------------------

# BetaRef = 0.875025655404
beta_list = []

nBetas = 10
BetaMin = 0.1
BetaMax = 1.0
for i in range(nBetas + 1):
    vb = BetaMin + i * ((BetaMax - BetaMin) / nBetas)
    beta_list.append(vb)

# -------------------------------------------------------------------------
# The impedence must be a double indexed list of lists.
# The first index corresponds to beta.
# The second the index corresponds the frequency.
# -------------------------------------------------------------------------

compbfz = []
bfz = []


# -------------------------------------------------------------------------
def lMult(n, list):
    """
    Returns a new list with each value multiplied by n.
    """
    newList = []
    for val in list:
        nv = val * n
        newList.append(nv)
    return newList


# I take the existing Z(freq) list and mulitply it by
# some factor for each beta value.
# -------------------------------------------------------------------------

for i in range(len(beta_list)):
    comp = lMult(1, Z)  # Use this one to compare fd and bfd nodes
    v = lMult((i / 10.0), Z)
    compbfz.append(comp)
    bfz.append(v)

bf_impeDict = {}
bf_impeDict["freqs"] = f_list
bf_impeDict["betas"] = beta_list
bf_impeDict["z_imp"] = bfz

# ------------------------------------------------------------------------
# Testing:
# Note: the f_ and compbf_ nodes should have the same Z,
# and bfd should have a lower Z.
# ------------------------------------------------------------------------

paramsDict = {}
paramsDict["bunch"] = b
b.dumpBunch("bunch_init.dat")


bf_impedancenode = BetFreqDep_LImpedance_Node(length, nMacrosMin, nBins, b, bf_impeDict)
addImpedanceNode(lattice, position, bf_impedancenode)


print("===========Lattice modified =======================================")
print("New Lattice = ", lattice.getName(), " length [m] = ", lattice.getLength(), " nodes = ", len(lattice.getNodes()))

print("Ready to track")

bf_impedancenode.trackBunch(b)

print("tracking done")

b.dumpBunch("bunch_final.dat")
bunch_final = read_values_from_file("bunch_final.dat")

bunch_pyorbit_to_orbit(lattice.getLength(), b, "pybunch_final.dat")
print("Stop.")


def test_final_bunch():
    expected_bunch_final = os.path.join(script_dir, "expected_bfd_LImp_final.dat")
    expected_bunch_final = read_values_from_file(expected_bunch_final)

    assert len(bunch_final) == len(expected_bunch_final)

    for a, e in zip(bunch_final, expected_bunch_final):
        assert a == pytest.approx(e, abs=0.000000001)


os.remove("bunch_final.dat")
os.remove("pybunch_final.dat")
os.remove("bunch_init.dat")
