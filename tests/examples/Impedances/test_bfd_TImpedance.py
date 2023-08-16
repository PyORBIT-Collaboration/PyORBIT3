import math
import os
import pytest

from orbit.core.bunch import Bunch
from orbit.teapot import teapot
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.utils import consts
from orbit.impedances import addImpedanceNode
from orbit.impedances import BetFreqDep_TImpedance_Node


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

total_macroSize = 1.0e16
b.mass(consts.mass_proton)

ERef = 1.0  # Gev
print("Reference energy is 1.0 (GeV).")
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
# Set up transverse impedance nodes:
# -------------------------------------------------------------------------

length = 248.0
nMacrosMin = 1
nBins = 64  # number of longitudinal slices
position = 0.0

# -------------------------------------------------------------------------
# Impedance data:
# SNS transverse impedance tables. EKicker impedance
# from measurements by H. Hahn and D. Davino.
# Impedance in Ohms/m. Kicker and RF impedances are inductive with
# real part positive and imaginary part negative by Chao definition.
# -------------------------------------------------------------------------

ZP = [
    complex(36904.387, -92036.662),
    complex(36456.161, -87439.318),
    complex(34716.164, -84162.832),
    complex(32207.922, -81707.918),
    complex(29454.361, -79576.055),
    complex(26851.512, -77430.568),
    complex(24512.290, -75295.880),
    complex(22511.300, -73245.276),
    complex(20923.146, -71352.039),
    complex(19822.135, -69689.310),
    complex(19219.099, -68299.947),
    complex(18983.254, -67159.242),
    complex(18964.653, -66233.345),
    complex(19013.351, -65488.407),
    complex(18979.725, -64890.578),
    complex(18782.876, -64405.997),
    complex(18495.238, -64000.784),
    complex(18209.990, -63641.055),
    complex(18020.314, -63292.924),
    complex(18019.158, -62922.586),
    complex(18250.869, -62512.992),
    complex(18651.347, -62084.468),
    complex(19141.817, -61662.404),
    complex(19643.506, -61272.186),
    complex(20077.786, -60939.140),
    complex(20396.751, -60675.860),
    complex(20621.040, -60466.526),
    complex(20780.566, -60291.473),
    complex(20905.245, -60131.038),
    complex(21024.960, -59965.589),
    complex(21163.150, -59782.544),
    complex(21328.879, -59585.044),
    complex(21529.263, -59378.360),
    complex(21771.418, -59167.759),
    complex(22062.438, -58958.503),
    complex(22404.471, -58753.963),
    complex(22788.628, -58553.296),
    complex(23204.528, -58355.089),
    complex(23641.788, -58157.927),
    complex(24090.024, -57960.399),
    complex(24537.942, -57761.603),
    complex(24972.225, -57561.770),
    complex(25379.283, -57361.285),
    complex(25745.526, -57160.534),
    complex(26057.403, -56959.901),
    complex(26309.939, -56759.625),
    complex(26517.289, -56559.624),
    complex(26696.196, -56359.770),
    complex(26863.405, -56159.937),
    complex(27035.628, -55959.998),
    complex(27223.327, -55759.896),
    complex(27423.005, -55559.734),
    complex(27629.277, -55359.634),
    complex(27836.760, -55159.718),
    complex(28040.083, -54960.108),
    complex(28236.753, -54760.789),
    complex(28430.693, -54561.439),
    complex(28626.696, -54361.694),
    complex(28829.555, -54161.192),
    complex(29044.039, -53959.570),
    complex(29269.662, -53756.948),
    complex(29494.225, -53554.511),
    complex(29703.939, -53353.589),
    complex(29885.018, -53155.515),
]

ZM = [
    complex(-36904.387, -92036.662),
    complex(-35537.319, -98454.148),
    complex(-31917.808, -107107.170),
    complex(-26247.086, -117791.080),
    complex(-19012.507, -130023.330),
    complex(-10702.778, -143320.070),
    complex(-1806.603, -157197.430),
    complex(7187.314, -148828.430),
    complex(15790.266, -135241.380),
    complex(23513.549, -122525.270),
    complex(29868.459, -111163.960),
    complex(34366.289, -101641.310),
    complex(36604.713, -94357.279),
    complex(36819.781, -89091.792),
    complex(35533.669, -85346.877),
    complex(33269.902, -82623.249),
    complex(30552.005, -80421.622),
    complex(27865.194, -78291.571),
    complex(25411.569, -76143.769),
    complex(23266.333, -74050.737),
    complex(21504.092, -72085.759),
    complex(20199.451, -70322.119),
    complex(19407.853, -68823.957),
    complex(19043.126, -67587.851),
    complex(18955.621, -66580.093),
    complex(18995.393, -65766.834),
    complex(19012.495, -65114.223),
    complex(18877.728, -64588.408),
    complex(18615.224, -64155.514),
    complex(18317.839, -63781.657),
    complex(18078.752, -63432.953),
    complex(17991.144, -63075.517),
    complex(18133.522, -62680.523),
    complex(18475.945, -62256.524),
    complex(18939.871, -61828.830),
    complex(19446.526, -61422.827),
    complex(19917.136, -61063.902),
    complex(20282.202, -60773.597),
    complex(20540.771, -60545.044),
    complex(20722.612, -60358.639),
    complex(20857.641, -60194.716),
    complex(20975.770, -60033.613),
    complex(21104.968, -59857.792),
    complex(21258.828, -59665.442),
    complex(21444.495, -59461.79),
    complex(21669.088, -59252.132),
    complex(21939.722, -59041.711),
    complex(22262.022, -58835.235),
    complex(22630.575, -58633.189),
    complex(23035.023, -58434.167),
    complex(23464.985, -58236.757),
    complex(23910.079, -58039.543),
    complex(24359.649, -57841.267),
    complex(24801.017, -57641.803),
    complex(25220.597, -57441.533),
    complex(25604.796, -57240.842),
    complex(25940.026, -57040.116),
    complex(26215.284, -56839.695),
    complex(26438.699, -56639.600),
    complex(26626.975, -56439.702),
    complex(26796.853, -56239.876),
    complex(26965.080, -56039.994),
    complex(27146.509, -55839.951),
    complex(27342.041, -55639.799),
]

ZF = [
    complex(0.000, -157197.430),
    complex(7187.314, -148828.430),
    complex(15790.266, -135241.380),
    complex(23513.549, -122525.270),
    complex(29868.459, -111163.960),
    complex(34366.289, -101641.310),
    complex(36604.713, -94357.279),
    complex(36819.781, -89091.792),
    complex(35533.669, -85346.877),
    complex(33269.902, -82623.249),
    complex(30552.005, -80421.622),
    complex(27865.194, -78291.571),
    complex(25411.569, -76143.769),
    complex(23266.333, -74050.737),
    complex(21504.092, -72085.759),
    complex(20199.451, -70322.119),
    complex(19407.853, -68823.957),
    complex(19043.126, -67587.851),
    complex(18955.621, -66580.093),
    complex(18995.393, -65766.834),
    complex(19012.495, -65114.223),
    complex(18877.728, -64588.408),
    complex(18615.224, -64155.514),
    complex(18317.839, -63781.657),
    complex(18078.752, -63432.953),
    complex(17991.144, -63075.517),
    complex(18133.522, -62680.523),
    complex(18475.945, -62256.524),
    complex(18939.871, -61828.830),
    complex(19446.526, -61422.827),
    complex(19917.136, -61063.902),
    complex(20282.202, -60773.597),
    complex(20540.771, -60545.044),
    complex(20722.612, -60358.639),
    complex(20857.641, -60194.716),
    complex(20975.770, -60033.613),
    complex(21104.968, -59857.792),
    complex(21258.828, -59665.442),
    complex(21444.495, -59461.79),
    complex(21669.088, -59252.132),
    complex(21939.722, -59041.711),
    complex(22262.022, -58835.235),
    complex(22630.575, -58633.189),
    complex(23035.023, -58434.167),
    complex(23464.985, -58236.757),
    complex(23910.079, -58039.543),
    complex(24359.649, -57841.267),
    complex(24801.017, -57641.803),
    complex(25220.597, -57441.533),
    complex(25604.796, -57240.842),
    complex(25940.026, -57040.116),
    complex(26215.284, -56839.695),
    complex(26438.699, -56639.600),
    complex(26626.975, -56439.702),
    complex(26796.853, -56239.876),
    complex(26965.080, -56039.994),
    complex(27146.509, -55839.951),
    complex(27342.041, -55639.799),
    complex(27542.041, -55439.799),
    complex(27742.041, -55239.799),
    complex(27942.041, -55039.799),
    complex(28142.041, -54839.799),
    complex(28342.041, -54639.799),
    complex(28542.041, -54439.799),
]

for i in range(len(ZP)):
    # Multiply by 100 to make the effect bigger for benchmark purpose.
    ZP[i] = 100.0 * ZP[i]

for i in range(len(ZM)):
    # Multiply by 100 to make the effect bigger for benchmark purpose.
    ZM[i] = 100.0 * ZM[i]

for i in range(len(ZF)):
    # Multiply by 100 to make the effect bigger for benchmark purpose.
    ZF[i] = 100.0 * ZF[i]

# -------------------------------------------------------------------------
# Impedance dictionaries:
# f_impeDict is for frequency-dependent impedance, FreqDep_TImpedance_Node.
# compbf_impeDict is for comparison between BetFreqDep_TImpedance_Node
# and FreqDep_TImpedance_Node. Hence, impedance is independent of beta.
# In bf_impeDict, impedance does change with beta for test of
# BetFreqDep_TImpedance_Node.
# -------------------------------------------------------------------------

f_list = []
for i in range(len(ZP)):
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
for i in range(nBetas):
    vb = BetaMin + i * ((BetaMax - BetaMin) / (nBetas - 1))
    beta_list.append(vb)

# -------------------------------------------------------------------------
# The impedence must be a double indexed list of lists.
# The first index corresponds to beta.
# The second index corresponds the frequency.
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


# Take the existing Z(freq) list and mulitply it by
# some factor for each beta value.
# -------------------------------------------------------------------------

for i in range(len(beta_list)):
    comp = lMult(1, ZF)  # Use this one to compare fd and bfd nodes
    v = lMult(((i + 1) / 10.0), ZF)
    compbfz.append(comp)
    bfz.append(v)

f_impeDict = {}
f_impeDict["freqs"] = f_list
f_impeDict["zx_imp"] = ZF
f_impeDict["zy_imp"] = ZF

compbf_impeDict = {}
compbf_impeDict["freqs"] = f_list
compbf_impeDict["betas"] = beta_list
compbf_impeDict["zx_imp"] = compbfz
compbf_impeDict["zy_imp"] = compbfz

bf_impeDict = {}
bf_impeDict["freqs"] = f_list
bf_impeDict["betas"] = beta_list
bf_impeDict["zx_imp"] = bfz
bf_impeDict["zy_imp"] = bfz

# ------------------------------------------------------------------------
# Testing:
# Note: the f_ and compbf_ nodes should have the same Z,
# and bfd should have a lower Z.
# ------------------------------------------------------------------------

paramsDict = {}
paramsDict["bunch"] = b
b.dumpBunch("bunch_init.dat")

qX = 6.2
alphaX = 0.0
betaX = 10.0
qY = 6.2
alphaY = 0.0
betaY = 10.0

useX = 1
useY = 1
bf_impedancenode = BetFreqDep_TImpedance_Node(length, nMacrosMin, nBins, useX, useY, b, bf_impeDict, qX, alphaX, betaX, qY, alphaY, betaY)
addImpedanceNode(lattice, position, bf_impedancenode)


print("===========Lattice modified =======================================")
print("New Lattice = ", lattice.getName(), " length [m] = ", lattice.getLength(), " nodes = ", len(lattice.getNodes()))

print("Ready to track")

# impedancenode.trackBunch(b)
# f_impedancenode.trackBunch(b)
# compbf_impedancenode.trackBunch(b)
bf_impedancenode.trackBunch(b)

print("tracking done")

b.dumpBunch("bunch_final.dat")
bunch_final = read_values_from_file("bunch_final.dat")

bunch_pyorbit_to_orbit(lattice.getLength(), b, "pybunch_final.dat")
print("Stop.")


def test_final_bunch():
    expected_bunch_final = os.path.join(script_dir, "expected_bfd_TImp_final.dat")
    expected_bunch_final = read_values_from_file(expected_bunch_final)

    assert len(bunch_final) == len(expected_bunch_final)

    for a, e in zip(bunch_final, expected_bunch_final):
        assert a == pytest.approx(e, abs=0.000000001)


os.remove("bunch_init.dat")
os.remove("bunch_final.dat")
os.remove("pybunch_final.dat")
