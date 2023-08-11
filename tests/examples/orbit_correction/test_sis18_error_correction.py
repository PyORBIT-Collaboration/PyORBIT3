# ------------------------------------------------
# pyORBIT error and correction example
# ------------------------------------------------

# import orbit.core
import pytest
import os

from orbit.teapot import TEAPOT_Lattice
from bunch import Bunch
from orbit.errors import AddErrorSet


# ATTENTION !!! The python packet numpy and scipy are required
from orbit.orbit_correction import orbit, correction


def read_tuples_from_file(filename):
    result = []
    with open(filename, "r") as file:
        for line in file:
            values = line.strip().split()
            if values:
                tuple_values = tuple(map(float, values))
                result.append(tuple_values)
    return result


script_dir = os.path.dirname(__file__)

print("Start.")
# ---------------------------------------------Bunch init---------------------------------------------
b = Bunch()
b.mass(0.93827231)
b.macroSize(1.0)

energy = 1.0  # Gev
b.getSyncParticle().kinEnergy(energy)
# ---------------------------------------------Bunch init---------------------------------------------

print("Generate Lattice.")
# ---------------------------------------------Make a Teapot Lattice----------------------------------
lattice = TEAPOT_Lattice("lattice")
lattice.readMAD(os.path.join(script_dir, "sis18.lat"), "SIS18")
# ---------------------------------------------Make a Teapot Lattice----------------------------------

print("INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES")
# ---------------------------------------------ORBIT ERRORS-------------------------------------------
# WE INTRODUCE MISALIGNEMENT IN THE QUADRUPOLES; dx, dy = HOR AND VER DISPLACEMENT OF QUADRUPOLES
setDict = {}
paramsDict = {}
positioni = 0.0
positionf = lattice.getLength()
paramsDict["errtype"] = "StraightError"
paramsDict["subtype"] = "TransDisp"
paramsDict["sample"] = "Uniform"
paramsDict["maximum"] = 0.5
paramsDict["minimum"] = 0.0
paramsDict["dx"] = 0.005
paramsDict["dy"] = 0.007

setDict["elementtype"] = "quad"
setDict["ringline"] = "ring"

ESet = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict, seed_value=50)
# ESet  = AddErrorSet(lattice, positioni, positionf, setDict, paramsDict) # Random
# ---------------------------------------------ORBIT ERRORS-------------------------------------------


print("CALCULATE DISTORTED ORBIT AND PLOT")
# ---------------------------------------------CALCULATE DISTORTED ORBIT------------------------------
OrbitX, OrbitY = orbit(lattice, b).get_orbit()
# ---------------------------------------------CALCULATE DISTORTED ORBIT------------------------------

# ---------------------------------------------PLOT DISTORTED ORBIT-----------------------------------
x = []
y = []
s = []
for i in range(len(OrbitX)):
    s.append(OrbitX[i][0])
    x.append(OrbitX[i][1])
    y.append(OrbitY[i][1])

# ---------------------------------------------PLOT DISTORTED ORBIT-----------------------------------

print("CORRECTED ORBIT")
# ---------------------------------------------CORRECTED ORBIT----------------------------------------
corr = correction(lattice, b)
corr.orbit_corr()
# ---------------------------------------------CORRECTED ORBIT----------------------------------------

print("CALCULATE CORRECTED ORBIT AND PLOT")
# ---------------------------------------------CALCULATE CORRECTED ORBIT------------------------------
OrbitX_corr, OrbitY_corr = orbit(lattice, b).get_orbit()


def test_OrbitX():
    expected_OrbitX_corr = os.path.join(script_dir, "OrbitX_corrected.dat")
    expected_OrbitX_corr = read_tuples_from_file(expected_OrbitX_corr)

    assert len(OrbitX_corr) == len(expected_OrbitX_corr)

    for a, e in zip(OrbitX_corr, expected_OrbitX_corr):
        assert a == pytest.approx(e, abs=0.0000001)


def test_OrbitY():
    expected_OrbitY_corr = os.path.join(script_dir, "OrbitY_corrected.dat")
    expected_OrbitY_corr = read_tuples_from_file(expected_OrbitY_corr)

    assert len(OrbitY_corr) == len(expected_OrbitY_corr)

    for a, e in zip(OrbitY_corr, expected_OrbitY_corr):
        assert a == pytest.approx(e, abs=0.0000001)
