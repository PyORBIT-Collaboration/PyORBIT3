##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting
# diagnotics nodes
##############################################################

import pytest
import os


from orbit.teapot import teapot
from orbit.core.bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit
from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotStatLatsNodeSet, addTeapotMomentsNodeSet
from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit


def read_lines(file):
    with open(file, "r") as f:
        lines = f.readlines()

    stripped_line = [line.strip() for line in lines if not line.startswith("%")]
    stripped_content = "\n".join(stripped_line)

    return stripped_content


def read_values_from_file(file_path):
    values = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line:
                values.extend(map(float, line.split()))

    return values


script_dir = os.path.dirname(__file__)
SNS_ring_file = os.path.join(script_dir, "SNSring_pyOrbitBenchmark.LAT")
print("Start.")

# =====Main bunch parameters============

macrosize = 1
teapot_latt = teapot.TEAPOT_Ring()
print("Read MAD.")
teapot_latt.readMAD(SNS_ring_file, "RING")
print("Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))


controlbunch_600 = os.path.join(script_dir, "Bunches/controlbunch_600.dat")
controlbunch_600_ORBIT = os.path.join(script_dir, "Bunches/controlbunch_600_ORBIT.dat")

b = Bunch()
b.mass(0.93827231)
b.macroSize(macrosize)
energy = 1.0  # Gev
b.getSyncParticle().kinEnergy(energy)
b.readBunch(controlbunch_600, 1000)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, controlbunch_600_ORBIT)
# bunch_orbit_to_pyorbit(teapot_latt.getLength(), energy, "Bunches/Bm_KV_Uniform_1000",b)

tunes = TeapotTuneAnalysisNode("tune_analysis")
# tunes.assignTwiss(3.25011, 0.811013, 1.45074, -0.636723, 10.6922, -2.13656)
# addTeapotDiagnosticsNode(teapot_latt, 38.788, tunes)

tunes.assignTwiss(9.335, -1.826, -0.065, -0.03, 8.089, 0.497)
addTeapotDiagnosticsNode(teapot_latt, 51.035, tunes)


addTeapotMomentsNodeSet(teapot_latt, "moments", 3)
addTeapotStatLatsNodeSet(teapot_latt, "statlats")

print("===========Lattice modified =======================================")
print("New Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

# print "============= nodes inside the region ==========="
# print all nodes around the specified position
for node in teapot_latt.getNodes():
    print("node=", node.getName(), " type=", node.getType(), " L=", node.getLength())

# =====track bunch ============
teapot_latt.trackBunch(b)
b.dumpBunch("final1.dat")
teapot_latt.trackBunch(b)
b.dumpBunch("final2.dat")
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "bunch2.dat")


print("Stop.")

final_1 = read_lines("final1.dat")
final_2 = read_lines("final2.dat")
actual_bunch_2 = read_values_from_file("bunch2.dat")
actual_statlats = read_values_from_file("statlats")


def test_diagnostics_final_1():
    expected_final_1 = os.path.join(script_dir, "expectedfinal1.dat")
    expected_final_1 = read_lines(expected_final_1)

    assert final_1 == expected_final_1


def test_diagnostics_final_2():
    expected_final_2 = os.path.join(script_dir, "expectedfinal2.dat")
    expected_final_2 = read_lines(expected_final_2)

    assert final_2 == expected_final_2


def test_pyorbit_to_orbit_bunch2():
    expected_bunch_2 = os.path.join(script_dir, "expectedbunch2.dat")
    expected_bunch_2 = read_values_from_file(expected_bunch_2)

    assert len(actual_bunch_2) == len(expected_bunch_2)

    for e, a in zip(expected_bunch_2, actual_bunch_2):
        assert e == pytest.approx(a, abs=0.0000000001)


os.remove("moments")
os.remove("statlats")
os.remove("final1.dat")
os.remove("final2.dat")
os.remove("bunch2.dat")
