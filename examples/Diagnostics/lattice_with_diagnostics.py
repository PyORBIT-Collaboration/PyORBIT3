##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting
# diagnotics nodes
##############################################################


from orbit.teapot import teapot
from orbit.core.bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit
from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotStatLatsNodeSet, addTeapotMomentsNodeSet
from orbit.utils.orbit_mpi_utils import bunch_pyorbit_to_orbit

print("Start.")

# =====Main bunch parameters============

macrosize = 1

teapot_latt = teapot.TEAPOT_Ring()
print("Read MAD.")
teapot_latt.readMAD("SNSring_pyOrbitBenchmark.LAT", "RING")
print("Lattice=", teapot_latt.getName(), " length [m] =", teapot_latt.getLength(), " nodes=", len(teapot_latt.getNodes()))

b = Bunch()
b.mass(0.93827231)
b.macroSize(macrosize)
energy = 1.0  # Gev
b.getSyncParticle().kinEnergy(energy)
b.readBunch("Bunches/controlbunch_600.dat", 1000)
bunch_pyorbit_to_orbit(teapot_latt.getLength(), b, "Bunches/controlbunch_600_ORBIT.dat")
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
