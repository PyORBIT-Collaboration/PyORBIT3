import math
import sys
import random
import sys
import orbit.core

from bunch import Bunch
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from spacecharge import Boundary2D
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D
from orbit.diagnostics import TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotDiagnosticsNode

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import GaussDist3D

print("Start.")

# =====Make a Teapot style lattice======

lattice = teapot.TEAPOT_Lattice()
print("Read MAD.")
lattice.readMAD("SNSring_pyOrbitBenchmark.LAT", "RING")
print("Lattice=", lattice.getName(), " length [m] =", lattice.getLength(), " nodes=", len(lattice.getNodes()))

# ------------------------------
# Main Bunch init
# ------------------------------

n_particles = 100000
twissX = TwissContainer(alpha=0.0046902, beta=10.207, emittance=3.0e-5)
twissY = TwissContainer(alpha=0.056823, beta=10.639, emittance=3.0e-5)
twissZ = TwissContainer(alpha=0.0, beta=100000.0, emittance=0.008)
dist = GaussDist3D(twissX, twissY, twissZ)

b = Bunch()
for i in range(n_particles):
    (x, xp, y, yp, z, zp) = dist.getCoordinates()
    b.addParticle(x, xp, y, yp, z, zp)

total_macroSize = 1.0e14
b.mass(0.93827231)
energy = 1.0  # Gev
# b.readBunch(distribution_file, n_particles)
print("Bunch Generated.")
b.getSyncParticle().kinEnergy(energy)
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(total_macroSize / nParticlesGlobal)

# -----------------------------------
# Add Tune Analysis node
# -----------------------------------

tunes = TeapotTuneAnalysisNode("tune_analysis")
tunes.assignTwiss(10.207, 0.0469, -0.05, 0.0061, 10.639, 0.056)
addTeapotDiagnosticsNode(lattice, 0, tunes)

# -----------------------------------
# Add Space Charge nodes
# -----------------------------------

nMacrosMin = 1
sc_path_length_min = 0.00000001
sizeX = 32  # number of grid points in horizontal direction
sizeY = 32  # number of grid points in vertical direction
sizeZ = 16  # number of longitudinal slices

calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)
scLatticeModifications.setSC2p5DAccNodes(lattice, sc_path_length_min, calc2p5d)

# -----------------------------------
# Tracking
# -----------------------------------

paramsDict = {}
paramsDict["bunch"] = b

n_turns = 2
for i in range(n_turns):
    lattice.trackBunch(b, paramsDict)

b.dumpBunch("bunch_2p5D_turn" + str(n_turns) + ".dat")

print("Stop.")
