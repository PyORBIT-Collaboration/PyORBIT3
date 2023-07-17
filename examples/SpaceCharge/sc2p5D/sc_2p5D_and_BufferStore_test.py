# -----------------------------------------------------
# This example tracks the particles through a drift.
# The space charge is calculated by using 2.5D SC node.
# -----------------------------------------------------
import sys
import math
import random

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.teapot import teapot
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalc2p5D, Boundary2D

from orbit.core import orbit_mpi

# ------------------------------
# Let's make everything random
# ------------------------------
random.seed(10)

print("Start.")
# ---------------------------------
# make an initial bunch
# ---------------------------------
b_init = Bunch()
bunch_radius = 0.05  # m
bunch_length = 100.0  # m
nParts = 500
total_macroSize = 1.4e14
energy = 1.0  # GeV

# ------------------------------------------------
# generate an uniform cylinder distribution
# no energy or momentum spread
# ------------------------------------------------
for ip in range(nParts):
    r = bunch_radius * math.sqrt(random.random())
    phi = 2 * math.pi * random.random()
    x = r * math.sin(phi)
    y = r * math.cos(phi)
    z = bunch_length * 0.5 * (1.0 - 2 * random.random())
    b_init.addParticle(x, 0.0, y, 0.0, z, 0.0)

# -------------------------------
# set bunch parameters
# -------------------------------
nParticlesGlobal = b_init.getSizeGlobal()
b_init.macroSize(total_macroSize / nParticlesGlobal)
syncPart = b_init.getSyncParticle()
syncPart.kinEnergy(energy)

# ------------------------------------------
#          Initial Bunch is ready
#
#          Now let's make a lattice
# ------------------------------------------


def getLattice(lattice_length, n_parts, spaceChargeCalculator2p5D):
    elem = teapot.DriftTEAPOT("my_drift")
    elem.setLength(lattice_length)
    elem.setnParts(n_parts)
    teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
    teapot_lattice.addNode(elem)
    teapot_lattice.initialize()
    # we will put SC nodes as frequently as possible.
    # In this case one for each part of the Drift node
    sc_path_length_min = 0.00000001
    scNodes_arr = scLatticeModifications.setSC2p5DAccNodes(teapot_lattice, sc_path_length_min, spaceChargeCalculator2p5D)
    return (teapot_lattice, scNodes_arr)


# ----------------------------------------------
# make 2.5D space charge calculator
# ----------------------------------------------
sizeX = 256  # number of grid points in horizontal direction
sizeY = 256  # number of grid points in vertical direction
sizeZ = 512  # number of longitudinal slices in the 2.5D space charge solver
calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)

lattice_length = 2.0  # the length of the drift
n_parts = 10  # number of parts on what the drift will be chopped, or the number of SC nodes
(lattice, scNodes_arr) = getLattice(lattice_length, n_parts, calc2p5d)

print("Number of SC nodes =", len(scNodes_arr))

# -------------------------------
#  Lattice is ready
# -------------------------------

nodes = lattice.getNodes()
for node in nodes:
    print("node=", node.getName(), " s start,stop = %4.3f %4.3f " % lattice.getNodePositionsDict()[node])
    print("There are ", node.getNumberOfBodyChildren(), " child nodes.")


print("=======================================================")

rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)

count = 0
while 1 < 2:
    # --------------------------------------------------------------------
    # Create a new bunch, copy everything from b_init to this new bunch
    # --------------------------------------------------------------------
    b = Bunch()
    b_init.copyBunchTo(b)

    # ---------------------------------------
    # track the bunch through the lattice
    # ---------------------------------------
    lattice.trackBunch(b)

    count += 1
    if count % 1000 and rank == 0:
        print("iter=", count)

print("Finish.")
