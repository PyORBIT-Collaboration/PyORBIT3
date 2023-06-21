# -----------------------------------------------------
# This example tracks the particles distributed
# uniformly in a cylinder through a drift.
# The space charge is calculated by using 2.5D SC node.
# The change in the particles transverse momentum compared
# with the analytical formula.
# -----------------------------------------------------
import sys
import math
import random
import orbit.core

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.teapot import teapot
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5D, Boundary2D

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
nParts = 500000
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
    scLatticeModifications.setSC2p5DAccNodes(teapot_lattice, sc_path_length_min, spaceChargeCalculator2p5D)
    return teapot_lattice


# ----------------------------------------------
# make 2.5D space charge calculator
# ----------------------------------------------
sizeX = 256  # number of grid points in horizontal direction
sizeY = 256  # number of grid points in vertical direction
sizeZ = 1  # number of longitudinal slices in the 2.5D space charge solver
calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)

lattice_length = 2.0  # the length of the drift
n_parts = 1  # number of parts on what the drift will be chopped, or the number of SC nodes
lattice = getLattice(lattice_length, n_parts, calc2p5d)

# -------------------------------
#  Lattice is ready
# -------------------------------

nodes = lattice.getNodes()
for node in nodes:
    print("node=", node.getName(), " s start,stop = %4.3f %4.3f " % lattice.getNodePositionsDict()[node])
    print("There are ", node.getNumberOfBodyChildren(), " child nodes.")


# --------------------------------------------------------------------
# Create a new bunch, copy everything from b_init to this new bunch
# Later we will compare momentum kicks between two bunches
# --------------------------------------------------------------------
b = Bunch()
b_init.copyBunchTo(b)

# ---------------------------------------
# track the bunch through the lattice
# ---------------------------------------
lattice.trackBunch(b)

# ---------------------------------------
#  Now compare bunches and the theory
#  The theory dp = r * 2 * r0 * total_macrosize * charge^2 * drift_length / (gamma^3*beta^2*bunch_radius^2*bunch_length)
#  r0 is a classical radius of the particle
# ---------------------------------------
n_graph_points = 1000
n_step = int(b.getSize() / n_graph_points) + 1

r_dp_arr = []
r_dp_theory_arr = []
for ip in range(0, b.getSize(), n_step):
    z = b.z(ip)
    x = b.x(ip)
    y = b.y(ip)
    xp = b.xp(ip)
    yp = b.yp(ip)
    x0 = b_init.x(ip)
    y0 = b_init.y(ip)
    xp0 = b_init.xp(ip)
    yp0 = b_init.yp(ip)
    r = math.sqrt(x0**2 + y0**2)
    dp = math.sqrt((xp - xp0) ** 2 + (yp - yp0) ** 2)
    dp_th = (
        2.0
        * r
        * b.classicalRadius()
        * total_macroSize
        * (b.charge()) ** 2
        * lattice_length
        / (b.getSyncParticle().gamma() ** 3 * b.getSyncParticle().beta() ** 2 * bunch_radius**2 * bunch_length)
    )
    # print " %12.5e %12.5e %12.5e  %12.5e "%(z,r,dp,dp_th)
    r_dp_arr.append([r, dp])
    r_dp_theory_arr.append([r, dp_th])

# this is the example of using the Gnuplot package
import Gnuplot

g = Gnuplot.Gnuplot(debug=1)
g.title("dP vs r ")
g("set data style points")
g.plot(r_dp_arr, r_dp_theory_arr)
input("Please press return to stop:\n")

print("Finish.")
