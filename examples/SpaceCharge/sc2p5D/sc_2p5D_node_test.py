# -----------------------------------------------------
# Creates Space Charge Calculator Node
# -----------------------------------------------------
import sys
import math
import random
import orbit.core

random.seed(10)

from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit.teapot import teapot
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5D, Boundary2D

print("Start.")

# make a bunch
b = Bunch()
bunch_radius = 0.005
bunch_length = 500e-9 * 3e8 * 0.87502565
nParts = 100000

# Use random radius
for ip in range(nParts):
    r = bunch_radius * math.sqrt(random.random())
    phi = 2 * math.pi * random.random()
    x = r * math.sin(phi)
    y = r * math.cos(phi)
    z = bunch_length * 0.5 * (1.0 - 2 * random.random())
    b.addParticle(x, 0.0, y, 0.0, z, 0.0)

# set bunch parameters
macroSize = 1.56e13
energy = 0.08
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(macroSize / nParticlesGlobal)
syncPart = b.getSyncParticle()
syncPart.kinEnergy(energy)

# make a Teapot lattice
elem1 = teapot.DriftTEAPOT("drift1")
elem2 = teapot.QuadTEAPOT("quad1")
elem3 = teapot.QuadTEAPOT("quad2")

teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
teapot_lattice.addNode(elem1)
# teapot_lattice.addNode(elem2)
# teapot_lattice.addNode(elem3)

# set node prameters
elem1.setLength(0.2)
elem1.setnParts(2)
elem2.setLength(1.0)
elem2.setnParts(5)
elem2.addParam("kq", -0.7)
elem3.setLength(2.0)
elem3.setnParts(5)
elem3.addParam("kq", +0.7)

teapot_lattice.initialize()

# make 2.5D space charge calculator
sizeX = 64
sizeY = 64
sizeZ = 20
calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)

# set boundary
nBoundaryPoints = 100
N_FreeSpaceModes = 20
R_Boundary = 0.008
boundary = Boundary2D(nBoundaryPoints, N_FreeSpaceModes)

for i in range(nBoundaryPoints):
    x = R_Boundary * math.cos((2.0 * math.pi / (nBoundaryPoints - 1)) * i)
    y = R_Boundary * math.sin((2.0 * math.pi / (nBoundaryPoints - 1)) * i)
    boundary.setBoundaryPoint(i, x, y)
boundary.initialize()

# -------------------------------------------------------------------------
"""
For the Boundary2D class it is possible to have the predefined shapes
Circle - the last parameter will be the diameter of the circle
Ellipse - there will be two parameters at the end - 2*a and 2*b where a,b are semi-axises
Rectangle - there will be two parameters at the end - horizontal and vertical sizes
"""
# -------------------------------------------------------------------------
# boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2*R_Boundary)

sc_path_length_min = 0.05

scNode_arr = scLatticeModifications.setSC2p5DAccNodes(teapot_lattice, sc_path_length_min, calc2p5d)
# scNode_arr = scLatticeModifications.setSC2p5DAccNodes(teapot_lattice, sc_path_length_min, calc2p5d, boundary)

# track lattice
# teapot_lattice.trackBunch(b)
slice_length = 0.1
calc2p5d.trackBunch(b, slice_length)

# -------------------------------------
# momentum change
# coeff is: 2*r0*L*lambda/e*(gamma^3*beta^2)
# r/a^2 and grad/macrosize
# -------------------------------------
xyp_coeff = 2 * b.classicalRadius() * b.charge() ** 2 * 0.1 / (b.getSyncParticle().gamma() ** 3 * b.getSyncParticle().beta() ** 2)
for ip in range(10):
    x = b.x(ip)
    y = b.y(ip)
    r = math.sqrt(x * x + y * y)
    theta = math.atan2(y, x)
    xp = b.xp(ip)
    yp = b.yp(ip)
    xp_calc = math.cos(theta) * r * xyp_coeff * macroSize / (bunch_length * bunch_radius**2)
    yp_calc = math.sin(theta) * r * xyp_coeff * macroSize / (bunch_length * bunch_radius**2)
    xp_error = (xp - xp_calc) * 100 / xp_calc
    yp_error = (yp - yp_calc) * 100 / yp_calc
    print("debug xyp_coeff=", xyp_coeff)
    print(
        "r=%10.5g" % r,
        "x=%10.5g" % x,
        "y=%10.5g" % y,
        "xp = %10.5g " % xp,
        " xp_theory = %10.5g " % xp_calc,
        " xp_error = %10.5g " % xp_error,
        "%",
    )
    print(
        "r=%10.5g" % r,
        "x=%10.5g" % x,
        "y=%10.5g" % y,
        "yp = %10.5g " % yp,
        " yp_theory = %10.5g " % yp_calc,
        " xp_error = %10.5g " % yp_error,
        "%",
    )

# printout lattice
nodes = teapot_lattice.getNodes()
for node in nodes:
    print("node=", node.getName(), " s start,stop = %4.3f %4.3f " % teapot_lattice.getNodePositionsDict()[node])
    print("There are ", node.getNumberOfBodyChildren(), " children nodes.")
    childnodes = node.getChildNodes(AccNode.BODY)
    for childnode in childnodes:
        print("node child=", childnode.getName(), " length=", childnode.getLengthOfSC())

print("Finish.")
