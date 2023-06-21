# -----------------------------------------------------
# Creates Space Charge Calculator
# test Boundary
# -----------------------------------------------------
import sys
import math
import random
import orbit.core

random.seed(10)
import orbit_mpi

from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5D
from spacecharge import Boundary2D

print("Start.")

# Make a SC solver
sizeX = 64
sizeY = 64
sizeZ = 20
calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)

# Make bunch
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

macroSize = 1.56e13
energy = 0.08
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(macroSize / nParticlesGlobal)
b.getSyncParticle().kinEnergy(energy)

# bunchExtremaCalc->getExtremaXYZ

# Make a boundary
nBoundaryPoints = 100
N_FreeSpaceModes = 20
R_Boundary = 0.005
boundary = Boundary2D(nBoundaryPoints, N_FreeSpaceModes)

for i in range(nBoundaryPoints):
    x = R_Boundary * math.cos((2.0 * math.pi / (nBoundaryPoints - 1)) * i)
    y = R_Boundary * math.sin((2.0 * math.pi / (nBoundaryPoints - 1)) * i)
    boundary.setBoundaryPoint(i, x, y)
boundary.initialize()
b_maxx = boundary.getMaxX()
b_minx = boundary.getMinX()
b_maxy = boundary.getMaxY()
b_miny = boundary.getMinY()
print("MaxX=", b_maxx, " MinX=", b_minx, " MaxY=", b_maxy, " MinY=", b_miny)


# -------------------------------------------------------------------------
"""
For the Boundary2D class it is possible to have the predefined shapes
Circle - the last parameter will be the diameter of the circle
Ellipse - there will be two parameters at the end - 2*a and 2*b where a,b are semi-axises
Rectangle - there will be two parameters at the end - horizontal and vertical sizes
"""
# -------------------------------------------------------------------------
# boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2*R_Boundary)


# Shape boundary test
# Choose from Circle/Ellipse/Rectangle
shapeboundary = Boundary2D(nBoundaryPoints, N_FreeSpaceModes, "Circle", 0.01, 0.01)
print("shape name=", boundary.getShapeName())

# Set SC node parameters
pipe_radius = 0.02
slice_length = 0.1

# Track and analysis
# analysis
rhoGrid = calc2p5d.getRhoGrid()
phiGrid = calc2p5d.getPhiGrid()
x = bunch_radius / 2.0
y = 0.0
r = math.sqrt(x * x + y * y)

rho_theory = macroSize * 4.0 / (math.pi * rhoGrid.getSizeX() * rhoGrid.getSizeY())
phi_theory = macroSize * r**2 / (bunch_radius**2)
grad_theory = macroSize * 2 * r / (bunch_radius**2)

# without boundary
# print "without boundary."
# calc2p5d.trackBunch(b,slice_length)

# with boundary
print("with boundary.")
# calc2p5d.trackBunch(b,slice_length,boundary)
calc2p5d.trackBunch(b, slice_length, shapeboundary)

rho = rhoGrid.getValue(x, y)
phi = 2 * (phiGrid.getValue(x, y) - phiGrid.getValue(0.0, 0.0))
(ex, ey) = phiGrid.calcGradient(x, y)
grad = 2 * math.sqrt(ex * ex + ey * ey)

print("r=", r, " rho  = %12.5g " % rho, "  rho_theory = %12.5g " % rho_theory)
print("r=", r, " phi  = %12.5g " % phi, "  phi_theory = %12.5g " % phi_theory)
print("r=", r, " grad = %12.5g " % grad, " grad_theory = %12.5g " % grad_theory)

# -------------------------------------
# momentum change
# coeff is: 2*r0*L*lambda/e*(gamma^3*beta^2)
# r/a^2 and grad/macrosize
# -------------------------------------
xyp_coeff = 2 * b.classicalRadius() * b.charge() ** 2 * slice_length / (b.getSyncParticle().gamma() ** 3 * b.getSyncParticle().beta() ** 2)
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
