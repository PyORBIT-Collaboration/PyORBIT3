import sys
import math
import orbit.core

import orbit_mpi

from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5D
from spacecharge import PoissonSolverFFT2D
from spacecharge import Boundary2D

print("Start.")

sizeX = 128
sizeY = 128
sizeZ = 2
calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ, 1)

charge = 1.0

b = Bunch()
b.addPartAttr("macrosize")

macroSize = 1.0e-13

energy = 1.0
b.getSyncParticle().kinEnergy(energy)

b.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
b.addParticle(0.001, 0.0, 0.0, 0.0, 0.0, 0.0)
b.addParticle(-0.001, 0.0, 0.0, 0.0, 0.0, 0.0)
b.addParticle(0.0, 0.0, 0.001, 0.0, 0.0, 0.0)
b.addParticle(0.0, 0.0, -0.001, 0.0, 0.0, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, 0.001, 0.0)
b.addParticle(0.0, 0.0, 0.0, 0.0, -0.001, 0.0)

b.partAttrValue("macrosize", 0, 0, macroSize)
b.partAttrValue("macrosize", 1, 0, 0.0)
b.partAttrValue("macrosize", 2, 0, 0.0)
b.partAttrValue("macrosize", 3, 0, 0.0)
b.partAttrValue("macrosize", 4, 0, 0.0)
b.partAttrValue("macrosize", 5, 0, 0.0)
b.partAttrValue("macrosize", 6, 0, 0.0)

nParts = b.getSize()
macroSize = macroSize / nParts
b.macroSize(macroSize)

print("bunchSize = ", b.getSize())
print("macroSize=", b.macroSize())
print("mass=", b.mass())

slice_length = 0.1

print("without boundary.")
b.dumpBunch("pyorbit_bunch_test.in")
calc2p5d.trackBunch(b, slice_length)
b.dumpBunch("pyorbit_bunch_test.out")

print("with boundary.")
nBoundaryPoints = 128
N_FreeSpaceModes = 32
R_Boundary = 0.001
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

calc2p5d.trackBunch(b, slice_length, boundary)
b.dumpBunch("pyorbit_bunch_boundary_test.out")

print("Stop.")
