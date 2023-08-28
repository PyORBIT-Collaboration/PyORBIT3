# -----------------------------------------------------
# Creates Space Charge Calculator based on Rick Baartman
# suggestion and tests the getRho, getPhi, getLong Grids
# methods to avoid memory leaks
# -----------------------------------------------------

import sys
import math
import random

from orbit.core import orbit_mpi

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalc2p5Drb

print("Start.")

sizeX = 128
sizeY = 128
sizeZ = 2
calc2p5d = SpaceChargeCalc2p5Drb(sizeX, sizeY, sizeZ)
rho = calc2p5d.getRhoGrid()
phi = calc2p5d.getPhiGrid()
longGrid = calc2p5d.getLongGrid()
print("sX=", rho.getSizeX(), " sY=", phi.getSizeY(), " sZ=", longGrid.getSizeZ())

count = 0
while 1 < 2:
    rho = calc2p5d.getRhoGrid()
    phi = calc2p5d.getPhiGrid()
    longGrid = calc2p5d.getLongGrid()
    count += 1
    if count % 100000 == 0:
        print("i=", count)

del calc2p5d

print("Stop.")
