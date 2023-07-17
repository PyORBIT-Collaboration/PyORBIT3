# -----------------------------------------------------
# Grid2D values binning test
# User can change number of grid points (sizeZ) and the type of binning
# (grid1D.binValue or grid1D.binValueSmoothed)
# The test should be performed for sizeZ = 1,2,3,  20
# You expect to see uniform distribution
# -----------------------------------------------------

import sys
import math
import random

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import Grid1D
from orbit.diagnostics import profiles

random.seed(100)

print("Start.")

sizeZ = 5
zMin = 1.0
zMax = +4.0

extent = 0.2

grid1D = Grid1D(sizeZ, zMin, zMax)

# -------------------------------------------
# Bin the uniformly distributed values along z
# -------------------------------------------
nPlotPoints = 100000
x_polt_min = zMin - extent
x_polt_max = zMax + extent

value = 1.0

bunch_test = Bunch()
macrosize = 1.0
bunch_test.macroSize(macrosize)

for ind_p in range(nPlotPoints):
    x = random.uniform(x_polt_min, x_polt_max)
    bunch_test.addParticle(x, 2 * x, 3 * x, 4 * x, 5 * x, 6 * x)
    grid1D.binValue(value, x)
    # grid1D.binValueSmoothed(value,x)

print("nParts in bunch=", bunch_test.getSizeGlobal())

# ---------------------------------------------------
# Plot the density - it should be almost uniform
# --------------------------------------------------

x_arr = []
y_arr = []
for ind in range(sizeZ):
    x = grid1D.getGridZ(ind)
    y = grid1D.getValueOnGrid(ind)
    y_value = grid1D.getValueOnGrid(ind)
    print("debug z=", x, " rho=", y, " rho_value=", y_value, " err=", math.sqrt(y_value))
    x_arr.append(x)
    y_arr.append(y)

# ---- dump histogram into the file. profiles uses grid1D.binBunch(bunch)
# ---- so, the histogram.dat should have the same information as previous
# ---- debug printing
coord = "x"
histogram = "histogram.dat"
profiles(bunch_test, coord, histogram, sizeZ, zMin, zMax)

# ---- check value interpolation
print("debug ---check value interpolation ----")
print("debug range zMin=", zMin, " zMax=", zMax)
nPoints = 13
step = (zMax - zMin) / (nPoints - 1)
zMinNew = zMin - step
for ind in range(nPoints + 2):
    x = zMinNew + ind * step
    y = grid1D.getValue(x)
    # y = grid1D.getValueSmoothed(x)
    print("debug ind=", ind, " z=", x, " rho=", y)

# ---------------------------
# Plot part
# ---------------------------
import matplotlib.pyplot as plt

plt.plot(x_arr, y_arr)
plt.ylabel("rho(x)")
plt.xlabel("x")

plt.show()
