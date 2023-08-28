# -----------------------------------------------------
# Grid1D gradient test
# User can use grid1D.calcGradient(z) or grid1D.calcGradientSmoothed(z)
# the smoothed version should be more precise.
# -----------------------------------------------------

import sys
import math

from orbit.core.spacecharge import Grid1D


print("Start.")

sizeZ = 10
zMin = 1.0
zMax = +4.0

grid1D = Grid1D(sizeZ, zMin, zMax)


def Func(z):
    return 1.0 / (z * z)


def FuncGrad(z):
    return -2.0 / (z * z * z)


for iz in range(sizeZ):
    z = grid1D.getGridZ(iz)
    val = Func(z)
    grid1D.setValue(val, iz)

diff_max = 0.0
max_point = 0.0
# --- we exclude the last point because Grid1D wraps the start and end of the array
# --- which is used for ring linear space charge calculations
for iz in range(sizeZ - 1):
    z = grid1D.getGridZ(iz)
    grad0 = FuncGrad(z)
    # grad = grid1D.calcGradient(z)
    grad = grid1D.calcGradientSmoothed(z)
    diff = math.sqrt((grad - grad0) ** 2)
    print("debug z= %10.4f " % z, " exactGrad= %12.5g  grad = %12.5g" % (grad0, grad))
    if diff > diff_max:
        diff_max = diff
        max_point = z


print("max diff =", diff_max)
print("max point =", max_point)
z = max_point
print("Func at max point =", Func(z))
print("Grad at max point =", FuncGrad(z))

print("Stop.")
