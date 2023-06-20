# -----------------------------------------------------
# Creates Grid3D for charge density and Poisson Solver
# for a pointlike charge (3D case). The field also includes
# the effect of neighboring bunches.
# Copmares with the exact result.
# -----------------------------------------------------
#

import sys
import math
import orbit.core

import orbit_mpi

from spacecharge import Grid3D
from spacecharge import PoissonSolverFFT3D

print("========= Start ===========")

sizeX = 64
sizeY = 64
sizeZ = 64
xMin = -5.0
xMax = +5.0
yMin = -5.5
yMax = +5.5
zMin = -6.0
zMax = +6.0
print("Solver 3D grid (x,y,z) : ", (sizeX, sizeY, sizeZ))
solver = PoissonSolverFFT3D(sizeX, sizeY, sizeZ, xMin, xMax, yMin, yMax, zMin, zMax)

# ---- this should be even number
nBunches = 2 * 1
lambda_dist = 15.0

solver.numExtBunches(nBunches)
solver.distBetweenBunches(lambda_dist)

# ---- we have to update Green function FFT after we changed nBunches or lambda
solver.updateGeenFunction()

print("debug nBunches=", solver.numExtBunches())
print("debug distance=", solver.distBetweenBunches())

gridRho = Grid3D(sizeX, sizeY, sizeZ)
gridRho.setGridX(xMin, xMax)
gridRho.setGridY(yMin, yMax)
gridRho.setGridZ(zMin, zMax)

gridPhi = Grid3D(sizeX, sizeY, sizeZ)
gridPhi.setGridX(xMin, xMax)
gridPhi.setGridY(yMin, yMax)
gridPhi.setGridZ(zMin, zMax)

chrage_pos_x = 1.0
chrage_pos_y = 1.0
chrage_pos_z = 1.0
charge = 1.0
gridRho.binValue(charge, chrage_pos_x, chrage_pos_y, chrage_pos_z)

print("charge position (x,y,z) =", (chrage_pos_x, chrage_pos_y, chrage_pos_z))


class Potential3D:
    """
    It calculates a potential from the central and several neighboring point-like charges
    if the charge distribution is periodic in z-direction with period lambda_dist.
    """

    def __init__(self, charge, nBunches, lambda_dist, chrage_pos_x, chrage_pos_y, chrage_pos_z):
        self.charge = charge
        self.nBunches = 2 * int(nBunches / 2)
        self.lambda_dist = lambda_dist
        self.r0_v = [chrage_pos_x, chrage_pos_y, chrage_pos_z]

    def getPotential(self, r_v):
        # print "debug =============== point=",r_v
        phi = 0.0
        for ind in range(-nBunches / 2, (nBunches / 2) + 1):
            x = r_v[0] - self.r0_v[0]
            y = r_v[1] - self.r0_v[1]
            z = r_v[2] - (self.r0_v[2] + ind * self.lambda_dist)
            phi += 1.0 / math.sqrt(x**2 + y**2 + z**2)
            # print "debug ind=",ind," (x,y,z)=",(x,y,z)," dist =",math.sqrt(x**2 + y**2 + z**2)
        phi *= self.charge
        return phi


potential3D = Potential3D(charge, nBunches, lambda_dist, chrage_pos_x, chrage_pos_y, chrage_pos_z)

print("Start solver.")

solver.findPotential(gridRho, gridPhi)

r_test = 3.0
n_angle_steps = 20
angle_step = 360.0 / (n_angle_steps - 1)
print("  i  theta      x       y        z         r          phi         phi_theory    ratio phi/theory  ")
for i in range(n_angle_steps):
    angle = math.pi * i * angle_step / 180.0
    x = r_test * math.cos(angle)
    y = r_test * math.sin(angle)
    z = -1.0
    phi = gridPhi.getValue(x, y, z)
    dist = math.sqrt(
        (chrage_pos_x - x) * (chrage_pos_x - x) + (chrage_pos_y - y) * (chrage_pos_y - y) + (chrage_pos_z - z) * (chrage_pos_z - z)
    )
    phi_th = potential3D.getPotential([x, y, z])
    theta = angle * 180.0 / math.pi
    print("", i, " %8.4f  %7.4f  %7.4f  %7.4f  %7.4f  %12.5g  %12.5g  %12.7g  " % (theta, x, y, z, dist, phi, phi_th, (phi / phi_th)))

print("Stop.")

sys.exit(0)

count = 0
while 1 < 2:
    solver.findPotential(gridRho, gridPhi)
    count += 1
    if count % 10 == 0:
        print("solved n=", count)
