# -----------------------------------------------------
# Creates Grid2D for charge density and potential
# Boundary2D instance for solving Poisson problem
# for a charged string (2D case) with a boundary induced potential
# Compares with the exact result.
# phi(r) = ln(abs(r-a)/abs(r-a')) where r,a,a' are vectors
#           abs(a') = R^2/abs(a)  a' is parallel to a
#           a - specifies the position of the charge
# -----------------------------------------------------
# The one very useful property of the FFT Poisson solver - scalability
# you can change the absolute value of steps X and Y and you will get
# the right results as soon you keep stepX/stepY constant.

import sys
import math
import orbit.core

import orbit_mpi

from spacecharge import Grid2D
from spacecharge import Boundary2D
from spacecharge import PoissonSolverFFT2D

print("Start.")

sizeX = 200
sizeY = 200

xMin = -5.0
xMax = +5.0
yMin = -5.0
yMax = +5.0

nBoundaryPoints = 100
N_FreeSpaceModes = 20
R_Boundary = 5.0

gridRho = Grid2D(sizeX, sizeY, xMin, xMax, yMin, yMax)
gridPhi = Grid2D(sizeX, sizeY, xMin, xMax, yMin, yMax)

# The one very useful property of the FFT Poisson solver - scalability
# you can change the absolute value of steps X and Y and you will get
# the right results as soon you keep stepX/stepY constant.
solver = PoissonSolverFFT2D(sizeX, sizeY, xMin / 2.0, xMax / 2.0, yMin / 2.0, yMax / 2.0)

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
# boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",xMax-xMin)

print("shape name=", boundary.getShapeName())

chrage_pos_x = 2.5
chrage_pos_y = 0.0
chrage_pos_a = math.sqrt(chrage_pos_x * chrage_pos_x + chrage_pos_y * chrage_pos_y)
charge = 1.0
gridRho.binValue(charge, chrage_pos_x, chrage_pos_y)

R = R_Boundary
a_prime = R * R / chrage_pos_a
a_prime_x = a_prime * chrage_pos_x / chrage_pos_a
a_prime_y = a_prime * chrage_pos_y / chrage_pos_a

solver.findPotential(gridRho, gridPhi)

boundary.addBoundaryPotential(gridRho, gridPhi)

# -----potential delta-------------------------------
# our potential on the wall is zero, but the exact
# solution is not zero. There is a constant difference
# between these potentials
# ---------------------------------------------------
dist = (chrage_pos_x - R) * (chrage_pos_x - R) + chrage_pos_y * chrage_pos_y
dist = math.sqrt(dist)
dist_prime = (a_prime_x - R) * (a_prime_x - R) + a_prime_y * a_prime_y
dist_prime = math.sqrt(dist_prime)
phi_delta = -math.log(dist / dist_prime)


r_test = 4.0
n_angle_steps = 10
angle_step = 360.0 / (n_angle_steps - 1)
print("  i    x       y         phi         phi_theory    ratio phi/theory  ")
for i in range(n_angle_steps):
    angle = math.pi * i * angle_step / 180.0
    x = r_test * math.cos(angle)
    y = r_test * math.sin(angle)
    phi = gridPhi.getValue(x, y)
    dist = (chrage_pos_x - x) * (chrage_pos_x - x) + (chrage_pos_y - y) * (chrage_pos_y - y)
    dist = math.sqrt(dist)
    dist_prime = (a_prime_x - x) * (a_prime_x - x) + (a_prime_y - y) * (a_prime_y - y)
    dist_prime = math.sqrt(dist_prime)
    phi_th = -math.log(dist / dist_prime) - phi_delta
    ratio = 0.0
    if phi_th != 0.0:
        ratio = phi / phi_th
    print("", i, " %7.4f  %7.4f  %12.5g  %12.5g  %12.7g  " % (x, y, phi, phi_th, ratio))

print("Stop.")
