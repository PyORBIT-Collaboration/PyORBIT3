# -----------------------------------------------------
# Creates Grid2D for charge density and Poisson Solver
# for a charged string (2D case)
# Copmares with the exact result.
# -----------------------------------------------------
# The one very useful property of the FFT Poisson solver - scalability
# you can change the absolute value of steps X and Y and you will get
# the right results as soon you keep stepX/stepY constant.
#

import sys
import math
import pytest

from orbit.core import orbit_mpi

from orbit.core.spacecharge import Grid2D, PoissonSolverFFT2D

sizeX = 200
sizeY = 200
xMin = -5.0
xMax = +5.0
yMin = -5.0
yMax = +5.0

scale_coeff = 3.0
solver = PoissonSolverFFT2D(sizeX, sizeY, xMin / scale_coeff, xMax / scale_coeff, yMin / scale_coeff, yMax / scale_coeff)

gridRho = Grid2D(sizeX, sizeY, xMin, xMax, yMin, yMax)
gridPhi = Grid2D(sizeX, sizeY, xMin, xMax, yMin, yMax)

chrage_pos_x = 2.5
chrage_pos_y = 1.0
charge = 1.0
gridRho.binValue(charge, chrage_pos_x, chrage_pos_y)

solver.findPotential(gridRho, gridPhi)

r_test = 4.0
n_angle_steps = 10
angle_step = 360.0 / (n_angle_steps - 1)
comparison_string = "  i    x       y         phi         phi_theory    ratio phi/theory\n"
for i in range(n_angle_steps):
    angle = math.pi * i * angle_step / 180.0
    x = r_test * math.cos(angle)
    y = r_test * math.sin(angle)
    phi = gridPhi.getValue(x, y)
    dist = (chrage_pos_x - x) * (chrage_pos_x - x) + (chrage_pos_y - y) * (chrage_pos_y - y)
    dist = math.sqrt(dist)
    phi_th = -math.log(dist)
    comparison_string += "{} {:7.4f} {:7.4f} {:12.5g} {:12.5g} {:12.7g}".format(i, x, y, phi, phi_th, (phi / phi_th)) + "\n"

print(comparison_string)


def test_Poisson():
    expected_out = """  i    x       y         phi         phi_theory    ratio phi/theory
0  4.0000  0.0000     -0.58933     -0.58933    0.9999996
1  3.0642  2.5712     -0.51245     -0.51245    0.9999998
2  0.6946  3.9392      -1.2382      -1.2382            1
3 -2.0000  3.4641      -1.6352      -1.6352            1
4 -3.7588  1.3681      -1.8357      -1.8357            1
5 -3.7588 -1.3681      -1.9009      -1.9009            1
6 -2.0000 -3.4641      -1.8467      -1.8467            1
7  0.6946 -3.9392      -1.6599      -1.6599            1
8  3.0642 -2.5712      -1.2852      -1.2852            1
9  4.0000 -0.0000     -0.58933     -0.58933    0.9999996\n"""
    assert comparison_string == expected_out
