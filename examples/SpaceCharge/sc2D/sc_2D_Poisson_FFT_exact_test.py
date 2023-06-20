# -----------------------------------------------------
# Creates Grid2D for charge density and Poisson Solver
# for a charged string (2D case)
# Compares with the exact result for the pint like (2D - string) charge
# at the one grid point. The theory and calculation should be
# the exactly same, because there is no smoothing of the grid or the
# gradient calculations. The maximal deviation from the theory should be
# on the level of the significant digits - around 10^-15
# -----------------------------------------------------
# The one very useful property of the FFT Poisson solver - scalability
# you can change the absolute value of steps X and Y and you will get
# the right results as soon you keep stepX/stepY constant.
#

import sys
import math
import orbit.core

import orbit_mpi

from spacecharge import Grid2D
from spacecharge import PoissonSolverFFT2D

print("Start.")

sizeX = 200
sizeY = 200
xMin = -5.0
xMax = +5.0
yMin = -5.0
yMax = +5.0

# let's make the solver
scale_coeff = 3.0
solver = PoissonSolverFFT2D(sizeX, sizeY, xMin / scale_coeff, xMax / scale_coeff, yMin / scale_coeff, yMax / scale_coeff)

# let's make the grids for the charge density and the potential
gridRho = Grid2D(sizeX, sizeY, xMin, xMax, yMin, yMax)
gridPhi = Grid2D(sizeX, sizeY, xMin, xMax, yMin, yMax)

# set the charge on the density grid point
ind_x = int(sizeX / 2)
ind_y = int(sizeY / 2)
charge = 1.0
gridRho.setValue(charge, ind_x, ind_y)

# solve the Poisson eqation. The results are in the potential Grid2D - gridPhi
solver.findPotential(gridRho, gridPhi)

# -------------------------------------------------------
# Now, let's test the results
# -------------------------------------------------------

x0 = gridRho.getGridX(ind_x)
y0 = gridRho.getGridY(ind_y)

x_step = (gridRho.getMaxX() - gridRho.getMinX()) / (sizeX - 1)
y_step = (gridRho.getMaxY() - gridRho.getMinY()) / (sizeY - 1)
step_total = math.sqrt(x_step**2 + y_step**2)

max_dev = 0.0
max_diff_ix = -1
max_diff_iy = -1
max_grad_diff = 0.0
max_grad_diff_ix = -1
max_grad_diff_iy = -1
max_diff_phy = 0
for ix in range(sizeX):
    for iy in range(sizeY):
        if ix != ind_x or ind_y != iy:
            x = gridRho.getGridX(ix)
            y = gridRho.getGridY(iy)
            dist = math.sqrt((x - x0) ** 2 + (y - y0) ** 2)
            # -------value---------------------------
            phi_th = -charge * math.log(dist)
            phi = gridPhi.getValueOnGrid(ix, iy)
            diff = math.fabs(phi_th - phi)
            if diff > max_dev:
                max_dev = diff
                max_diff_ix = ix
                max_diff_iy = iy
                max_diff_phy = phi_th
            # -------gradient---------------------------
            if dist > 2 * step_total:
                (gradX, gradY) = gridPhi.calcGradient(x, y)
                gradX_Th = -charge * (x - x0) / (dist * dist)
                gradY_Th = -charge * (y - y0) / (dist * dist)
                grad_diff = 100.0 * math.sqrt((gradX_Th - gradX) ** 2 + (gradY_Th - gradY) ** 2) / math.sqrt(gradX_Th**2 + gradY_Th**2)
                if grad_diff > max_grad_diff:
                    max_grad_diff = grad_diff
                    max_grad_diff_ix = ix
                    max_grad_diff_iy = iy

print("charge position ix =", ind_x)
print("charge position iy =", ind_y)
print("potential at the cahrge position phi =", gridPhi.getValueOnGrid(ind_x, ind_y))
print("================================")
print("max deviation          =", max_dev)
print("phi theory =", max_diff_phy)
print("max diff position ix = ", max_diff_ix)
print("max diff position iy = ", max_diff_iy)
print("================================")
print("max gradient deviation [%] =", "%5.3f" % max_grad_diff)
print("max grad diff position ix = ", max_grad_diff_ix)
print("max grad diff position iy = ", max_grad_diff_iy)
print("Stop.")
