# -----------------------------------------------------
# Creates Grid3D for charge density and Poisson Solver
# for a charged point (3D case)
# Compares with the exact result for the pint like charge
# at the one grid point. The theory and calculation should be
# the exactly same, because there is no smoothing of the grid or the
# gradient calculations. The maximal deviation from the theory should be
# on the level of the significant digits - around 10^-15
# -----------------------------------------------------
# We can scale the grids with respect to the solver by
# a constant coefficient. It is result of scaling of the Green
# function which is 1/r

import sys
import math

from orbit.core import orbit_mpi

from orbit.core.spacecharge import Grid3D, PoissonSolverFFT3D

print("Start.")

sizeX = 64
sizeY = 64
sizeZ = 64
xMin = -5.0
xMax = +5.0
yMin = -5.5
yMax = +5.5
zMin = -6.0
zMax = +6.0
print(" x,y,z sizes: ", sizeX, " ", sizeY, " ", sizeZ)
solver = PoissonSolverFFT3D(sizeX, sizeY, sizeZ, xMin, xMax, yMin, yMax, zMin, zMax)

scale_coeff = 1.3

gridRho = Grid3D(sizeX, sizeY, sizeZ)
gridRho.setGridX(scale_coeff * xMin, scale_coeff * xMax)
gridRho.setGridY(scale_coeff * yMin, scale_coeff * yMax)
gridRho.setGridZ(scale_coeff * zMin, scale_coeff * zMax)

gridPhi = Grid3D(sizeX, sizeY, sizeZ)
gridPhi.setGridX(scale_coeff * xMin, scale_coeff * xMax)
gridPhi.setGridY(scale_coeff * yMin, scale_coeff * yMax)
gridPhi.setGridZ(scale_coeff * zMin, scale_coeff * zMax)

chrage_ind_x = int(sizeX / 2)
chrage_ind_y = int(sizeY / 2)
chrage_ind_z = int(sizeZ / 2)

charge = 1.0
gridRho.setValue(charge, chrage_ind_x, chrage_ind_y, chrage_ind_z)

print("Start solver.")

solver.findPotential(gridRho, gridPhi)

# -------------------------------------------------------
# Now, let's test the results
# -------------------------------------------------------

x0 = gridRho.getGridX(chrage_ind_x)
y0 = gridRho.getGridY(chrage_ind_y)
z0 = gridRho.getGridZ(chrage_ind_z)

x_step = (gridRho.getMaxX() - gridRho.getMinX()) / (sizeX - 1)
y_step = (gridRho.getMaxY() - gridRho.getMinY()) / (sizeY - 1)
z_step = (gridRho.getMaxZ() - gridRho.getMinZ()) / (sizeZ - 1)
step_total = math.sqrt(x_step**2 + y_step**2 + z_step**2)


max_dev = 0.0
max_diff_ix = -1
max_diff_iy = -1
max_diff_iz = -1
max_grad_diff = 0.0
max_grad_diff_ix = -1
max_grad_diff_iy = -1
max_grad_diff_iz = -1
max_diff_phi = 0.0
max_diff_solv_phi = 0.0
for ix in range(sizeX):
    for iy in range(sizeY):
        for iz in range(sizeZ):
            if ix != chrage_ind_x or chrage_ind_y != iy or chrage_ind_z != iz:
                x = gridRho.getGridX(ix)
                y = gridRho.getGridY(iy)
                z = gridRho.getGridZ(iz)
                dist = math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)
                phi_th = charge / dist
                phi = gridPhi.getValueOnGrid(ix, iy, iz)
                diff = math.fabs(phi_th - phi)
                if diff > max_dev:
                    max_dev = diff
                    max_diff_ix = ix
                    max_diff_iy = iy
                    max_diff_iz = iz
                    max_diff_phi = phi_th
                    max_diff_solv_phi = phi
                # -------gradient---------------------------
                if dist > 3 * step_total:
                    (gradX, gradY, gradZ) = gridPhi.calcGradient(x, y, z)
                    gradX_Th = -charge * (x - x0) / (dist * dist * dist)
                    gradY_Th = -charge * (y - y0) / (dist * dist * dist)
                    gradZ_Th = -charge * (z - z0) / (dist * dist * dist)
                    grad_diff = 100.0 * math.sqrt((gradX_Th - gradX) ** 2 + (gradY_Th - gradY) ** 2 + (gradZ_Th - gradZ) ** 2)
                    grad_diff /= math.sqrt(gradX_Th**2 + gradY_Th**2 + gradZ_Th**2)
                    if grad_diff > max_grad_diff:
                        max_grad_diff = grad_diff
                        max_grad_diff_ix = ix
                        max_grad_diff_iy = iy
                        max_grad_diff_iz = iz

print("charge position ix =", chrage_ind_x)
print("charge position iy =", chrage_ind_y)
print("charge position iz =", chrage_ind_z)
print("potential at the cahrge position phi =", gridPhi.getValueOnGrid(chrage_ind_x, chrage_ind_y, chrage_ind_z))
print("================================")
print("phi theory at max deviation =", max_diff_phi)
print("phi Solver at max deviation =", max_diff_solv_phi)
print("max Solver-Theory deviation =", max_dev)
print("max diff position ix = ", max_diff_ix)
print("max diff position iy = ", max_diff_iy)
print("max diff position iz = ", max_diff_iz)
print("================================")
print("max gradient deviation [%] =", "%5.3f" % max_grad_diff)
print("max grad diff position ix = ", max_grad_diff_ix)
print("max grad diff position iy = ", max_grad_diff_iy)
print("max grad diff position iz = ", max_grad_diff_iz)
print("Stop.")
sys.exit()

count = 0
while 1 < 2:
    solver.findPotential(gridRho, gridPhi)
    count += 1
    if count % 10 == 0:
        print("solved n=", count)
