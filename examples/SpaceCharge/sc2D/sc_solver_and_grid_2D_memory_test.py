# -----------------------------------------------------
# Creates Grid2D and Boundary2D to check the memory leak
# -----------------------------------------------------
import sys
import orbit.core

from spacecharge import Grid2D
from spacecharge import PoissonSolverFFT2D

print("Start.")

count = 0
while 1 < 2:
    grid0 = Grid2D(100, 200)
    grid1 = Grid2D(100, 200)
    solver = PoissonSolverFFT2D(100, 200)
    solver.findPotential(grid0, grid1)
    count = count + 1
    if count % 100 == 0:
        print("count=", count)

print("Stop.")

sys.exit(1)
