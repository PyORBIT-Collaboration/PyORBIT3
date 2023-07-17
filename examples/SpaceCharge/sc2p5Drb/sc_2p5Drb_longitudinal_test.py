# -----------------------------------------------------
# Creates Space Charge Calculator based on Rick Baartman
# suggestion and tracks the test bunch through the calculator.
# This script will check the longitudinal kick.
# -----------------------------------------------------

import sys
import math
import random


from orbit.core import orbit_mpi

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalc2p5Drb

print("Start.")

# -----------------------------------------------------
# Make the Space Charge solver based on Rick Baartman
# suggestion.
# -----------------------------------------------------

sizeX = 64
sizeY = 64
sizeZ = 20
long_avg_n = 3
calc2p5d = SpaceChargeCalc2p5Drb(sizeX, sizeY, sizeZ)
calc2p5d.setLongAveragingPointsN(long_avg_n)
# ----------------------------------------------------
# make the bunch with a triangle shaped longitudinal distribution
# ----------------------------------------------------
b = Bunch()
bunch_radius = 0.005
bunch_length = 200.0

nParts = 100000

for ip in range(nParts):
    r = bunch_radius * math.sqrt(random.random())
    phi = 2 * math.pi * random.random()
    x = r * math.sin(phi)
    y = r * math.cos(phi)
    z = 0.5 * bunch_length * math.sqrt(random.random())
    if random.random() > 0.5:
        z = z - bunch_length / 2
    else:
        z = bunch_length / 2 - z
    b.addParticle(x, 0.0, y, 0.0, z, 0.0)

macroSize = 1.0e13
energy = 1.0
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(macroSize / nParticlesGlobal)
b.getSyncParticle().kinEnergy(energy)
gamma = b.getSyncParticle().gamma()
beta = b.getSyncParticle().beta()

print("gamma=", b.getSyncParticle().gamma())
print("beta=", b.getSyncParticle().beta())
print("bunchSize = ", b.getSize(), " global size=", nParticlesGlobal)
print("macroSize=", b.macroSize())
print("globalMacrosize=", macroSize)
print("mass=", b.mass())

# ---------------------------------------------------
# parameters of the Space Charge node:
# pipe radius and the length of space charge path
# ---------------------------------------------------

pipe_radius = 0.010
slice_length = 0.1

# ----------------------------------------------------
# Track the bunch through the SC calculator
# ----------------------------------------------------

# b.dumpBunch("pyorbit_bunch_test_in.dat")
print("Start Poisson Solver.")
calc2p5d.trackBunch(b, slice_length, pipe_radius)
print("Stop Poisson Solver.")
# b.dumpBunch("pyorbit_bunch_test_out.dat")

# -------------------------------------------------------
# Start of the analysis
# -------------------------------------------------------
longGrid = calc2p5d.getLongGrid()
longDerivGrid = calc2p5d.getLongDerivativeGrid()
# ----------------------------------------------------------------
# theoretical  and simulated dPz should be constant with z
# ----------------------------------------------------------------
slope_z_theory = (1.534698e-18 * slice_length / (gamma**2 * beta)) * b.mass() * (4 * macroSize / (bunch_length**2))

slope_z_avg = 0.0
slope2_z_avg = 0.0
count = 0
for ip in range(b.getSize()):
    z = b.z(ip)
    dP = b.dE(ip) / beta
    x = b.x(ip)
    y = b.y(ip)
    r = math.sqrt(x * x + y * y)
    coef_r = 1.0
    if r > bunch_radius:
        coef_r = 2 * math.log(pipe_radius / r)
    else:
        coef_r = 1.0 + 2 * math.log(pipe_radius / bunch_radius) - (r / bunch_radius) ** 2
    if z < -0.25 * bunch_length or z > 0.25 * bunch_length:
        slope_z_avg += dP * (math.fabs(z) / z) / coef_r
        slope2_z_avg += (dP / coef_r) ** 2
        count += 1
slope_z_avg /= count
slope2_z_avg /= count
slope_z_err = math.sqrt((slope2_z_avg - slope_z_avg * slope_z_avg))
print(" ========= theory and simulation should agree! ============== ")
print("particles slope delta(pz)            = %12.5g +- %12.5g " % (slope_z_avg, slope_z_err))
print("particles slope delta(pz) from theory= %12.5g" % slope_z_theory)

# ----------------------------------------------------------
# Check the long_avg_n -points quadratic interpolations for
# longitudinal derivative
# longGrid - Grid1D density distribution
# longDerivGrid - Grid1D density distribution derivative (normalized with 1/z_step)
# --- Grid1D has its own derivative method
# -----------------------------------------------------------
print("=========================================================")
z_step = longGrid.getGridZ(1) - longGrid.getGridZ(0)
for iz in range(sizeZ):
    z_1 = longGrid.getGridZ(iz)
    z_2 = longDerivGrid.getGridZ(iz)
    dz = z_2 - z_1
    val = longGrid.getValue(z_1)
    grad_1 = longGrid.calcGradient(z_1)
    grad_2 = longDerivGrid.getValue(z_1) * z_step
    st = "ind =", iz, " z[m]= %12.6g   dz= %12.6g    val= %12.6g    grad1,2 = ( %12.6g , %12.6g) " % (z_1, dz, val, grad_1, grad_2)
    print(st)
print("=========================================================")

# --------------------------------------------
# plot the graphs for a longitudinal density and derivative
# --------------------------------------------
nStep = 300

long_arr = []
for ix in range(nStep + 1):
    z = ix * bunch_length / nStep - bunch_length / 2.0
    long_rho = longGrid.getValue(z) / longGrid.getStepZ()
    long_arr.append((z, long_rho))

long_grad_arr = []
for ix in range(nStep + 1):
    z = ix * bunch_length / nStep - bunch_length / 2.0
    long_grad_rho = longDerivGrid.getValue(z)
    long_grad_arr.append((z, long_grad_rho))

# -------------------------------------------------
# this is the example of using the Gnuplot package
# Comment this part if you do not have this package installed.
# -------------------------------------------------
import Gnuplot

gLong = Gnuplot.Gnuplot()
gLong.title("Long Density as a function of z")
gLong("set style data lines")
gLong.plot(long_arr)

gLongGrad = Gnuplot.Gnuplot()
gLongGrad.title("Derivative of Long. Density ")
gLongGrad("set style data lines")
gLongGrad.plot(long_grad_arr)

input("Please press return to stop:\n")
# -------------------------------------------------

print("Stop.")
# orbit_mpi.finalize("Test is done!")
