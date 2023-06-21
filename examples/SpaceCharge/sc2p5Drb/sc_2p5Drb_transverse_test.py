# -----------------------------------------------------
# Creates Space Charge Calculator based on Rick Baartman
# suggestion and tracks the test bunch through the calculator.
# This is a test for the transverse space charge.
# The distribution is uniform in the longitudinal direction.
# -----------------------------------------------------

import sys
import math
import random
import orbit.core

import orbit_mpi

from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5Drb

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
# make the bunch
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
    z = bunch_length * 0.5 * (1.0 - 2 * random.random())
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
print("Stop Poisson Solver 0.")
# b.dumpBunch("pyorbit_bunch_test_out.dat")

# ----------------------------------------------
# check if Space Charge Calculator has a memory
# ----------------------------------------------
for ip in range(nParts):
    b.px(ip, 0.0)
    b.py(ip, 0.0)
    b.dE(ip, 0.0)

calc2p5d.trackBunch(b, slice_length, pipe_radius)
print("Stop Poisson Solver 1.")
for ip in range(nParts):
    b.px(ip, 0.0)
    b.py(ip, 0.0)
    b.dE(ip, 0.0)

calc2p5d.trackBunch(b, slice_length, pipe_radius)
print("Stop Poisson Solver 2.")

# -------------------------------------------------------
# Start of the analysis
# -------------------------------------------------------
rhoGrid = calc2p5d.getRhoGrid()
phiGrid = calc2p5d.getPhiGrid()

x = bunch_radius / 2.0
y = 0.0
r = math.sqrt(x * x + y * y)

# --------------------------------------------------------
# rhoGrid.getSizeX() and rhoGrid.getSizeY() - number of grid points in X and Y
# --------------------------------------------------------
rho = rhoGrid.getValue(x, y)
rho_theory = macroSize * 4.0 / (math.pi * rhoGrid.getSizeX() * rhoGrid.getSizeY())
print("r=", r, " rho  = %12.5g " % rho, "  rho_theory = %12.5g " % rho_theory)

# --------------------------------------------------------------------------
# remember that the potential is an abstract potential: lambda*ln(r)
# The potential for a charged string in CGS is 2*lambda*ln(r)
# --------------------------------------------------------------------------
phi = 2 * (phiGrid.getValue(x, y) - phiGrid.getValue(0.0, 0.0))
phi_theory = macroSize * r**2 / (bunch_radius**2)
print("r=", r, " phi  = %12.5g " % phi, "  phi_theory = %12.5g " % phi_theory)

(ex, ey) = phiGrid.calcGradient(x, y)
grad = 2 * math.sqrt(ex * ex + ey * ey)
grad_theory = macroSize * 2 * r / (bunch_radius**2)
print("r=", r, " grad = %12.5g " % grad, " grad_theory = %12.5g " % grad_theory)

# ----------------------------------------------------------------
# theoretical  and simulated coeff delta(r_prime/r)
# ----------------------------------------------------------------
slope_theory = (2.0 * 1.534698e-18 * slice_length / (gamma**3 * beta**2)) * (macroSize / bunch_length) / (bunch_radius**2)

slope_avg = 0.0
slope2_avg = 0.0
count = 0
for ip in range(b.getSize()):
    x = b.x(ip)
    y = b.y(ip)
    xp = b.xp(ip)
    yp = b.yp(ip)
    z = b.z(ip)
    dE = b.dE(ip)
    r = math.sqrt(x * x + y * y)
    p = math.sqrt(xp * xp + yp * yp)
    scalar_product = x * xp + y * yp
    if scalar_product != 0.0:
        scalar_product = math.fabs(scalar_product) / scalar_product
    p = p * scalar_product
    if r > 0.5 * bunch_radius and r < 0.9 * bunch_radius:
        slope_avg += p / r
        slope2_avg += (p / r) ** 2
        count += 1
slope_avg /= count
slope2_avg /= count
slope_err = math.sqrt((slope2_avg - slope_avg * slope_avg))
print("particles slope delta(p)/r            = %12.5g +- %12.5g " % (slope_avg, slope_err))
print("particles slope delta(p)/r from theory= %12.5g" % slope_theory)


nStep = 300
rho_arr = []
for ix in range(nStep + 1):
    x = ix * bunch_radius / nStep
    y = 0.0
    r = math.sqrt(x * x + y * y)
    rho = rhoGrid.getValue(x, y)
    rho_arr.append((r, rho))

phi_arr = []
phi_00 = phiGrid.getValue(0.0, 0.0)
for ix in range(2 * nStep + 1):
    x = (ix - nStep) * bunch_radius / nStep
    y = 0.0
    r = math.sqrt(x * x + y * y)
    phi = phiGrid.getValue(x, y) - phi_00
    phi_arr.append((x, phi))

# -------------------------------------------------
# this is the example of using the Gnuplot package
# -------------------------------------------------
import Gnuplot

gRho = Gnuplot.Gnuplot()
gRho.title("Transverse Density vs. radius")
gRho("set style data lines")
gRho.plot(rho_arr)

gPhi = Gnuplot.Gnuplot()
gPhi.title("Potential vs. distance from center")
gPhi("set style data lines")
gPhi.plot(phi_arr)


input("Please press return to stop:\n")
# -------------------------------------------------


print("Stop.")
