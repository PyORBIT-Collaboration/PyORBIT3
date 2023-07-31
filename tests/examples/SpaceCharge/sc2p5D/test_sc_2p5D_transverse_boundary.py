# -----------------------------------------------------
# Creates Space Charge Calculator
# test Boundary
# -----------------------------------------------------
import sys
import math
import random
import pytest

random.seed(10)
from orbit.core import orbit_mpi

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalc2p5D, Boundary2D

print("Start.")

# Make a SC solver
sizeX = 64
sizeY = 64
sizeZ = 20
calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)

# Make bunch
b = Bunch()
bunch_radius = 0.005
bunch_length = 500e-9 * 3e8 * 0.87502565
nParts = 100000

# Use random radius
for ip in range(nParts):
    r = bunch_radius * math.sqrt(random.random())
    phi = 2 * math.pi * random.random()
    x = r * math.sin(phi)
    y = r * math.cos(phi)
    z = bunch_length * 0.5 * (1.0 - 2 * random.random())
    b.addParticle(x, 0.0, y, 0.0, z, 0.0)

macroSize = 1.56e13
energy = 0.08
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(macroSize / nParticlesGlobal)
b.getSyncParticle().kinEnergy(energy)

# bunchExtremaCalc->getExtremaXYZ

# Make a boundary
nBoundaryPoints = 100
N_FreeSpaceModes = 20
R_Boundary = 0.005
boundary = Boundary2D(nBoundaryPoints, N_FreeSpaceModes)

for i in range(nBoundaryPoints):
    x = R_Boundary * math.cos((2.0 * math.pi / (nBoundaryPoints - 1)) * i)
    y = R_Boundary * math.sin((2.0 * math.pi / (nBoundaryPoints - 1)) * i)
    boundary.setBoundaryPoint(i, x, y)
boundary.initialize()
b_maxx = boundary.getMaxX()
b_minx = boundary.getMinX()
b_maxy = boundary.getMaxY()
b_miny = boundary.getMinY()

result_string = f"MaxX={b_maxx:.11g} MinX={b_minx:.12g} MaxY={b_maxy:.12g} MinY={b_miny:.12g}"
print(result_string)

# -------------------------------------------------------------------------
"""
For the Boundary2D class it is possible to have the predefined shapes
Circle - the last parameter will be the diameter of the circle
Ellipse - there will be two parameters at the end - 2*a and 2*b where a,b are semi-axises
Rectangle - there will be two parameters at the end - horizontal and vertical sizes
"""
# -------------------------------------------------------------------------
# boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2*R_Boundary)


# Shape boundary test
# Choose from Circle/Ellipse/Rectangle
shapeboundary = Boundary2D(nBoundaryPoints, N_FreeSpaceModes, "Circle", 0.01, 0.01)
print("shape name=", boundary.getShapeName())

# Set SC node parameters
pipe_radius = 0.02
slice_length = 0.1

# Track and analysis
# analysis
rhoGrid = calc2p5d.getRhoGrid()
phiGrid = calc2p5d.getPhiGrid()
x = bunch_radius / 2.0
y = 0.0
r = math.sqrt(x * x + y * y)

rho_theory = macroSize * 4.0 / (math.pi * rhoGrid.getSizeX() * rhoGrid.getSizeY())
phi_theory = macroSize * r**2 / (bunch_radius**2)
grad_theory = macroSize * 2 * r / (bunch_radius**2)

# without boundary
# print "without boundary."
# calc2p5d.trackBunch(b,slice_length)

# with boundary
print("with boundary.")
# calc2p5d.trackBunch(b,slice_length,boundary)
calc2p5d.trackBunch(b, slice_length, shapeboundary)

rho = rhoGrid.getValue(x, y)
phi = 2 * (phiGrid.getValue(x, y) - phiGrid.getValue(0.0, 0.0))
(ex, ey) = phiGrid.calcGradient(x, y)
grad = 2 * math.sqrt(ex * ex + ey * ey)

with_boundary = (
    f"r={r} rho = {rho:12.5g} rho_theory = {rho_theory:12.5g}\n"
    f"r={r} phi = {phi:12.5g} phi_theory = {phi_theory:12.5g}\n"
    f"r={r} grad = {grad:12.5g} grad_theory = {grad_theory:12.5g}"
)

# Print the result
print(with_boundary)
output_string = " "
# -------------------------------------
# momentum change
# coeff is: 2*r0*L*lambda/e*(gamma^3*beta^2)
# r/a^2 and grad/macrosize
# -------------------------------------
xyp_coeff = 2 * b.classicalRadius() * b.charge() ** 2 * slice_length / (b.getSyncParticle().gamma() ** 3 * b.getSyncParticle().beta() ** 2)
for ip in range(10):
    x = b.x(ip)
    y = b.y(ip)
    r = math.sqrt(x * x + y * y)
    theta = math.atan2(y, x)
    xp = b.xp(ip)
    yp = b.yp(ip)
    xp_calc = math.cos(theta) * r * xyp_coeff * macroSize / (bunch_length * bunch_radius**2)
    yp_calc = math.sin(theta) * r * xyp_coeff * macroSize / (bunch_length * bunch_radius**2)
    xp_error = (xp - xp_calc) * 100 / xp_calc
    yp_error = (yp - yp_calc) * 100 / yp_calc
    # 	print("debug xyp_coeff=",xyp_coeff +"\n")
    # 	print("r=%10.5g"%r, "x=%10.5g"%x, "y=%10.5g"%y, "xp = %10.5g "%xp," xp_theory = %10.5g "%xp_calc," xp_error = %10.5g "%xp_error,"%")
    # 	print("r=%10.5g"%r, "x=%10.5g"%x, "y=%10.5g"%y, "yp = %10.5g "%yp," yp_theory = %10.5g "%yp_calc," xp_error = %10.5g "%yp_error,"%")
    # Round the floating-point values to 5 decimal places
    xyp_coeff_str = "{:.12g}".format(xyp_coeff)
    r_str = "{:.5g}".format(r)
    x_str = "{:.5g}".format(x)
    y_str = "{:.5g}".format(y)
    xp_str = "{:.5g}".format(xp)
    xp_calc_str = "{:.5g}".format(xp_calc)
    xp_error_str = "{:.5g}".format(xp_error)
    yp_str = "{:.5g}".format(yp)
    yp_calc_str = "{:.5g}".format(yp_calc)
    yp_error_str = "{:.5g}".format(yp_error)

    # Create the output string
    output_string += (
        "debug xyp_coeff= " + xyp_coeff_str + "\n"
        "r= "
        + r_str
        + " x= "
        + x_str
        + " y= "
        + y_str
        + " xp= "
        + xp_str
        + " xp_theory= "
        + xp_calc_str
        + " xp_error= "
        + xp_error_str
        + "%\n"
        + "r= "
        + r_str
        + " x= "
        + x_str
        + " y= "
        + y_str
        + " yp= "
        + yp_str
        + " yp_theory= "
        + yp_calc_str
        + " xp_error= "
        + yp_error_str
        + "%\n"
    )

# Print the result
print(output_string)


# pytest
def test_boundary():
    expected_boundary_values = "MaxX=0.005 MinX=-0.00499748271192 MaxY=0.00499937063837 MinY=-0.00499937063837"

    assert result_string == expected_boundary_values


def test_smth():
    expected_output = """r=0.0025 rho =    5.254e+09 rho_theory =   4.8493e+09
r=0.0025 phi =  -3.8368e+12 phi_theory =      3.9e+12
r=0.0025 grad =   3.1468e+15 grad_theory =     3.12e+15"""
    assert with_boundary == expected_output


def test_smthelse():
    expected = """ debug xyp_coeff= 1.59072750489e-18
r= 0.0037796 x= 0.0016331 y= -0.0034085 xp= 1.2341e-05 xp_theory= 1.235e-05 xp_error= -0.079547%
r= 0.0037796 x= 0.0016331 y= -0.0034085 yp= -2.5848e-05 yp_theory= -2.5777e-05 xp_error= 0.2749%
debug xyp_coeff= 1.59072750489e-18
r= 0.0022699 x= -0.0020926 y= 0.00087946 xp= -1.536e-05 xp_theory= -1.5825e-05 xp_error= -2.9424%
r= 0.0022699 x= -0.0020926 y= 0.00087946 yp= 6.4421e-06 yp_theory= 6.651e-06 xp_error= -3.1405%
debug xyp_coeff= 1.59072750489e-18
r= 0.0040419 x= 0.0034158 y= 0.0021608 xp= 2.6372e-05 xp_theory= 2.5832e-05 xp_error= 2.0908%
r= 0.0040419 x= 0.0034158 y= 0.0021608 yp= 1.6533e-05 yp_theory= 1.6341e-05 xp_error= 1.172%
debug xyp_coeff= 1.59072750489e-18
r= 0.0028626 x= 0.0028626 y= 5.9773e-08 xp= 2.1731e-05 xp_theory= 2.1648e-05 xp_error= 0.38167%
r= 0.0028626 x= 0.0028626 y= 5.9773e-08 yp= -5.5562e-07 yp_theory= 4.5204e-10 xp_error= -1.2301e+05%
debug xyp_coeff= 1.59072750489e-18
r= 0.0049914 x= 0.0013792 y= 0.0047971 xp= 9.4637e-06 xp_theory= 1.043e-05 xp_error= -9.2658%
r= 0.0049914 x= 0.0013792 y= 0.0047971 yp= 3.2994e-05 yp_theory= 3.6278e-05 xp_error= -9.0514%
debug xyp_coeff= 1.59072750489e-18
r= 0.0038833 x= 0.0026296 y= -0.0028575 xp= 1.9623e-05 xp_theory= 1.9886e-05 xp_error= -1.3258%
r= 0.0038833 x= 0.0026296 y= -0.0028575 yp= -2.1444e-05 yp_theory= -2.161e-05 xp_error= -0.76728%
debug xyp_coeff= 1.59072750489e-18
r= 0.0041078 x= 0.0011006 y= -0.0039576 xp= 7.8835e-06 xp_theory= 8.3232e-06 xp_error= -5.2833%
r= 0.0041078 x= 0.0011006 y= -0.0039576 yp= -2.9201e-05 yp_theory= -2.993e-05 xp_error= -2.434%
debug xyp_coeff= 1.59072750489e-18
r= 0.0040677 x= 0.0030168 y= 0.0027286 xp= 2.1983e-05 xp_theory= 2.2815e-05 xp_error= -3.6476%
r= 0.0040677 x= 0.0030168 y= 0.0027286 yp= 2.0707e-05 yp_theory= 2.0635e-05 xp_error= 0.34899%
debug xyp_coeff= 1.59072750489e-18
r= 0.0049558 x= -0.00094734 y= 0.0048645 xp= -6.4238e-06 xp_theory= -7.1643e-06 xp_error= -10.337%
r= 0.0049558 x= -0.00094734 y= 0.0048645 yp= 3.4551e-05 yp_theory= 3.6788e-05 xp_error= -6.0786%
debug xyp_coeff= 1.59072750489e-18
r= 0.0010519 x= 2.6799e-05 y= 0.0010516 xp= 9.3747e-08 xp_theory= 2.0267e-07 xp_error= -53.744%
r= 0.0010519 x= 2.6799e-05 y= 0.0010516 yp= 7.6252e-06 yp_theory= 7.9525e-06 xp_error= -4.1166%\n"""
    assert output_string == expected
