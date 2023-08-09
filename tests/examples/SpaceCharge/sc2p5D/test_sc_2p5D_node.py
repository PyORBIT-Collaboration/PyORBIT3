# -----------------------------------------------------
# Creates Space Charge Calculator Node
# -----------------------------------------------------
import math
import random

random.seed(10)

from orbit.teapot import teapot
from orbit.space_charge.sc2p5d import scLatticeModifications
from orbit.core.bunch import Bunch
from spacecharge import SpaceChargeCalc2p5D, Boundary2D

print("Start.")

# make a bunch
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

# set bunch parameters
macroSize = 1.56e13
energy = 0.08
nParticlesGlobal = b.getSizeGlobal()
b.macroSize(macroSize / nParticlesGlobal)
syncPart = b.getSyncParticle()
syncPart.kinEnergy(energy)

# make a Teapot lattice
elem1 = teapot.DriftTEAPOT("drift1")
elem2 = teapot.QuadTEAPOT("quad1")
elem3 = teapot.QuadTEAPOT("quad2")

teapot_lattice = teapot.TEAPOT_Lattice("teapot_lattice")
teapot_lattice.addNode(elem1)
# teapot_lattice.addNode(elem2)
# teapot_lattice.addNode(elem3)

# set node prameters
elem1.setLength(0.2)
elem1.setnParts(2)
elem2.setLength(1.0)
elem2.setnParts(5)
elem2.addParam("kq", -0.7)
elem3.setLength(2.0)
elem3.setnParts(5)
elem3.addParam("kq", +0.7)

teapot_lattice.initialize()

# make 2.5D space charge calculator
sizeX = 64
sizeY = 64
sizeZ = 20
calc2p5d = SpaceChargeCalc2p5D(sizeX, sizeY, sizeZ)

# set boundary
nBoundaryPoints = 100
N_FreeSpaceModes = 20
R_Boundary = 0.008
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
# boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2*R_Boundary)

sc_path_length_min = 0.05

scNode_arr = scLatticeModifications.setSC2p5DAccNodes(teapot_lattice, sc_path_length_min, calc2p5d)
# scNode_arr = scLatticeModifications.setSC2p5DAccNodes(teapot_lattice, sc_path_length_min, calc2p5d, boundary)

# track lattice
# teapot_lattice.trackBunch(b)
slice_length = 0.1
calc2p5d.trackBunch(b, slice_length)

# -------------------------------------
# momentum change
# coeff is: 2*r0*L*lambda/e*(gamma^3*beta^2)
# r/a^2 and grad/macrosize
# -------------------------------------
xyp_coeff = 2 * b.classicalRadius() * b.charge() ** 2 * 0.1 / (b.getSyncParticle().gamma() ** 3 * b.getSyncParticle().beta() ** 2)
output_string = " "
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
    # 	print("debug xyp_coeff=",xyp_coeff)
    # 	print("r=%10.5g"%r, "x=%10.5g"%x, "y=%10.5g"%y, "xp = %10.5g "%xp," xp_theory = %10.5g "%xp_calc," xp_error = %10.5g "%xp_error,"%")
    # 	print("r=%10.5g"%r, "x=%10.5g"%x, "y=%10.5g"%y, "yp = %10.5g "%yp," yp_theory = %10.5g "%yp_calc," xp_error = %10.5g "%yp_error,"%")

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

    output_string += (
        "debug xyp_coeff= " + xyp_coeff_str + "\n"
        "r="
        + r_str
        + " x="
        + x_str
        + " y="
        + y_str
        + " xp="
        + xp_str
        + " xp_theory="
        + xp_calc_str
        + " xp_error="
        + xp_error_str
        + "%\n"
        + "r="
        + r_str
        + " x="
        + x_str
        + " y="
        + y_str
        + " yp="
        + yp_str
        + " yp_theory="
        + yp_calc_str
        + " xp_error="
        + yp_error_str
        + "%\n"
    )
print(output_string)

print("Finish.")


def test_output():
    expected_output = """ debug xyp_coeff= 1.59072750489e-18
r=0.0037796 x=0.0016331 y=-0.0034085 xp=1.2422e-05 xp_theory=1.235e-05 xp_error=0.57668%
r=0.0037796 x=0.0016331 y=-0.0034085 yp=-2.6098e-05 yp_theory=-2.5777e-05 xp_error=1.2461%
debug xyp_coeff= 1.59072750489e-18
r=0.0022699 x=-0.0020926 y=0.00087946 xp=-1.5318e-05 xp_theory=-1.5825e-05 xp_error=-3.2089%
r=0.0022699 x=-0.0020926 y=0.00087946 yp=6.2988e-06 yp_theory=6.651e-06 xp_error=-5.2956%
debug xyp_coeff= 1.59072750489e-18
r=0.0040419 x=0.0034158 y=0.0021608 xp=2.6262e-05 xp_theory=2.5832e-05 xp_error=1.6651%
r=0.0040419 x=0.0034158 y=0.0021608 yp=1.6457e-05 yp_theory=1.6341e-05 xp_error=0.71075%
debug xyp_coeff= 1.59072750489e-18
r=0.0028626 x=0.0028626 y=5.9773e-08 xp=2.1662e-05 xp_theory=2.1648e-05 xp_error=0.064156%
r=0.0028626 x=0.0028626 y=5.9773e-08 yp=-7.144e-07 yp_theory=4.5204e-10 xp_error=-1.5814e+05%
debug xyp_coeff= 1.59072750489e-18
r=0.0049914 x=0.0013792 y=0.0047971 xp=9.5657e-06 xp_theory=1.043e-05 xp_error=-8.2886%
r=0.0049914 x=0.0013792 y=0.0047971 yp=3.3172e-05 yp_theory=3.6278e-05 xp_error=-8.5623%
debug xyp_coeff= 1.59072750489e-18
r=0.0038833 x=0.0026296 y=-0.0028575 xp=1.9649e-05 xp_theory=1.9886e-05 xp_error=-1.1953%
r=0.0038833 x=0.0026296 y=-0.0028575 yp=-2.1742e-05 yp_theory=-2.161e-05 xp_error=0.6124%
debug xyp_coeff= 1.59072750489e-18
r=0.0041078 x=0.0011006 y=-0.0039576 xp=8.0136e-06 xp_theory=8.3232e-06 xp_error=-3.7193%
r=0.0041078 x=0.0011006 y=-0.0039576 yp=-2.942e-05 yp_theory=-2.993e-05 xp_error=-1.7028%
debug xyp_coeff= 1.59072750489e-18
r=0.0040677 x=0.0030168 y=0.0027286 xp=2.1882e-05 xp_theory=2.2815e-05 xp_error=-4.0866%
r=0.0040677 x=0.0030168 y=0.0027286 yp=2.0659e-05 yp_theory=2.0635e-05 xp_error=0.1156%
debug xyp_coeff= 1.59072750489e-18
r=0.0049558 x=-0.00094734 y=0.0048645 xp=-6.197e-06 xp_theory=-7.1643e-06 xp_error=-13.502%
r=0.0049558 x=-0.00094734 y=0.0048645 yp=3.451e-05 yp_theory=3.6788e-05 xp_error=-6.1921%
debug xyp_coeff= 1.59072750489e-18
r=0.0010519 x=2.6799e-05 y=0.0010516 xp=1.1332e-07 xp_theory=2.0267e-07 xp_error=-44.086%
r=0.0010519 x=2.6799e-05 y=0.0010516 yp=7.5121e-06 yp_theory=7.9525e-06 xp_error=-5.5382%\n"""
    assert output_string == expected_output
