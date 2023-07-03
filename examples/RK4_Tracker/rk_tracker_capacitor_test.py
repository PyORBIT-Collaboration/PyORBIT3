# -----------------------------------------------------
# Track the ORBIT bunch through the external field
# of the capacitor Ex , Ey , Ez = E0 , all in [V/m].
# To avoid to act on the synch particle the if-conditions is
# used inside the field source.
# This example uses the ORBIT Bunch which has the coordinates
# relative to the SyncPart synchronous
# particle instance. The coordinates are (x,x',y,y',z,dE)
# The tracking results are compared to the nalytical calculations.
# -----------------------------------------------------
import sys
import math
import random
import orbit.core

from bunch import Bunch

from trackerrk4 import RungeKuttaTracker
from orbit_utils import PyBaseFieldSource

random.seed(100)


# the implementation of the field source
class FieldSource(PyBaseFieldSource):
    """
    Quad. Bx = K*y , By = K*x, K in [T/m]
    """

    def __init__(self, Ex, Ey, E0):
        PyBaseFieldSource.__init__(self)
        self.Ex = Ex
        self.Ey = Ey
        self.E0 = E0

    def getElectricMagneticField(self, x, y, z, t):
        # print "Efield x,y,z, t=",(x,y,z,t)
        ex = Ex
        ey = Ey
        ez = E0
        if x**2 + y**2 < (0.00001) ** 2:
            ex = 0.0
            ey = 0.0
            ez = 0.0
        bx = 0.0
        by = 0.0
        bz = 0.0
        return (ex, ey, ez, bx, by, bz)


print("Start.")


b_init = Bunch()
print("Part. m=", b_init.mass())
print("Part. q=", b_init.charge())
TK = 0.1856  # in [GeV]
syncPart = b_init.getSyncParticle()
syncPart.kinEnergy(TK)

b_init.addParticle(0.001, 0.0, 0.0, 0.0, 0.0, 0.0)
b_init.compress()

print("============before================================================================")
(x, xp, y, yp, z, dE) = (b_init.x(0), b_init.xp(0), b_init.y(0), b_init.yp(0), b_init.z(0), b_init.dE(0))
print("(x,xp,y,yp,z,dE)= (%12.9f,%12.9f,%12.9f,%12.9f,%12.9f,%12.9f)" % (x, xp, y, yp, z, dE))

# -------------------------------------------
# z- acceleration
# -------------------------------------------
b = Bunch()
b_init.copyBunchTo(b)

Ex = 0.0
Ey = 0.0
E0 = 100000.0  # [V/m]
length = 1.0  # [m]
fieldSource = FieldSource(Ex, Ey, E0)
tracker = RungeKuttaTracker(length)
print("Tracker Entrance plane (a,b,c,d)=", tracker.entrancePlane())
print("Tracker Exit     plane (a,b,c,d)=", tracker.exitPlane())
tracker.spatialEps(0.0000001)
# we can to specify the number of time steps
tracker.stepsNumber(40)

tracker.trackBunch(b, fieldSource)

print("============after RK4 Ez non-zero tracker===========================================")
(x, xp, y, yp, z, dE) = (b.x(0), b.xp(0), b.y(0), b.yp(0), b.z(0), b.dE(0))
print("(x,xp,y,yp,z,dE)= (%12.9f,%12.9f,%12.9f,%12.9f,%12.9f,%12.9f)" % (x, xp, y, yp, z, dE))
print("theory said we should have dE [GeV] =", E0 * length / 1.0e9)
# b.dumpBunch()

# -------------------------------------------
# x- acceleration
# -------------------------------------------
b = Bunch()
b_init.copyBunchTo(b)

Ex = 100000.0
Ey = 0.0
E0 = 0.0  # [V/m]
length = 1.0  # [m]
fieldSource = FieldSource(Ex, Ey, E0)
tracker = RungeKuttaTracker(length)
tracker.spatialEps(0.0000001)
# we can to specify the number of time steps
tracker.stepsNumber(40)

tracker.trackBunch(b, fieldSource)

print("============after RK4 Ex non-zero tracker===========================================")
(x, xp, y, yp, z, dE) = (b.x(0), b.xp(0), b.y(0), b.yp(0), b.z(0), b.dE(0))

print("(x,xp,y,yp,z,dE)= (%12.9f,%12.9f,%12.9f,%12.9f,%12.9f,%12.9f)" % (x, xp, y, yp, z, dE))
c = 2.99792458e8
d_time = length / (c * syncPart.beta())
d_px = d_time * Ex
xp = c * d_px * 1.0e-9 / syncPart.momentum()

print("theory said we should have xp [rad] =", xp)
# b.dumpBunch()

# -------------------------------------------
# y- acceleration
# -------------------------------------------
b = Bunch()
b_init.copyBunchTo(b)

Ex = 0.0
Ey = 100000.0
E0 = 0.0  # [V/m]
length = 1.0  # [m]
fieldSource = FieldSource(Ex, Ey, E0)
tracker = RungeKuttaTracker(length)
tracker.spatialEps(0.0000001)
# we can to specify the number of time steps
tracker.stepsNumber(40)

tracker.trackBunch(b, fieldSource)

print("============after RK4 Ey non-zero tracker===========================================")
(x, xp, y, yp, z, dE) = (b.x(0), b.xp(0), b.y(0), b.yp(0), b.z(0), b.dE(0))

print("(x,xp,y,yp,z,dE)= (%12.9f,%12.9f,%12.9f,%12.9f,%12.9f,%12.9f)" % (x, xp, y, yp, z, dE))
c = 2.99792458e8
d_time = length / (c * syncPart.beta())
d_py = d_time * Ey
yp = c * d_py * 1.0e-9 / syncPart.momentum()

print("theory said we should have yp [rad] =", yp)
# b.dumpBunch()

print("=====================================================================================")
print("Stop Calculations.")
