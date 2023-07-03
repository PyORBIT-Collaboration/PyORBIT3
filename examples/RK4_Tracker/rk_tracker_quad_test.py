# -----------------------------------------------------
# Track ORBIT bunch through the external field
# of the quadrupole Bx = K*y , By = K*x, K in [T/m].
# This example uses the ORBIT Bunch which has the coordinates
# relative to the SyncPart synchronous
# particle instance. The coordinates are (x,x',y,y',z,dE)
# The tracking results are compared to the TEAPOT Quad
# tracking.
# -----------------------------------------------------
import sys
import math
import random
import orbit.core

from bunch import Bunch
from bunch import BunchTwissAnalysis

from trackerrk4 import RungeKuttaTracker
from orbit_utils import PyBaseFieldSource

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import KVDist2D, KVDist3D
from orbit.bunch_generators import GaussDist2D, GaussDist3D
from orbit.bunch_generators import WaterBagDist2D, WaterBagDist3D
from orbit.bunch_generators import TwissAnalysis

from orbit.sns_linac import Quad, Drift

random.seed(100)


# the implementation of the field source
class FieldSource(PyBaseFieldSource):
    """
    Quad. Bx = K*y , By = K*x, K in [T/m]
    """

    def __init__(self, K):
        PyBaseFieldSource.__init__(self)
        self.K = K

    def getElectricMagneticField(self, x, y, z, t):
        # print "Efield x,y,z, t=",(x,y,z,t)
        ex = 0.0
        ey = 0.0
        ez = 0.0
        bx = self.K * y
        by = self.K * x
        bz = 0.0
        return (ex, ey, ez, bx, by, bz)


print("Start.")


b = Bunch()
print("Part. m=", b.mass())
print("Part. q=", b.charge())
TK = 0.1856  # in [GeV]
syncPart = b.getSyncParticle()
syncPart.kinEnergy(TK)

# ---- let's generate particles ---------
alphaX = 3
alphaY = 2
alphaZ = 1

betaX = 1.0
betaY = 2.0
betaZ = 3.0

emittX = 4.0 * 1.0e-6
emittY = 5.0 * 1.0e-6
emittZ = 6.0 * 1.0e-6

print(" aplha beta emitt X=", alphaX, betaX, emittX)
print(" aplha beta emitt Y=", alphaY, betaY, emittY)
print(" aplha beta emitt Z=", alphaZ, betaZ, emittZ)

twissX = TwissContainer(alphaX, betaX, emittX)
twissY = TwissContainer(alphaY, betaY, emittY)
twissZ = TwissContainer(alphaZ, betaZ, emittZ)

distGen = GaussDist3D(twissX, twissY, twissZ)
distGen = WaterBagDist3D(twissX, twissY, twissZ)
distGen = KVDist3D(twissX, twissY, twissZ)
n_parts = 10000
for i in range(n_parts):
    (x, xp, y, yp, z, dE) = distGen.getCoordinates()
    b.addParticle(x, xp, y, yp, z, dE)

b.compress()
syncPart = b.getSyncParticle()
syncPart.kinEnergy(TK)

# copy the initial bunch to another to track through TEAPOT Quad
b1 = Bunch()
b.copyBunchTo(b1)


twiss_analysis = BunchTwissAnalysis()

twiss_analysis.analyzeBunch(b)
print("============before==================")
print("X Twiss =", twiss_analysis.getTwiss(0))
print("Y Twiss =", twiss_analysis.getTwiss(1))
print("Z Twiss =", twiss_analysis.getTwiss(2))
# b.dumpBunch()

G = 30.0  # [T/m]
length = 0.1  # [m]
fieldSource = FieldSource(G)
tracker = RungeKuttaTracker(length)
print("Tracker Entrance plane (a,b,c,d)=", tracker.entrancePlane())
print("Tracker Exit     plane (a,b,c,d)=", tracker.exitPlane())
# the spatial eps is useless because we have quad field = 0 on axis for the syncPart
# tracker.spatialEps(0.00000000001)
# we have to specify the number of time steps
tracker.stepsNumber(40)
tracker.trackBunch(b, fieldSource)

twiss_analysis.analyzeBunch(b)
print("============after RK4 quad tracker===========")
print("Runge-Kutta N steps=", tracker.stepsNumber())
print("X Twiss =", twiss_analysis.getTwiss(0))
print("Y Twiss =", twiss_analysis.getTwiss(1))
print("Z Twiss =", twiss_analysis.getTwiss(2))
# b.dumpBunch()

# make TEAPOT quad
quad = Quad("Test_Quad")
quad.setParam("dB/dr", G)
quad.setLength(length)
quad.trackBunch(b1)

twiss_analysis.analyzeBunch(b1)
print("============after TEAPOT quad=================")
print("X Twiss =", twiss_analysis.getTwiss(0))
print("Y Twiss =", twiss_analysis.getTwiss(1))
print("Z Twiss =", twiss_analysis.getTwiss(2))
# b1.dumpBunch()

print("Stop Calculations.")

"""
phaseXXp = []
for i in range(b1.getSize()):
	x = b1.x(i)
	xp = b1.xp(i)
	phaseXXp.append([x,xp])

#this is the example of using the Gnuplot package
print "Start plot."

import Gnuplot
g = Gnuplot.Gnuplot(debug=0)
g.title('XXp Phase Space')
g('set data style points')
g.plot(phaseXXp)
raw_input('Please press return to stop:\n')
"""
