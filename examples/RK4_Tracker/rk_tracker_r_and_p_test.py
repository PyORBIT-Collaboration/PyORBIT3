# -----------------------------------------------------
# Track bunch with r and p through the external field
# The field is 1 T and has direction (0,1,0).
# This example uses the Bunch as a container of
# 6D absolute coordinates without SyncPart synchronous
# particle instance. The coordinates are (x,Px,y,Py,z,Pz)
# -----------------------------------------------------
import sys
import math
import orbit.core

from bunch import Bunch

from trackerrk4 import RungeKuttaTracker
from trackerrk4 import PyExternalEffects
from orbit_utils import PyBaseFieldSource


# the implementation of the field source
class FieldSource(PyBaseFieldSource):
    def __init__(self):
        PyBaseFieldSource.__init__(self)

    def getElectricMagneticField(self, x, y, z, t):
        # print "Efield x,y,z, t=",(x,y,z,t)
        ex = 0.0
        ey = 0.0
        ez = 0.0
        bx = 0.0
        by = 1.0
        bz = 0.0
        return (ex, ey, ez, bx, by, bz)


class ExternalEffects(PyExternalEffects):
    def __init__(self):
        PyExternalEffects.__init__(self)
        self.count = 0
        self.trajectoryXYArr = []

    def setupEffects(self, bunch):
        pass

    def prepareEffects(self, bunch, t):
        pass

    def finalizeEffects(self, bunch):
        pass

    def applyEffects(self, bunch, t, t_step, field_source, tracker):
        (x, y, z) = (bunch.x(0), bunch.y(0), bunch.z(0))
        self.trajectoryXYArr.append([x, z])
        print(" %5d " % self.count, " %8.6f  %8.6f  %8.6f " % (x, y, z))
        self.count = self.count + 1
        pass

    def applyEffectsForEach(self, bunch, index, inVct, outVct, t, t_step, field_source, tracker):
        print("  macro-part. index = %5d " % index)
        self.count = self.count + 1
        pass


print("Start.")


b = Bunch()
print("Part. m=", b.mass())
print("Part. q=", b.charge())

TK = 1.0
E = b.mass() + TK
P = math.sqrt(E * E - b.mass() * b.mass())
c_light = 2.99792458e8

print("TK[GeV] = ", TK)
print("P[GeV/c] = ", P)

b.addParticle(0.0, 0.0, 0.0, 0.0, 0.0, P)
b.compress()

# radius estimation
R = P * 1.0e9 / (c_light * b.charge() * 1.0)
print("R[m] = ", R)
n_step = 1000
time_step = (2 * 3.1415926 * R / (c_light * P / E)) / n_step

fS = FieldSource()
extEff = ExternalEffects()
print("ExternalEffects name=", extEff.name())

tracker = RungeKuttaTracker(10000.0)
print("Entrance plane (a,b,c,d)=", tracker.entrancePlane())
print("Exit     plane (a,b,c,d)=", tracker.exitPlane())
print("Length[m]=", tracker.length())

print("Start tracking.")
time_start = 0.0
time_max = time_step * n_step
print("================================================")
print("Step_Index    x     y   z ")
tracker.track(b, time_start, time_max, time_step, fS, extEff)
print("================================================")
print("Stop tracking.")

print("time step=", tracker.timeStep())
print("Stop.")

"""
#this is the example of using the Gnuplot package
import Gnuplot
g = Gnuplot.Gnuplot(debug=1)
g.title('XZ trajectory')
g('set data style linespoints')
g.plot(extEff.trajectoryXYArr)
raw_input('Please press return to stop:\n')
"""
