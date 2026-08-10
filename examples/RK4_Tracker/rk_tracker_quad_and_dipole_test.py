# -----------------------------------------------------
# Track ORBIT bunch through the external field
# of the dipole with By != 0 and quadrupole Bx = K*y , By = K*x, 
# where K in [T/m].
# This example uses the ORBIT Bunch which has the coordinates
# relative to the SyncPart synchronous particle instance. 
# The coordinates are (x,x',y,y',z,dE) The tracking results are compared to 
# the tracking through TEAPOT thick kickers and quad.
# -----------------------------------------------------
import sys
import math
import random

from orbit.core.bunch import Bunch

from orbit.py_linac.lattice import LinacAccLattice, Sequence
from orbit.py_linac.lattice import ThickKick
from orbit.py_linac.lattice import Drift, Quad, Bend

from orbit.lattice import AccActionsContainer

from orbit.core.orbit_utils import Matrix, PhaseVector

#---- PyORBIT RK4 3D tracking classes
from orbit.core.trackerrk4 import RungeKuttaTracker
from orbit.core.orbit_utils import PyBaseFieldSource

from orbit.core.orbit_utils import FieldSourceContainer
from orbit.core.field_sources import QuadFieldSource
from orbit.core.field_sources import DipoleFieldSource

#--------------------------------------------------------
#---- Start the example
#--------------------------------------------------------

bunch_init = Bunch()
print("Part. m=", bunch_init.mass())
print("Part. q=", bunch_init.charge())
eKin = 1.3  # in [GeV]
bunch_init.getSyncParticle().kinEnergy(eKin)

#---- Add particle to the bunch
(x, xp, y, yp, z, dE) = (0.002 , 0. , 0.003 , 0. , 0. , 0. )
bunch_init.addParticle(x, xp, y, yp, z, dE)

#-----------------------------------------------
#---- Let's make sequence and add acc. nodes
#-----------------------------------------------
sequence = Sequence("Test_Sequence")

dipole_kicker = ThickKick("Dipole")
dipole_length = 0.3
dipole_field = 0.2
dipole_kicker.setLength(dipole_length)
dipole_kicker.setFieldBy(dipole_field)

drift1 = Drift("Drift1")
drift_length = 0.2
drift1.setLength(drift_length)

quad = Quad("Quad")
quad_length = 0.5
quad_field = 5.0
quad.setLength(quad_length)
quad.setParam("dB/dr",quad_field)

drift2 = Drift("Drift2")
drift2.setLength(drift_length)

sequence.addNode(dipole_kicker)
sequence.addNode(drift1)
sequence.addNode(quad)
sequence.addNode(drift2)

#---- Now the lattice with nodes
lattice = LinacAccLattice("test_lattice")
lattice.setNodes(sequence.getNodes())
lattice.initialize()

print ("==========================================================")
print ("Lattice total length = ", lattice.getLength()," [m]")

bunch = Bunch()
bunch_init.copyBunchTo(bunch)

lattice.trackBunch(bunch)
#bunch.dumpBunch()
print ("==========================================================")
print ("TEAPOT tracking final (x,xp) = %+9.5f , %+9.5f ) [mm,mrad]"%(bunch.x(0)*1000,bunch.xp(0)*1000))
print ("TEAPOT tracking final (y,yp) = %+9.5f , %+9.5f ) [mm,mrad]"%(bunch.y(0)*1000,bunch.yp(0)*1000))
print ("==========================================================")

#----------------------------------------------------------
# Now we construct the 3D fields model for the same lattice
#----------------------------------------------------------

dipole_3d = DipoleFieldSource()
dipole_3d.sizesXYZ(10.,10.,dipole_length)
dipole_3d.fieldsXYZ(0.,dipole_field,0.)

quad_3d = QuadFieldSource()
quad_3d.length(quad_length)
quad_3d.gradient(quad_field)

#----------------------------------------------
# Our 3D volume starts with z = 0 plane
# The dipole center position: z = dipole_length/2
# The quad center position: z = dipole_length + drift_length + quad_length/2
# The exit plane z = dipole_length + drift_length + quad_length + drift_length
#----------------------------------------------

#---- Coordinate transformation for dipole_3d
transfCoordsMatrix = Matrix(4,4)
transfCoordsMatrix.unit()
dipole_center_pos = dipole_length/2
transfCoordsMatrix.set(2,3,-dipole_center_pos)
dipole_3d.transormfMatrix(transfCoordsMatrix)


#---- Coordinate transformation for dipole_3d
transfCoordsMatrix = Matrix(4,4)
transfCoordsMatrix.unit()
quad_center_pos = dipole_length + drift_length + quad_length/2
transfCoordsMatrix.set(2,3,-quad_center_pos)
quad_3d.transormfMatrix(transfCoordsMatrix)

#----- field source with all magnets
field_container = FieldSourceContainer()

#---- We have only 2 magnets. There is no need to define empty space (drifts)
field_container.addFieldSource(dipole_3d)
field_container.addFieldSource(quad_3d)

#---------------------------------------------------------
#---- Let's make RK4 tracker
#---------------------------------------------------------

#---- Length definition is used only for apparoximate estimation 
#---- of integration steps for the region from entance to exit planes.
length = dipole_length + drift_length + quad_length + drift_length
tracker = RungeKuttaTracker(length)
#tracker.spatialEps(0.001)
tracker.entrancePlane(0.,0.,-1.,0.)
tracker.exitPlane(0.,0.,1.,-length )
print ("Tracker Entrance plane (a,b,c,d)=",tracker.entrancePlane())
print ("Tracker Exit     plane (a,b,c,d)=",tracker.exitPlane())

bunch = Bunch()
bunch_init.copyBunchTo(bunch)

#---- User can play with the number of time integration steps
n_steps = 40000
tracker.stepsNumber(n_steps)
tracker.trackBunch(bunch,field_container)

synch_part_rvector = bunch.getSyncParticle().rVector()
synch_part_pvector = bunch.getSyncParticle().pVector()
print ("==========================================================")
print ("The synch. particle here has (x,y,z) != (0.,0.,0.) and ")
print ("(xp,yp,pz) != (0.,0.,pz):")
print ("Synch_part. position-vector=",synch_part_rvector)
print ("Synch_part. momentum-vector=",synch_part_pvector)
print ("Angle px/pz [mrad] = ",1000*synch_part_pvector[0]/synch_part_pvector[2])
print ("==========================================================")
print ("The particle coordinates are relative to the synch. particle vector:")
print ("3D Filed tracking final (x,xp) = %+9.3f , %+9.3f ) [mm,mrad]"%(bunch.x(0)*1000,bunch.xp(0)*1000))
print ("3D Filed tracking final (y,yp) = %+9.3f , %+9.3f ) [mm,mrad]"%(bunch.y(0)*1000,bunch.yp(0)*1000))
print ("==========================================================")
print ("After correction to synch. particle with coordinates (0,0,0)")
print ("and momentum vector = (0.,0.,pz), particle coordinates will be:")
x  = bunch.x(0) + synch_part_rvector[0]
xp = bunch.xp(0) + synch_part_pvector[0]/synch_part_pvector[2]
y  = bunch.y(0) + synch_part_rvector[1]
yp = bunch.yp(0) + synch_part_pvector[1]/synch_part_pvector[2]
print ("(x,xp) = %+9.5f , %+9.5f ) [mm,mrad]"%(x*1000,xp*1000))
print ("(y,yp) = %+9.5f , %+9.5f ) [mm,mrad]"%(y*1000,yp*1000))
print ("==========================================================")
print ("Stop.")