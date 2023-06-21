import sys
import math
import orbit.core

import orbit_mpi

from spacecharge import Grid1D
from bunch import Bunch
from orbit_mpi import mpi_comm

print("Start.")

b = Bunch()

nParts = 5
sizeZ = 200
zMin = -5.0
zMax = +5.0
pos_z = 1.0
charge = 1.0

for i in range(nParts):
    b.addParticle(0.1 + i, 0.2 + i, 0.3 + i, 0.4 + i, 0.5 + i, 0.6 + i)
    # b.macroSize(charge)
    b.addPartAttr("macrosize")
    b.partAttrValue("macrosize", i, i, charge)
b.compress()

comm_local = b.getMPIComm()
print("comm name=", orbit_mpi.MPI_Comm_get_name(comm_local))

gridZ = Grid1D(sizeZ)
gridZ.setGridZ(zMin, zMax)
zGrid = gridZ.getGridZ(1)
print("zGrid = ", zGrid)

# ez = gridZ.calcGradient(pos_z)
# print "Ez = ",ez

value = gridZ.getValue(pos_z)
print("value = %5.5f" % value)

# gridZ.binValue(charge,pos_z)
gridZ.binBunch(b)
print("value = ", value)
# gridZ.synchronizeMPI(comm_local)
print("value = ", value)

b.dumpBunch("bunch_dump_test.dat")
print("Stop.")
