from orbit.core.bunch import Bunch
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice


lattice = TEAPOT_Lattice()
lattice.readMADX("inputs/sns_ring_dual_solenoid.lat", "rnginjsol")

for node in lattice.getNodes():
    print(node)
