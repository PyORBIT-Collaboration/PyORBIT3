from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT


# Parameters
length = 5.0
length_quad = 0.25 * length
length_drift = 0.25 * length
kq = 0.65
nparts = 4


# Method 1
lattice = TEAPOT_Lattice()

node = DriftTEAPOT("drift1")
node.setLength(length_drift * 0.5)
node.setnParts(nparts)
lattice.addNode(node)

node = QuadTEAPOT("quad1")
node.setLength(length_quad)
node.setnParts(nparts)
node.setParam("kq", +kq)
lattice.addNode(node)

node = DriftTEAPOT("drift2")
node.setLength(length_drift)
node.setnParts(nparts)
lattice.addNode(node)

node = QuadTEAPOT("quad2")
node.setLength(length_quad)
node.setnParts(nparts)
node.setParam("kq", -kq)
lattice.addNode(node)

node = DriftTEAPOT("drift3")
node.setLength(length_drift * 0.5)
node.setnParts(nparts)
lattice.addNode(node)


# Method 2
nodes = [
    DriftTEAPOT(length=length_drift * 0.5, nparts=nparts, name="drift1"),
    QuadTEAPOT(length=length_quad, kq=+kq, nparts=nparts, name="quad1"),
    DriftTEAPOT(length=length_drift, nparts=nparts, name="drift2"),
    QuadTEAPOT(length=length_quad, kq=-kq, nparts=nparts, name="quad2"),
    DriftTEAPOT(length=length_drift * 0.5, nparts=nparts, name="drift3"),
]
lattice = TEAPOT_Lattice()
for node in nodes:
    lattice.addNode(node)
