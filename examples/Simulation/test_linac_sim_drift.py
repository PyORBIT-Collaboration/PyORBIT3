import argparse
import math
import os
import pathlib
import random
import sys

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalc3D
from orbit.lattice import AccActionsContainer
from orbit.lattice import AccNode
from orbit.lattice import AccLattice
from orbit.sim.linac import BunchMonitor
from orbit.sim.linac import BunchWriter
from orbit.space_charge.sc3d import setSC3DAccNodes
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.utils.consts import mass_proton


# Setup
# --------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Lattice
# --------------------------------------------------------------------------------------

distance = 5.0
nsteps = 200
delta_s = distance / nsteps

lattice = TEAPOT_Lattice()
node = DriftTEAPOT()
node.setLength(distance)
node.setnParts(nsteps)
lattice.addNode(node)
lattice.initialize()

sc_calc = SpaceChargeCalc3D(64, 64, 64)
sc_nodes = setSC3DAccNodes(lattice, delta_s, sc_calc)


# Bunch
# --------------------------------------------------------------------------------------

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(0.0025)

nparts = 128_000
for i in range(nparts):
    x = 0.010 * random.gauss()
    y = 0.010 * random.gauss()
    z = 0.010 * random.gauss()
    xp = 0.0
    yp = 0.0
    dE = 0.0
    bunch.addParticle(x, xp, y, yp, z, dE)

intensity = 2.00e+09
size_global = bunch.getSizeGlobal()
macro_size = intensity / size_global
bunch.macroSize(macro_size)


# Simulation
# --------------------------------------------------------------------------------------

bunch_writer = BunchWriter(output_dir=output_dir)

bunch_monitor = BunchMonitor(
    output_dir=output_dir,
    stride=0.01,
    stride_write=2.0,
    rf_frequency=402.5e+06,
    bunch_writer=bunch_writer,
    verbose=True,
)


action_container = AccActionsContainer()
action_container.addAction(bunch_monitor, 0) # entrance
action_container.addAction(bunch_monitor, 2) # exit

lattice.trackBunch(bunch, actionContainer=action_container)
