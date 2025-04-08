import argparse
import copy
import os
import pathlib
from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.danilov_envelope import DanilovEnvelope22
from orbit.danilov_envelope import DanilovEnvelopeTracker22
from orbit.danilov_envelope import DanilovEnvelopeSolverNode22
from orbit.danilov_envelope import add_danilov_envelope_solver_nodes_22
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.space_charge.sc2p5d import SC2p5D_AccNode
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice
from utils import BunchMonitor


# Parse arguments
# --------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--intensity", type=float, default=0.0)
parser.add_argument("--max-part-length", type=float, default=0.1)
parser.add_argument("--mismatch", type=float, default=0.0)
parser.add_argument("--periods", type=int, default=5)
args = parser.parse_args()


# Setup
# --------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Set up simulation
# --------------------------------------------------------------------------------------

envelope = DanilovEnvelope22(
    intrinsic_emittance=20.00e-06,
    eps_x_frac=0.5,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    intensity=args.intensity * 1.0e+14,
    mode=1,
    params=None
)

lattice = make_fodo_lattice(
    phase_adv_x=np.radians(85.0),
    phase_adv_y=np.radians(85.0),
    length=5.0,
    mass=envelope.mass,
    kin_energy=envelope.kin_energy,
    max_part_length=args.max_part_length,
    verbose=1,
)

tracker = DanilovEnvelopeTracker22(lattice, path_length_max=args.max_part_length)
tracker.match_zero_sc(envelope, method="2d")    
envelope_init = envelope.copy()

history = tracker.track(envelope_init.copy(), periods=args.periods, history=True)

figwidth = 4.0 * args.periods
figwidth = min(figwidth, 10.0)

fig, ax = plt.subplots(figsize=(figwidth, 2.5), constrained_layout=True)
ax.plot(history["s"], history["xrms"] * 1000.0)
ax.plot(history["s"], history["yrms"] * 1000.0)
ax.set_ylim(0.0, ax.get_ylim()[1])
ax.set_xlabel("Distance [m]")
ax.set_ylabel("Size [mm]")
plt.show()



