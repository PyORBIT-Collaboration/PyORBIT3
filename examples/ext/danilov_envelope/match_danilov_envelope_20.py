import argparse
import copy
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt

from orbit.danilov_envelope import DanilovEnvelope20
from orbit.danilov_envelope import DanilovEnvelopeTracker20
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice


# Parse arguments
# --------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--intensity", type=float, default=100.0)
parser.add_argument("--eps_x", type=float, default=10.00e-06)
parser.add_argument("--eps_y", type=float, default=10.00e-06)
parser.add_argument("--max-part-length", type=float, default=0.1)
parser.add_argument("--periods", type=int, default=5)
args = parser.parse_args()


# Setup
# --------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Set up simulation
# --------------------------------------------------------------------------------------

envelope = DanilovEnvelope20(
    eps_x=args.eps_x,
    eps_y=args.eps_y,
    mass=mass_proton,
    kin_energy=1.000,
    length=100.0,
    intensity=(args.intensity * 1.00e+14),
    params=None,
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

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, envelope.to_bunch())
lattice_params = matrix_lattice.getRingParametersDict()

tracker = DanilovEnvelopeTracker20(lattice, path_length_max=args.max_part_length)


# Find periodic envelope
# --------------------------------------------------------------------------------------

tracker.match_zero_sc(envelope)

envelope_unmatched = envelope.copy()

tracker.match(envelope, periods=args.periods, verbose=2)
    

# Plot results bunch
# --------------------------------------------------------------------------------------

history_unmatched = tracker.track(envelope_unmatched, periods=args.periods, history=True)
history = tracker.track(envelope, periods=args.periods, history=True)

figwidth = 4.0 * args.periods
figwidth = min(figwidth, 10.0)

fig, ax = plt.subplots(figsize=(figwidth, 2.5), constrained_layout=True)
ax.plot(history["s"], history["xrms"] * 1000.0, color="blue", alpha=1.0)
ax.plot(history["s"], history["yrms"] * 1000.0, color="red",  alpha=1.0)
ax.plot(history_unmatched["s"], history_unmatched["xrms"] * 1000.0, color="blue", alpha=0.2)
ax.plot(history_unmatched["s"], history_unmatched["yrms"] * 1000.0, color="red",  alpha=0.2)
ax.set_ylim(0.0, ax.get_ylim()[1])
ax.set_xlabel("Distance [m]")
ax.set_ylabel("Size [mm]")

filename = "fig_match_rms.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()