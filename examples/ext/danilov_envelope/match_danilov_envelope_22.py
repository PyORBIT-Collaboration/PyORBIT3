import argparse
import copy
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt

from orbit.danilov_envelope import DanilovEnvelope22
from orbit.danilov_envelope import DanilovEnvelopeTracker22
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice


# Parse arguments
# --------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--phase-adv-x", type=float, default=85.0)
parser.add_argument("--phase-adv-y", type=float, default=85.0)
parser.add_argument("--intensity", type=float, default=50.0)
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

envelope = DanilovEnvelope22(
    intrinsic_emittance=20.00e-06,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    intensity=args.intensity * 1.0e+14,
    mode=1,
    params=None
)

lattice = make_fodo_lattice(
    phase_adv_x=np.radians(args.phase_adv_x),
    phase_adv_y=np.radians(args.phase_adv_y),
    length=5.0,
    mass=envelope.mass,
    kin_energy=envelope.kin_energy,
    max_part_length=args.max_part_length,
    verbose=1,
)

tracker = DanilovEnvelopeTracker22(lattice, path_length_max=args.max_part_length)


# Find periodic envelope
# --------------------------------------------------------------------------------------

envelopes = {}

tracker.match_zero_sc(envelope, method="2d")
envelopes["mismatched"] = envelope.copy()

tracker.match(envelope, periods=args.periods, method="replace_avg", verbose=2)
envelopes["matched"] = envelope.copy()


# Plot results bunch
# --------------------------------------------------------------------------------------

histories = {}

for key, envelope in envelopes.items():
    history = tracker.track(envelope, periods=args.periods, history=True)
    histories[key] = copy.deepcopy(history)


figwidth = 3.0 * args.periods
figwidth = min(figwidth, 6.0)

fig, axs = plt.subplots(figsize=(figwidth, 3.0), nrows=2, sharex=True, constrained_layout=True)
for i, key in enumerate(["matched", "mismatched"]):
    alpha = [1.0, 0.15][i]
    axs[0].plot(histories[key]["s"], histories[key]["xrms"] * 1000.0, alpha=alpha, color="blue")
    axs[0].plot(histories[key]["s"], histories[key]["yrms"] * 1000.0, alpha=alpha, color="red")
    axs[1].plot(histories[key]["s"], histories[key]["rxy"], alpha=alpha, color="black")

axs[0].set_ylim(0.0, axs[0].get_ylim()[1])
axs[1].set_ylim(-1.0, 1.0)

for ax in axs:
    ax.set_xlabel("Distance [m]")
axs[0].set_ylabel("Size [mm]")
axs[1].set_ylabel("rxy")

filename = "fig_match.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()