import argparse
import copy
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt

from orbit.envelope import DanilovEnvelope
from orbit.envelope import DanilovEnvelopeMonitor
from orbit.envelope import DanilovEnvelopeTracker
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice

plt.style.use("style.mplstyle")


# Parse arguments
# --------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--phase-adv-x", type=float, default=85.0)
parser.add_argument("--phase-adv-y", type=float, default=85.0)
parser.add_argument("--intensity", type=float, default=50.0e14)
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

envelope = DanilovEnvelope(
    eps_1=20.00e-06,
    eps_2=0.0,
    mass=0.938,
    kin_energy=1.0,
    length=100.0,
    line_density=(args.intensity / 100.0),
    params=None,
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

tracker = DanilovEnvelopeTracker(lattice, path_length_max=args.max_part_length)


# Find periodic envelope
# --------------------------------------------------------------------------------------

envelopes = {}

tracker.match_zero_sc(envelope, method="2d")
envelopes["mismatched"] = envelope.copy()

tracker.match(envelope, periods=args.periods, method="replace_avg", verbose=2)
# tracker.match(
#     envelope, 
#     periods=args.periods, 
#     method="least_squares",
#     xtol=1e-15, 
#     ftol=1e-15, 
#     gtol=1e-15, 
#     verbose=2
# )
envelopes["matched"] = envelope.copy()


# Plot results bunch
# --------------------------------------------------------------------------------------

_, history = tracker.track(envelopes["matched"], periods=args.periods, history=True)
_, history_unmatched = tracker.track(envelopes["mismatched"], periods=args.periods, history=True)


figwidth = 3.0 * args.periods
figwidth = min(figwidth, 10.0)

fig, ax = plt.subplots(figsize=(figwidth, 2.5), constrained_layout=True)
ax.plot(history["s"], history["xrms"] * 1000.0, color="blue", alpha=1.0)
ax.plot(history["s"], history["yrms"] * 1000.0, color="red", alpha=1.0)
ax.plot(history_unmatched["s"], history_unmatched["xrms"] * 1000.0, color="blue", alpha=0.2)
ax.plot(history_unmatched["s"], history_unmatched["yrms"] * 1000.0, color="red", alpha=0.2)
ax.set_ylim(0.0, ax.get_ylim()[1])
ax.set_xlabel("Distance [m]")
ax.set_ylabel("Size [mm]")

filename = "fig_match_rms.png"
filename = os.path.join(output_dir, filename)
plt.savefig(filename, dpi=300)
plt.show()