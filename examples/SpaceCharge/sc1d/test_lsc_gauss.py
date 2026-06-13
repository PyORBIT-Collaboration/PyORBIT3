import argparse

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import LSpaceChargeCalc
from orbit.core.spacecharge import SpaceChargeCalcSliceBySlice2D
from orbit.core.spacecharge import Boundary2D
from orbit.bunch_utils import collect_bunch
from orbit.utils.consts import mass_proton

parser = argparse.ArgumentParser()
parser.add_argument("--b-a", type=float, default=3.0)
parser.add_argument("--ring-length", type=float, default=250.0)
parser.add_argument("--nmacros-min", type=int, default=0)
parser.add_argument("--use_sc", type=int, default=1)
parser.add_argument("--nbins", type=int, default=32)
parser.add_argument("--bunch-radius", type=float, default=0.010)
parser.add_argument("--bunch-length-rms", type=float, default=50.0)
parser.add_argument("--nparts", type=int, default=10_000)
parser.add_argument("--intensity", type=float, default=2e14)

parser.add_argument("--grad", type=int, default=0)
parser.add_argument("--modes", type=int, default=None)
args = parser.parse_args()

# Create bunch
bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.0)
bunch.macroSize(args.intensity / args.nparts)

# Generate particles
rng = np.random.default_rng(0)
r = args.bunch_radius * np.sqrt(rng.uniform(0.0, 1.0, size=args.nparts))
t = rng.uniform(0.0, 2.0 * np.pi, size=args.nparts)
x = r * np.cos(t)
y = r * np.sin(t)
xp = np.zeros_like(x)
yp = np.zeros_like(y)
z = rng.normal(scale=args.bunch_length_rms * 0.5, size=args.nparts)
dE = np.zeros_like(z)

coords = np.column_stack([x, xp, y, yp, z, dE])
for i in range(args.nparts):
    bunch.addParticle(*coords[i, :])

# Create space charge calculator
sc_calc = LSpaceChargeCalc(
    args.b_a, args.ring_length, args.nmacros_min, args.use_sc, args.nbins
)
if args.modes:
    sc_calc.setNumModes(args.modes)
if args.grad:
    sc_calc.setUseGrad(args.grad)
    sc_calc.setSmoothGrad(1)

# Track bunch
coords_in = collect_bunch(bunch)["coords"]
sc_calc.trackBunch(bunch)
coords_out = collect_bunch(bunch)["coords"]

# Plot position-energy distribution. The energy gain is proportional
# to the longitudinal electric field.
import matplotlib.pyplot as plt

xmax = [0.5 * args.ring_length, 1.1 * np.max(coords_out[:, 5])]
xmax = np.array(xmax)
limits = list(zip(-xmax, xmax))

fig, ax = plt.subplots()
ax.scatter(coords_out[:, 4], coords_out[:, 5], ec="none", s=1.0, c="black")
ax.set_xlim(limits[0])
ax.set_ylim(limits[1])
plt.show()
