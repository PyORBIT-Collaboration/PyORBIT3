"""Test envelope tracker speed in SNS ring."""

import argparse
import time
import cProfile
import pstats
import sys

import numpy as np
from tqdm import trange

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.envelope import Envelope
from orbit.envelope import EnvelopeTracker
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

sys.path.append("..")
from utils import gen_dist

parser = argparse.ArgumentParser()
parser.add_argument("--bunch-length", type=float, default=120.0)
parser.add_argument("--kin-energy", type=float, default=1.300)
parser.add_argument("--intensity", type=float, default=2e14)

parser.add_argument("--nparts", type=int, default=10_000)
parser.add_argument("--turns", type=int, default=25)
parser.add_argument("--sc", type=int, default=0)
parser.add_argument("--sc-grid", type=int, default=64)
args = parser.parse_args()


lattice = TEAPOT_Ring()
lattice.readMADX("inputs/sns_ring.lat", "rnginjsol")
lattice.initialize()

for node in lattice.getNodes():
    try:
        node.setUsageFringeFieldIN(False)
        node.setUsageFringeFieldOUT(False)
    except:
        pass

for node in lattice.getNodes():
    max_length = 1.0
    if node.getLength() > max_length:
        node.setnParts(1 + int(node.getLength() / max_length))

bunch = Bunch()
bunch.mass(mass_proton)
sync_part = bunch.getSyncParticle()
sync_part.kinEnergy(args.kin_energy)

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
matrix_lattice_params = matrix_lattice.getRingParametersDict()
alpha_x = matrix_lattice_params["alpha x"]
alpha_y = matrix_lattice_params["alpha y"]
beta_x = matrix_lattice_params["beta x [m]"]
beta_y = matrix_lattice_params["beta y [m]"]
eps_x = 25.0e-06
eps_y = eps_x

cov_matrix = np.zeros((6, 6))
cov_matrix[0, 0] = eps_x * beta_x
cov_matrix[2, 2] = eps_y * beta_y
cov_matrix[0, 1] = cov_matrix[1, 0] = -eps_x * alpha_x
cov_matrix[2, 3] = cov_matrix[3, 2] = -eps_y * alpha_y
cov_matrix[1, 1] = eps_x * (1.0 + alpha_x**2) / beta_x
cov_matrix[3, 3] = eps_y * (1.0 + alpha_y**2) / beta_y
cov_matrix[4, 4] = (args.bunch_length / 4.0) ** 2
cov_matrix[5, 5] = 0.0

cov_matrix_init = np.copy(cov_matrix)

# Track envelope
print("ENVELOPE")

envelope = Envelope(
    bunch=bunch,
    cov_matrix=cov_matrix_init,
    intensity=args.intensity,
)
envelope_tracker = EnvelopeTracker(lattice, space_charge=("2d" if args.sc else None))

start_time = time.time()

profiler = cProfile.Profile()
profiler.enable()

for turn in trange(args.turns):
    envelope_tracker.track(envelope)

time_per_turn = (time.time() - start_time) / args.turns

profiler.disable()

print("Time per turn:", time_per_turn)

profiler_stats = pstats.Stats(profiler)
profiler_stats.sort_stats(pstats.SortKey.TIME)
profiler_stats.print_stats(10)

# Track bunch
print("BUNCH")

rng = np.random.default_rng()

bunch_coords = np.zeros((args.nparts, 6))
bunch_coords[:, :4] = gen_dist(
    size=args.nparts, cov_matrix=cov_matrix_init[0:4, 0:4], name="kv"
)
bunch_coords[:, 4] = args.bunch_length * rng.uniform(-0.5, 0.5, size=args.nparts)

for i in range(bunch_coords.shape[0]):
    bunch.addParticle(*bunch_coords[i])

if args.sc:
    sc_calc = SpaceChargeCalc2p5D(64, 64, 1)
    sc_path_length_min = 1.00e-06
    sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

    bunch_size = bunch.getSizeGlobal()
    bunch.macroSize(args.intensity / bunch_size)

start_time = time.time()

profiler = cProfile.Profile()
profiler.enable()

for turn in trange(args.turns):
    lattice.trackBunch(bunch)

time_per_turn = (time.time() - start_time) / args.turns

profiler.disable()

print("Time per turn:", time_per_turn)

profiler_stats = pstats.Stats(profiler)
profiler_stats.sort_stats(pstats.SortKey.TIME)
profiler_stats.print_stats(10)
