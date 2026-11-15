"""Profile envelope tracking in FODO lattice."""

import argparse
import time
import cProfile
import pstats

import numpy as np
from tqdm import trange

from orbit.core.bunch import Bunch
from orbit.envelope import Envelope
from orbit.teapot import QuadTEAPOT
from orbit.teapot import DriftTEAPOT
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton


parser = argparse.ArgumentParser()
parser.add_argument("--bunch-length", type=float, default=5.0)
parser.add_argument("--kin-energy", type=float, default=1.300)
parser.add_argument("--intensity", type=float, default=2e14)

parser.add_argument("--nslice", type=int, default=10)
parser.add_argument("--kq", type=float, default=0.25)

parser.add_argument("--turns", type=int, default=500)
parser.add_argument("--sc", type=int, default=0)
parser.add_argument("--fit", type=int, default=1)
parser.add_argument("--static", type=int, default=1)
args = parser.parse_args()

nodes = [
    QuadTEAPOT(length=0.5, kq=+args.kq),
    DriftTEAPOT(length=1.0),
    QuadTEAPOT(length=1.0, kq=-args.kq),
    DriftTEAPOT(length=1.0),
    QuadTEAPOT(length=0.5, kq=+args.kq),
]

lattice = TEAPOT_Lattice()
for node in nodes:
    node.setnParts(args.nslice)
    node.setUsageFringeFieldIN(False)
    node.setUsageFringeFieldOUT(False)
    lattice.addNode(node)
lattice.initialize()

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
cov_matrix[4, 4] = args.bunch_length**2 / 12.0
cov_matrix[5, 5] = 0.0

envelope = Envelope(
    cov_matrix=cov_matrix,
    kin_energy=sync_part.kinEnergy(),
    mass=bunch.mass(),
    charge=bunch.charge(),
    intensity=args.intensity,
)
sc = "2d" if args.sc else None

profiler = cProfile.Profile()
profiler.enable()
for turn in trange(args.turns):
    lattice.trackEnvelope(envelope, sc=sc, fit=args.fit, static=args.static)
profiler.disable()

stats = pstats.Stats(profiler)
stats.sort_stats(pstats.SortKey.TIME)
stats.print_stats(10)
