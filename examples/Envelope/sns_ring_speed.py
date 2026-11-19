"""Profile envelope tracking in SNS ring lattice."""

import argparse
import time
import cProfile
import pstats
import sys

import numpy as np
from tqdm import trange

from orbit.core.bunch import Bunch
from orbit.envelope import Envelope
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import BendTEAPOT
from orbit.utils.consts import mass_proton


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bunch-length", type=float, default=120.0)
    parser.add_argument("--kin-energy", type=float, default=1.300)
    parser.add_argument("--intensity", type=float, default=2e14)

    parser.add_argument("--nparts", type=int, default=25_000)
    parser.add_argument("--turns", type=int, default=25)
    parser.add_argument("--sc", type=int, default=0)
    parser.add_argument("--sc-grid", type=int, default=64)
    parser.add_argument("--fringe", type=int, default=0)
    parser.add_argument("--fit", type=int, default=0)
    parser.add_argument("--static", type=int, default=1)
    parser.add_argument("--fallback", type=int, default=0)
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
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

    envelope = Envelope(
        kin_energy=sync_part.kinEnergy(),
        mass=bunch.mass(),
        charge=bunch.charge(),
        time=sync_part.time(),
        cov_matrix=cov_matrix_init,
        intensity=args.intensity,
    )

    profiler = cProfile.Profile()
    profiler.enable()
    for turn in trange(args.turns):
        lattice.trackEnvelope(
            envelope,
            static=args.static,
            fit=args.fit,
            fallback=args.fallback,
            sc=("2d" if args.sc else None),
        )
    profiler.disable()

    profiler_stats = pstats.Stats(profiler)
    profiler_stats.sort_stats(pstats.SortKey.TIME)
    profiler_stats.print_stats(10)


if __name__ == "__main__":
    main(parse_args())
