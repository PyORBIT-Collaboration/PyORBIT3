import argparse
import math
import sys
from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.danilov_envelope import DanilovEnvelope20
from orbit.danilov_envelope import DanilovEnvelopeTracker20
from orbit.danilov_envelope import DanilovEnvelopeSolverNode20
from orbit.danilov_envelope import add_danilov_envelope_solver_nodes_20
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.space_charge.sc2p5d import SC2p5Drb_AccNode
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

from utils import make_fodo_lattice
from utils import get_bunch_cov


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--periods", type=int, default=10)
    parser.add_argument("--intensity", type=float, default=0.0)
    parser.add_argument("--eps_x", type=float, default=10.00e-06)
    parser.add_argument("--eps_y", type=float, default=10.00e-06)
    parser.add_argument("--max-part-length", type=float, default=0.1)
    parser.add_argument("--mismatch", type=float, default=0.0)
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    # Create envelope
    envelope = DanilovEnvelope20(
        eps_x=args.eps_x,
        eps_y=args.eps_y,
        mass=mass_proton,
        kin_energy=1.000,
        length=100.0,
        intensity=(args.intensity * 1.00e+14),
        params=None,
    )

    # Create lattice
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
    pprint(lattice_params)

    # Create envelope tracker and fine periodic solution
    tracker = DanilovEnvelopeTracker20(lattice)
    tracker.match_zero_sc(envelope)
    tracker.match(envelope)
    pprint(envelope.twiss())

    # Apply random mismatch
    if args.mismatch:
        envelope.params[0] *= (1.0 + args.mismatch)
        envelope.params[2] *= (1.0 + args.mismatch)
    pprint(envelope.twiss())

    # Track envelope and print turn-by-turn beam size
    data = {}
    for name in ["envelope", "bunch"]:
        data[name] = {}
        for key in ["x_rms", "y_rms", "xp_rms", "yp_rms"]:
            data[name][key] = []

    envelope_copy = envelope.copy()
    for period in range(args.periods + 1):
        if period > 0:
            tracker.track(envelope_copy)
        
        cov_matrix = envelope_copy.cov()
        x_rms = np.sqrt(cov_matrix[0, 0]) * 1000.0
        y_rms = np.sqrt(cov_matrix[2, 2]) * 1000.0
        xp_rms = np.sqrt(cov_matrix[1, 1]) * 1000.0
        yp_rms = np.sqrt(cov_matrix[3, 3]) * 1000.0
        print("{:0.3f} {:0.3f} {:0.3f} {:0.3f}".format(x_rms, y_rms, xp_rms, yp_rms))

        data["envelope"]["x_rms"].append(x_rms)
        data["envelope"]["y_rms"].append(y_rms)
        data["envelope"]["xp_rms"].append(xp_rms)
        data["envelope"]["yp_rms"].append(yp_rms)

    # Track bunch
    lattice = make_fodo_lattice(
        phase_adv_x=np.radians(85.0),
        phase_adv_y=np.radians(85.0),
        length=5.0,
        mass=envelope.mass,
        kin_energy=envelope.kin_energy,
        max_part_length=args.max_part_length,
        verbose=1,
    )
    sc_calc = SpaceChargeCalc2p5D(64, 64, 1)
    sc_path_length_min = 1.00e-06
    sc_nodes = setSC2p5DAccNodes(lattice, sc_path_length_min, sc_calc)

    bunch = envelope.to_bunch(env=False, size=100_000)
    if envelope.intensity > 0.0:
        bunch.macroSize(envelope.intensity / bunch.getSize())

    ############## DEBUG
    cov_matrix_samp = get_bunch_cov(bunch)
    cov_matrix_pred = envelope.cov()
    for cov_matrix in [cov_matrix_pred, cov_matrix_samp]:
        print("debug cov")
        print(np.round(cov_matrix[:4, :4] * 1.00e+06))
    ##############

    for period in range(args.periods + 1):
        if period > 0:
            lattice.trackBunch(bunch)

        cov_matrix = get_bunch_cov(bunch)
        x_rms = np.sqrt(cov_matrix[0, 0]) * 1000.0
        y_rms = np.sqrt(cov_matrix[2, 2]) * 1000.0
        xp_rms = np.sqrt(cov_matrix[1, 1]) * 1000.0
        yp_rms = np.sqrt(cov_matrix[3, 3]) * 1000.0
        print("{:0.3f} {:0.3f} {:0.3f} {:0.3f}".format(x_rms, y_rms, xp_rms, yp_rms))

        data["bunch"]["x_rms"].append(x_rms)
        data["bunch"]["y_rms"].append(y_rms)
        data["bunch"]["xp_rms"].append(xp_rms)
        data["bunch"]["yp_rms"].append(yp_rms)

    # Plot comparison
    fig, axs = plt.subplots(ncols=2, figsize=(7, 3), constrained_layout=True)
    for ax, key in zip(axs, ["x_rms", "y_rms"]):
        ax.plot(data["bunch"][key], color="black")
        ax.plot(data["envelope"][key], marker=".", ms=1.0, lw=0.0, color="red")
    axs[0].set_ylabel("Beam size [mm]")
    for ax in axs:
        ax.set_xlabel("Period")
        ax.set_ylim(0.0, ax.get_ylim()[1] * 2.0)
    plt.show()


if __name__ == "__main__":#
    main(parse_args())