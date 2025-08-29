import argparse
import math
import os
import pathlib
import random
from pprint import pprint

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTuneAnalysis
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import GaussDist2D
from orbit.bunch_generators import WaterBagDist2D
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import QuadTEAPOT
from orbit.utils.consts import mass_proton

from utils import make_lattice


def main(args: argparse.Namespace) -> None:

    # Setup
    # ------------------------------------------------------------------------------------

    path = pathlib.Path(__file__)
    output_dir = os.path.join("outputs", path.stem)
    os.makedirs(output_dir, exist_ok=True)


    # Lattice
    # ------------------------------------------------------------------------------------
    
    lattice = make_lattice()

    bunch = Bunch()
    bunch.mass(mass_proton)
    bunch.getSyncParticle().kinEnergy(1.000)

    # Compute lattice parameters from one-turn transfer matrix
    matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
    lattice_params = matrix_lattice.getRingParametersDict()
    pprint(lattice_params)

    # Store some parameters as variables
    lattice_alpha_x = lattice_params["alpha x"]
    lattice_alpha_y = lattice_params["alpha y"]
    lattice_beta_x = lattice_params["beta x [m]"]
    lattice_beta_y = lattice_params["beta y [m]"]
    lattice_eta_x = lattice_params["dispersion x [m]"]
    lattice_etap_x = lattice_params["dispersion x [m]"]


    # Tune diagnostics node
    # ------------------------------------------------------------------------------------

    # The following class is adapted from `orbit.diagnostics` module.

    class TuneAnalysisNode(DriftTEAPOT):
        def __init__(self, name: str = "tuneanalysis no name") -> None:
            DriftTEAPOT.__init__(self, name)
            self.setType("tune calculator teapot")
            self.setLength(0.0)
            self.calc = BunchTuneAnalysis()

        def track(self, params_dict: dict) -> None:
            length = self.getLength(self.getActivePartIndex())
            bunch = params_dict["bunch"]
            self.calc.analyzeBunch(bunch)

        def set_twiss(
            self, 
            beta_x: float,
            alpha_x: float,
            eta_x: float,
            etap_x: float,
            beta_y: float,
            alpha_y: float,
        ) -> None:
            self.calc.assignTwiss(beta_x, alpha_x, eta_x, etap_x, beta_y, alpha_y)


    tune_node = TuneAnalysisNode()
    tune_node.set_twiss(
        beta_x=lattice_beta_x,
        beta_y=lattice_beta_y,
        alpha_x=lattice_alpha_x,
        alpha_y=lattice_alpha_y,
        eta_x=lattice_eta_x,
        etap_x=lattice_etap_x,
    )
    lattice.getNodes()[0].addChildNode(tune_node, 0)


    # Bunch
    # ------------------------------------------------------------------------------------

    # Generate a matched transverse phase space distribution. The longitudinal 
    # distribution will be uniform in position (z) and a delta function in energy 
    # deviation (dE).
    emittance_x = 25.0e-06
    emittance_y = 25.0e-06
    bunch_twiss_x = TwissContainer(lattice_alpha_x, lattice_beta_x, emittance_x)
    bunch_twiss_y = TwissContainer(lattice_alpha_y, lattice_beta_y, emittance_y)
    bunch_dist = GaussDist2D(bunch_twiss_x, bunch_twiss_y)

    for index in range(args.n):
        (x, xp, y, yp) = bunch_dist.getCoordinates()
        z = random.uniform(-25.0, 25.0)
        dE = 0.0
        bunch.addParticle(x, xp, y, yp, z, dE)


    # Tracking
    # ------------------------------------------------------------------------------------

    for turn in range(args.turns):
        lattice.trackBunch(bunch)

        twiss_calc = BunchTwissAnalysis()
        twiss_calc.analyzeBunch(bunch)
        xrms = math.sqrt(twiss_calc.getCorrelation(0, 0)) * 1000.0
        yrms = math.sqrt(twiss_calc.getCorrelation(2, 2)) * 1000.0

        print("turn={} xrms={:0.3f} yrms={:0.3f}".format(turn + 1, xrms, yrms))

    filename = "bunch.dat"
    filename = os.path.join(output_dir, filename)
    bunch.dumpBunch(filename)


    # Analysis    
    # ------------------------------------------------------------------------------------

    particles = np.loadtxt(filename, comments="%")
    particles = pd.DataFrame(
        particles, 
        columns=[  # https://github.com/PyORBIT-Collaboration/PyORBIT3/issues/78
            "x", 
            "xp", 
            "y", 
            "yp", 
            "z",
            "dE",   
            "phase_x",
            "phase_y",
            "tune_x",
            "tune_y",
            "action_x",
            "action_y",
        ] 
    )

    print(particles.iloc[:, 6:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=1000, help="# of macroparticles")
    parser.add_argument("--turns", type=int, default=10)
    args = parser.parse_args()

    main(args)
    