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
from orbit.core.bunch import BunchTuneAnalysis4D
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


# Setup
# ------------------------------------------------------------------------------------

path = pathlib.Path(__file__)
output_dir = os.path.join("outputs", path.stem)
os.makedirs(output_dir, exist_ok=True)


# Initialize lattice and bunch
# ------------------------------------------------------------------------------------

lattice = make_lattice()

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.000)


# Analyze transfer matrix
# ------------------------------------------------------------------------------------

def build_norm_matrix_from_twiss_2d(alpha: float, beta: float) -> np.ndarray:
    norm_matrix_inv = np.array([[beta, 0.0], [-alpha, 1.0]]) * np.sqrt(1.0 / beta)
    norm_matrix = np.linalg.inv(norm_matrix_inv)
    print(norm_matrix)


matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
lattice_params = matrix_lattice.getRingParametersDict()

lattice_alpha_x = lattice_params["alpha x"]
lattice_alpha_y = lattice_params["alpha y"]
lattice_beta_x = lattice_params["beta x"]
lattice_beta_y = lattice_params["beta y"]

norm_matrix = np.zeros((4, 4))
norm_matrix[0:2, 0:2] = build_norm_matrix_from_twiss_2d(lattice_alpha_x, lattice_beta_x)
norm_matrix[2:4, 2:4] = build_norm_matrix_from_twiss_2d(lattice_alpha_x, lattice_beta_y)


# Add tune diagnostic node
# ------------------------------------------------------------------------------------

class TuneAnalysisNode(DriftTEAPOT):
    def __init__(self, name: str = "tune_analysis_4d_no_name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.setLength(0.0)
        self._tune_analysis = BunchTuneAnalysis4D()

    def track(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        self._tune_analysis.analyzeBunch(bunch)

    def set_norm_matrix(self, norm_matrix: list[list[float]]) -> None:
        norm_matrix_list = list(norm_matrix)
        for i in range(4):
            for j in range(4):
                self._tune_analysis.setMatrixElement(i, j, norm_matrix_list[i][j])

    def get_norm_matrix(self) -> list[list[float]]:
        norm_matrix = [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
        for i in range(4):
            for j in range(4):
                norm_matrix[i][j] = self._tune_analysis.getMatrixElement(i, j)


tune_node = TuneAnalysisNode()
tune_node.set_matrix(norm_matrix)
lattice.getNodes()[0].addChildNode(tune_node, 0)


# Generate phase space distribution
# ------------------------------------------------------------------------------------

bunch_emitt_x = 25.0e-06
bunch_emitt_y = 25.0e-06
bunch_twiss_x = TwissContainer(lattice_alpha_x, lattice_beta_x, bunch_emitt_x)
bunch_twiss_y = TwissContainer(lattice_alpha_y, lattice_beta_y, bunch_emitt_y)
bunch_dist_xy = GaussDist2D(bunch_twiss_x, bunch_twiss_y)

n_parts = 1000
for index in range(n_parts):
    (x, xp, y, yp) = bunch_dist_xy.getCoordinates()
    z = random.uniform(-25.0, 25.0)
    dE = 0.0
    bunch.addParticle(x, xp, y, yp, z, dE)


# Tracking
# ------------------------------------------------------------------------------------

n_turns = 10
for turn in range(n_turns):
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

# Collect phase information from bunch
phase_info = {}
for j, key in enumerate(["phase_1", "phase_2", "tune_1", "tune_2", "action_1", "action_2"]):
    phase_info[key] = []
    for i in range(bunch.getSize()):
        phase_info[key].append(bunch.partAttrValue("ParticlePhaseAttributes", i, j))

phase_info = pd.DataFrame(phase_info)
print(phase_info)

# Read phase information from file
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
        "phase_1",
        "phase_2",
        "tune_1",
        "tune_2",
        "action_1",
        "action_2",
    ] 
)
print(particles.iloc[:, 6:])