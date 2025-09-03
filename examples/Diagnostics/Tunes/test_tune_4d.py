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

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
lattice_params = matrix_lattice.getRingParametersDict()

norm_matrix = None


# Add tune diagnostic node
# ------------------------------------------------------------------------------------

class TuneAnalysisNode(DriftTEAPOT):
    def __init__(self, name: str = "tune_analysis_no_name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.setLength(0.0)
        self._tune_analysis = BunchTuneAnalysis4D()

    def track(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        self._tune_analysis.analyzeBunch(bunch)

    def set_matrix(self, norm_matrix: np.ndarray) -> None:
        self._tune_analysis.setMatrix(norm_matrix.tolist())


tune_node = TuneAnalysisNode()
# tune_node.set_matrix(norm_matrix)
lattice.getNodes()[0].addChildNode(tune_node, 0)


# Generate phase space distribution
# ------------------------------------------------------------------------------------

emittance_x = 25.0e-06
emittance_y = 25.0e-06
bunch_twiss_x = TwissContainer(lattice_params["alpha x"], lattice_params["beta x [m]"], emittance_x)
bunch_twiss_y = TwissContainer(lattice_params["alpha y"], lattice_params["beta y [m]"], emittance_y)
bunch_dist = GaussDist2D(bunch_twiss_x, bunch_twiss_y)

n_parts = 1000
for index in range(n_parts):
    (x, xp, y, yp) = bunch_dist.getCoordinates()
    z = random.uniform(-25.0, 25.0)
    dE = 0.0
    bunch.addParticle(x, xp, y, yp, z, dE)


# Track bunch
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

phase_info = {}
for j, key in enumerate(["phase_x", "phase_y", "tune_x", "tune_y", "action_x", "action_y"]):
    phase_info[key] = []
    for i in range(bunch.getSize()):
        phase_info[key].append(bunch.partAttrValue("ParticlePhaseAttributes", i, j))

phase_info = pd.DataFrame(phase_info)
print(phase_info)