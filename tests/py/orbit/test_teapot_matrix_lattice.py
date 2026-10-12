from pathlib import Path
from pprint import pprint

from orbit.core.bunch import Bunch
from orbit.parsers import MADX_Parser
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.utils.consts import mass_proton

madx_file = Path(__file__).parent / "inputs" / "sns_ring.lat"

lattice = TEAPOT_Lattice()
lattice.readMADX(str(madx_file), "rnginjsol")
lattice.initialize()

bunch = Bunch()
bunch.mass(mass_proton)
bunch.getSyncParticle().kinEnergy(1.3)

matrix_lattice = TEAPOT_MATRIX_Lattice(lattice, bunch)
ring_params = matrix_lattice.getRingParametersDict()

pprint(ring_params)

tune_x, tune_y = matrix_lattice.getTunes()

print("tune_x (from integration):", tune_x)
print("tune_y (from integration):", tune_y)
print("frac_tune_x (from one-turn matrix):", ring_params["fractional tune x"])
print("frac_tune_y (from one-turn matrix):", ring_params["fractional tune y"])
print()
