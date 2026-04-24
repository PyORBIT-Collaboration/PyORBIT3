import math
import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import SyncParticle
from orbit.envelope import Envelope
from orbit.utils.consts import mass_proton


mass = mass_proton  # [GeV]
kin_energy = 1.000  # [GeV]

bunch = Bunch()
bunch.mass(mass)

sync_part = bunch.getSyncParticle()
sync_part.kinEnergy(kin_energy)

cov_matrix = np.identity(6)
cov_matrix[0, 0] = 10.0e-3 ** 2
cov_matrix[2, 2] = 10.0e-3 ** 2
cov_matrix[4, 4] = 50.0

envelope = Envelope(sync_part=sync_part, cov_matrix=cov_matrix)


