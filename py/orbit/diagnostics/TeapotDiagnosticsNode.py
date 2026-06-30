"""TEAPOT-style bunch diagnostic nodes."""
from typing import IO

import numpy as np

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTuneAnalysis
from orbit.teapot import DriftTEAPOT

from .diagnostics import StatLats
from .diagnostics import StatLatsSetMember
from .diagnostics import Moments
from .diagnostics import MomentsSetMember
from .diagnostics import BPMSignal
from .matrix import build_norm_matrix_from_cov
from .matrix import build_norm_matrix_from_tmat
from .matrix import TransferMatrixAnalysis


class TeapotStatLatsNode(DriftTEAPOT):
    def __init__(self, filename: str, name: str = "statlats no name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.statlats = StatLats(filename)
        self.setType("statlats teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattice_length = 0.0
        self.file_out = open(filename, "w")

    def track(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        self.statlats.writeStatLats(self.position, bunch, self.lattice_length)

    def setPosition(self, position: float) -> None:
        self.position = position

    def closeStatLats(self) -> None:
        self.file_out.close()

    def setLatticeLength(self, length: float) -> None:
        self.lattice_length = length


class TeapotStatLatsNodeSetMember(DriftTEAPOT):
    def __init__(self, file: IO[str], name: str = "statlats no name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.statlats = StatLatsSetMember(file)
        self.setType("statlats teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattice_length = 0.0
        self.active = True
        self.file = file

    def track(self, params_dict: dict) -> None:
        if self.active:
            bunch = params_dict["bunch"]
            self.statlats.writeStatLats(self.position, bunch, self.lattice_length)

    def setPosition(self, position: float) -> None:
        self.position = position

    def setLatticeLength(self, length: float) -> None:
        self.lattice_length = length

    def activate(self) -> None:
        self.active = True

    def deactivate(self) -> None:
        self.active = False

    def resetFile(self, file: IO[str]) -> None:
        self.file = file
        self.statlats.resetFile(self.file)


class TeapotMomentsNode(DriftTEAPOT):
    def __init__(self, filename: str, order: int, no_dispersion: bool = True, emit_norm: bool = False, name: str = "moments no name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.moments = Moments(filename, order, no_dispersion, emit_norm)
        self.setType("moments teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattice_length = 0.0
        self.file_out = open(filename, "w")

    def track(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        self.moments.writeMoments(self.position, bunch, self.lattice_length)

    def setPosition(self, position: float) -> None:
        self.position = position

    def closeMoments(self) -> None:
        self.file_out.close()

    def setLatticeLength(self, length: float) -> None:
        self.lattice_length = length


class TeapotMomentsNodeSetMember(DriftTEAPOT):
    def __init__(self, file: IO[str], order: int, no_dispersion: bool = True, emit_norm: bool = False, name: str = "moments no name") -> None:
        DriftTEAPOT.__init__(self, str(name))

        self.file = file
        self.moments = MomentsSetMember(self.file, order, no_dispersion, emit_norm)
        self.setType("moments teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattice_length = 0.0
        self.active = True

    def track(self, params_dict: dict) -> None:
        if self.active:
            bunch = params_dict["bunch"]
            self.moments.writeMoments(self.position, bunch, self.lattice_length)

    def setPosition(self, position: float) -> None:
        self.position = position

    def setLatticeLength(self, length: float) -> None:
        self.lattice_length = length

    def activate(self) -> None:
        self.active = True

    def deactivate(self) -> None:
        self.active = False

    def resetFile(self, file: IO[str]) -> None:
        self.file = file
        self.moments.resetFile(self.file)


class TeapotTuneAnalysisNode(DriftTEAPOT):
    """Estimates tunes from coordinates on neighboring turns.

    This node computes the tunes and actions of each particle in the bunch.
    We use the Average Phase Advance (APA) method to estimate the tunes [1].
    We use only a single turn rather than the average of multiple turns.

    [1] https://cds.cern.ch/record/292773/files/p147.pdf
    [2] https://arxiv.org/pdf/1207.5526
    [3] S. Y. Lee, *Accelerator Physics*
    """
    def __init__(self, name: str = "TeapotTuneAnalysis no name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.tune_calc = BunchTuneAnalysis()
        self.setType("tune calculator teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.active = True
        self.keys = ["phase_1", "phase_2", "tune_1", "tune_2", "action_1", "action_2"]

    def track(self, params_dict: dict) -> None:
        if self.active:
            self.tune_calc.analyzeBunch(params_dict["bunch"])

    def activate(self) -> None:
        self.active = True

    def deactivate(self) -> None:
        self.active = False

    def setPosition(self, position: float) -> None:
        self.position = position

    def setNormMatrix(self, norm_matrix: np.ndarray) -> None:
        ndim = norm_matrix.shape[0]
        for i in range(ndim):
            for j in range(ndim):
                self.tune_calc.setNormMatrixElement(i, j, norm_matrix[i, j])

    def getNormMatrix(self) -> np.ndarray:
        norm_matrix = np.zeros((6, 6))
        for i in range(6):
            for j in range(6):
                norm_matrix[i][j] = self.tune_calc.getNormMatrixElement(i, j)
        return norm_matrix

    def setNormMatrixFromTwiss(
        self, 
        betax: float, 
        alphax: float, 
        etax: float, 
        etapx: float, 
        betay: float, 
        alphay: float,
    ) -> None:
        """Set normalization matrix from Twiss parameters (x, y) and dispersion.
        
        betax{y}: Beta parameter in x{y} plane.
        alphax{y}: Alpha parameter in x{y} plane.
        etax: Dispersion in x plane.
        etapx: Disperion prime in x plane.
        """
        self.tune_calc.setNormMatrixFromTwiss(betax, alphax, etax, etapx, betay, alphay)

    def setNormMatrixFromTransferMatrix(self, transfer_matrix: np.ndarray) -> None:
        norm_matrix = build_norm_matrix_from_tmat(transfer_matrix)
        self.setNormMatrix(norm_matrix)

    def setNormMatrixFromCovMatrix(self, cov_matrix: np.ndarray) -> None:
        # Assume that S = M S M^T, where S is the covariance matrix
        # and M is the transfer matrix. Then M and SU (U is the Poisson matrix)
        # have different eigenvalues but the same eigenvectors So we can compute
        # the normalization matrix directly from SU, without knowing M.
        #
        # However, I'm not sure how to order the eigenvectors of SU. there is
        # no guaranteed ordering from np.linalg.eig. By default, we sort the
        # eigenvectors of SU by their eigenvalues (eigenemittances), so the
        # smallest eigenemittance is mode 1, the next is mode 2, and so on.
        # So if you compare this method to `setNormMatrixFromTransferMatrix`,
        # you may get {nu1, nu2} -> {nu2, nu1}.
        #
        # This is only a problem in coupled lattices with 4D normalization.
        # With 2D normalization there is no ambiguity. This function will
        # check if there are off-block-diagonal terms in the covariance
        # matrix to determine whether to use 2D or 4D normalization.
        norm_matrix = build_norm_matrix_from_cov(cov_matrix)
        self.setNormMatrix(norm_matrix)

    def getData(self, bunch: Bunch, index: int = None) -> dict[str, float] | dict[str, np.ndarray]:
        """Return tune and action data.
        
        Args:
            bunch: A Bunch object.
            index: Particle index. If None, return data for all particles.
        
        Returns:
            data: Dictionary with the following keys:
                - "phase_1"
                - "phase_2"
                - "tune_1"
                - "tune_2"
                - "action_1"
                - "action_2"
                - "action_3"

                If `index` is provided, each value is a float. Otherwise each value
                is a list of floats. If the lattice is uncoupled, 1->x and 2->y.
        """
        data = {}
        if index is None:
            for j, key in enumerate(self.keys):
                data[key] = []
                for index in range(bunch.getSize()):            
                    value = bunch.partAttrValue("ParticlePhaseAttributes", index, j)
                    data[key].append(value)
                data[key] = np.array(data[key])
        else:
            index = int(index)
            bunch_size = bunch.getSize()
            if (index < bunch_size):
                raise ValueError("particle index < 0")
            if (index > bunch_size - 1):
                raise ValueError("particle index > bunch.getSize() - 1")
            for j, key in enumerate(self.keys):
                data[key] = bunch.partAttrValue("ParticlePhaseAttributes", index, j)
        return data
    
    def getTunes(self, bunch: Bunch, index: int = None) -> tuple[float, float] | tuple[np.ndarray, np.ndarray]:
        """Return fractional tunes (nu_1, nu_2).
        
        Args:
            bunch: A Bunch object.
            index: Particle index (not ID).
        
        Returns:
            tune_1: Fractional tune (mode 1).
            tune_2: Fractional tune (mode 2).
        """
        data = self.getData(bunch, index)
        return tuple([data[key] for key in ["tune_1", "tune_2"]])

    def getActions(self, bunch: Bunch, index: int) -> tuple[float, float] | tuple[np.ndarray, np.ndarray]:
        """Return actions (J_1, J_2).
        
        Args:
            bunch: Bunch object.
            index: Particle index (not ID).
        
        Returns:
            J_1: Action (mode 1).
            J_2: Action (mode 2).
        """
        data = self.getData(bunch, index)
        return tuple([data[key] for key in ["action_1", "action_2"]])


class TeapotBPMSignalNode(DriftTEAPOT):
    def __init__(self, name: str = "BPMSignal no name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.bpm = BPMSignal()
        self.setType("BPMSignal")
        self.setLength(0.0)
        self.position = 0.0

    def track(self, params_dict: dict) -> None:
        bunch = params_dict["bunch"]
        self.bpm.analyzeSignal(params_dict["bunch"])

    def setPosition(self, position: float) -> None:
        self.position = position

    def getSignal(self) -> tuple[float, float]:
        x_avg = self.bpm.getSignalX()
        y_avg = self.bpm.getSignalY()
        return x_avg, y_avg
