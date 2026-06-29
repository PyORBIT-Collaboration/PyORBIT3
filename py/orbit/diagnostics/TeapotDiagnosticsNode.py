"""TEAPOT-style bunch diagnostic nodes."""

import numpy as np

from ..utils import orbitFinalize
from ..utils import NamedObject
from ..utils import ParamsDictObject

from ..lattice import AccNode
from ..lattice import AccActionsContainer
from ..lattice import AccNodeBunchTracker

from ..teapot import DriftTEAPOT

from .diagnostics import StatLats
from .diagnostics import StatLatsSetMember
from .diagnostics import Moments
from .diagnostics import MomentsSetMember
from .diagnostics import BPMSignal

from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTuneAnalysis


class TeapotStatLatsNode(DriftTEAPOT):
    """
    The statlats node class for TEAPOT lattice
    """

    def __init__(self, filename, name="statlats no name"):
        """
        Constructor. Creates the StatLats TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.statlats = StatLats(filename)
        self.setType("statlats teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattlength = 0.0
        self.file_out = open(filename, "w")

    def track(self, paramsDict):
        """
        The statlats-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.statlats.writeStatLats(self.position, bunch, self.lattlength)

    def setPosition(self, pos):
        self.position = pos

    def closeStatLats(self):
        self.file_out.close()

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength


class TeapotStatLatsNodeSetMember(DriftTEAPOT):
    """
    The statlats node class for TEAPOT lattice
    """

    def __init__(self, file, name="statlats no name"):
        """
        Constructor. Creates the StatLats TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.statlats = StatLatsSetMember(file)
        self.setType("statlats teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattlength = 0.0
        self.active = True
        self.file = file

    def track(self, paramsDict):
        """
        The statlats-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        if self.active:
            length = self.getLength(self.getActivePartIndex())
            bunch = paramsDict["bunch"]
            self.statlats.writeStatLats(self.position, bunch, self.lattlength)

    def setPosition(self, pos):
        self.position = pos

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength

    def activate(self):
        self.active = True

    def deactivate(self):
        self.active = False

    def resetFile(self, file):
        self.file = file
        self.statlats.resetFile(self.file)


class TeapotMomentsNode(DriftTEAPOT):
    """
    The moments node class for TEAPOT lattice
    """

    def __init__(self, filename, order, nodispersion=True, emitnorm=False, name="moments no name"):
        """
        Constructor. Creates the StatLats TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.moments = Moments(filename, order, nodispersion, emitnorm)
        self.setType("moments teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattlength = 0.0
        self.file_out = open(filename, "w")

    def track(self, paramsDict):
        """
        The moments-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.moments.writeMoments(self.position, bunch, self.lattlength)

    def setPosition(self, pos):
        self.position = pos

    def closeMoments(self):
        self.file_out.close()

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength


class TeapotMomentsNodeSetMember(DriftTEAPOT):
    """
    The moments node class for TEAPOT lattice
    """

    def __init__(self, file, order, nodispersion=True, emitnorm=False, name="moments no name"):
        """
        Constructor. Creates the Moments TEAPOT element.
        """
        DriftTEAPOT.__init__(self, str(name))

        self.file = file
        self.moments = MomentsSetMember(self.file, order, nodispersion, emitnorm)
        self.setType("moments teapot")
        self.setLength(0.0)
        self.position = 0.0
        self.lattlength = 0.0
        self.active = True

    def track(self, paramsDict):
        """
        The moments-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        if self.active:
            length = self.getLength(self.getActivePartIndex())
            bunch = paramsDict["bunch"]
            self.moments.writeMoments(self.position, bunch, self.lattlength)

    def setPosition(self, pos):
        self.position = pos

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength

    def activate(self):
        self.active = True

    def deactivate(self):
        self.active = False

    def resetFile(self, file):
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
    def __init__(self, name: str = "tuneanalysis no name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.bunchtune = BunchTuneAnalysis()
        self.setType("tune calculator teapot")
        self.lattlength = 0.0
        self.setLength(0.0)
        self.position = 0.0
        self.active = True

        self.keys = ["phase_1", "phase_2", "tune_1", "tune_2", "action_1", "action_2"]

    def track(self, paramsDict: dict) -> None:
        """Implementation of the AccNodeBunchTracker class track(probe) method."""
        if not self.active:
            return
        self.bunchtune.analyzeBunch(paramsDict["bunch"])

    def setActive(self, active: bool) -> None:
        self.active = active

    def setPosition(self, position: float) -> None:
        self.position = position

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength

    def setNormMatrix(self, norm_matrix: np.ndarray) -> None:
        """Set the normalization matrix.
        
        Args:
            norm_matrix: Normalization matrix of shape (4, 4) or (6, 6).
        """
        ndim = norm_matrix.shape[0]
        for i in range(ndim):
            for j in range(ndim):
                self.bunchtune.setNormMatrixElement(i, j, norm_matrix[i, j])

    def getNormMatrix(self) -> np.ndarray:
        """Return normalization matrix of shape (6, 6)."""
        norm_matrix = np.zeros((6, 6))
        for i in range(6):
            for j in range(6):
                norm_matrix[i][j] = self.bunchtune.getNormMatrixElement(i, j)
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
        self.bunchtune.assignTwiss(betax, alphax, etax, etapx, betay, alphay)

    def getData(self, bunch: Bunch, index: int = None) -> dict[str, float] | dict[str, list[float]]:
        """Return tune and action data.
        
        Args:
            bunch: A Bunch object.
            index: Particle index. If None, return data for all particles.
        
        Returns:
            data: Dictionary with keys {"phase_1", "phase_2", "tune_1", "tune_2", 
                "action_1", "action_2" "action_3"}. (If the lattice is uncoupled, 
                1->x and 2->y.) If `index` is None, each value of `data` is a list 
                of floats; otherwise each value is a float.
        """
        data = {}
        if index is None:
            for j, key in enumerate(self.keys):
                data[key] = []
                for index in range(bunch.getSize()):            
                    value = bunch.partAttrValue("ParticlePhaseAttributes", index, j)
                    data[key].append(value)
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
    
    def getTunes(self, bunch: Bunch, index: int = None) -> dict[str, float] | dict[str, list[float]]:
        """Return fractional tunes (nu_1, nu_2).
        
        Args:
            bunch: A Bunch object.
            index: Particle index in bunch (not ID).
        
        Returns:
            tune_1: Fractional tune (mode 1).
            tune_2: Fractional tune (mode 2).
        """
        data = self.getData(bunch, index)
        return tuple([data[key] for key in ["tune_1", "tune_2"]])
    
    
    def getActions(self, bunch: Bunch, index: int) -> dict[str, float] | dict[str, list[float]]:
        """Return actions (J_1, J_2).
        
        Args:
            bunch: A Bunch object.
            index: Particle index in bunch (not ID).
        
        Returns:
            J_1: Action (mode 1).
            J_2: Action (mode 2).
        """
        data = self.getData(bunch, index)
        return tuple([data[key] for key in ["action_1", "action_2"]])


class TeapotBPMSignalNode(DriftTEAPOT):
    def __init__(self, name="BPMSignal no name"):
        """
        Constructor. Creates the StatLats TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.bpm = BPMSignal()
        self.setType("BPMSignal")
        self.lattlength = 0.0
        self.setLength(0.0)
        self.position = 0.0

    def track(self, paramsDict):
        """
        The bunchtuneanalysis-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.bpm.analyzeSignal(bunch)

    def setPosition(self, pos):
        self.position = pos

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength

    def getSignal(self):
        xAvg = self.bpm.getSignalX()
        yAvg = self.bpm.getSignalY()
        return xAvg, yAvg
