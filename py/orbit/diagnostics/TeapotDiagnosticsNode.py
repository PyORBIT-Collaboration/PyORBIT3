"""TEAPOT-style bunch diagnostic nodes."""

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

from orbit.core.bunch import BunchTuneAnalysis
from orbit.core.bunch import BunchTuneAnalysis4D


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
    """Tune analysis node.
    
    The tunes are estimated from the particle phase space coordinates on
    neighboring turns. The coordinates are "normalized" in each 2D phase 
    plane (x-x', y-y') using provided twiss parameters (alpha, beta) as 
    well as the dispersion in the x plane. In the normalized frame, the 
    turn-by-turn coordinates advance in phase around a circle. The fractional
    tune is defined by the change in angle for a single turn.

    The calculation will be incorrect if the normalization matrix does not
    map the turn-by-turn coordinates to a circle. In a linear lattice
    without space charge, the normalization will be correct if the provided
    twiss parameters are the same as the periodic lattice twiss parameters
    at the location of the node.
    """
    def __init__(self, name: str = "tuneanalysis no name") -> None:
        """Constructor."""
        DriftTEAPOT.__init__(self, name)
        self.bunchtune = BunchTuneAnalysis()
        self.setType("tune calculator teapot")
        self.lattlength = 0.0
        self.setLength(0.0)
        self.position = 0.0

    def track(self, paramsDict: dict) -> None:
        """Implementation of the AccNodeBunchTracker class track(probe) method."""
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.bunchtune.analyzeBunch(bunch)

    def setPosition(self, pos):
        self.position = pos

    def setLatticeLength(self, lattlength):
        self.lattlength = lattlength

    def assignTwiss(
        self, 
        betax: float, 
        alphax: float, 
        etax: float, 
        etapx: float, 
        betay: float, 
        alphay: float,
    ) -> None:
        """Set the twiss parameters for the coordinate normalization.
        
        betax{y}: Beta parameter in x{y} plane.
        alphax{y}: Alpha parameter in x{y} plane.
        etax: Dispersion in x plane.
        etapx: Disperion prime in x plane.
        """
        self.bunchtune.assignTwiss(betax, alphax, etax, etapx, betay, alphay)


class TeapotTuneAnalysis4DNode(DriftTEAPOT):
    """Tune analysis node with 4D normalization.
    
    The tunes are estimated from the particle phase space coordinates on
    neighboring turns. Each set of 4D phase space coordinates is normalized
    with the provided 4 x 4 matrix.

    In the normalized frame (if the matrix is properly chosen), the 
    turn-by-turn coordinates advance in phase around a circle. The fractional
    tune is defined by the change in angle on a single turn, divided by
    2 pi.

    Note that any normalization matrix is allowed. It is up to the user 
    to ensure proper normalization, i.e., based on the one-turn transfer
    matrix.
    """
    def __init__(self, name: str = "tuneanalysis4d no name") -> None:
        DriftTEAPOT.__init__(self, name)
        self.bunchtune = BunchTuneAnalysis4D()
        self.setType("tune calculator 4d teapot")
        self.lattlength = 0.0
        self.setLength(0.0)
        self.position = 0.0

    def track(self, paramsDict: dict) -> None:
        bunch = paramsDict["bunch"]
        self.bunchtune.analyzeBunch(bunch)

    def setPosition(self, pos: float) -> None:
        self.position = pos

    def setLatticeLength(self, lattlength: float) -> None:
        self.lattlength = lattlength

    def setNormMatrix(self, norm_matrix: list[list[float]]) -> None:
        norm_matrix_list = list(norm_matrix)
        for i in range(4):
            for j in range(4):
                value = float(norm_matrix_list[i][j])
                self.bunchtune.setNormMatrixElement(i, j, value)

    def getNormMatrix(self) -> list[list[float]]:
        norm_matrix = [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
        for i in range(4):
            for j in range(4):
                norm_matrix[i][j] = self.bunchtune.getNormMatrixElement(i, j)
        return norm_matrix


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
