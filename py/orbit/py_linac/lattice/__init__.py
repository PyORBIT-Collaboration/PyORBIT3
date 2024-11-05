## \namespace orbit::py_linac::lattice
## \brief The base classes of ORBIT Linac lattice structure.
##
## Classes:
## - LinacAcclattice       - Class. The linac lattice.
## - LinacAccNodes         - Module. Collection of the linac accelerator nodes: drifts, quads, RF gaps etc..
## - LinacRfGapNodes       - Module. Collection of RF Gap models

from orbit.py_linac.lattice.LinacAccLatticeLib import LinacAccLattice, RF_Cavity, Sequence

from orbit.py_linac.lattice.LinacAccNodes import BaseLinacNode, LinacNode, LinacMagnetNode
from orbit.py_linac.lattice.LinacAccNodes import MarkerLinacNode, Drift, Quad, AbstractRF_Gap, Bend
from orbit.py_linac.lattice.LinacAccNodes import Solenoid
from orbit.py_linac.lattice.LinacAccNodes import DCorrectorH, DCorrectorV
from orbit.py_linac.lattice.LinacAccNodes import ThickKick

from orbit.py_linac.lattice.LinacRfGapNodes import BaseRF_Gap, AxisFieldRF_Gap, RF_AxisFieldsStore

from orbit.py_linac.lattice.LinacApertureNodes import LinacApertureNode
from orbit.py_linac.lattice.LinacApertureNodes import CircleLinacApertureNode
from orbit.py_linac.lattice.LinacApertureNodes import EllipseLinacApertureNode
from orbit.py_linac.lattice.LinacApertureNodes import RectangleLinacApertureNode
from orbit.py_linac.lattice.LinacApertureNodes import LinacPhaseApertureNode
from orbit.py_linac.lattice.LinacApertureNodes import LinacEnergyApertureNode

from orbit.py_linac.lattice.LinacFieldOverlappingNodes import AxisField_and_Quad_RF_Gap
from orbit.py_linac.lattice.LinacFieldOverlappingNodes import OverlappingQuadsNode

from orbit.py_linac.lattice.LinacAccLatticeFunc import GetGlobalQuadGradient
from orbit.py_linac.lattice.LinacAccLatticeFunc import GetGlobalQuadGradientDerivative
from orbit.py_linac.lattice.LinacAccLatticeFunc import GetGlobalRF_AxisField
from orbit.py_linac.lattice.LinacAccLatticeFunc import getNodeForNameFromWholeLattice
from orbit.py_linac.lattice.LinacAccLatticeFunc import getNodePosDictForWholeLattice
from orbit.py_linac.lattice.LinacAccLatticeFunc import getAllNodesInLattice
from orbit.py_linac.lattice.LinacAccLatticeFunc import getAllMagnetsInLattice

from orbit.py_linac.lattice.LinacTransportMatrixGenNodes import LinacTrMatrixGenNode
from orbit.py_linac.lattice.LinacTransportMatrixGenNodes import LinacTrMatricesController

from orbit.py_linac.lattice.LinacDiagnosticsNodes import LinacBPM

__all__ = []
__all__.append("LinacAccLattice")

# AccNodes
__all__.append("BaseLinacNode")
__all__.append("LinacNode")
__all__.append("LinacMagnetNode")
__all__.append("MarkerLinacNode")
__all__.append("Drift")
__all__.append("Quad")
__all__.append("AbstractRF_Gap")
__all__.append("DCorrectorH")
__all__.append("DCorrectorV")
__all__.append("ThickKick")
__all__.append("Bend")
__all__.append("Solenoid")


__all__.append("LinacApertureNode")
__all__.append("CircleLinacApertureNode")
__all__.append("EllipseLinacApertureNode")
__all__.append("RectangleLinacApertureNode")
__all__.append("LinacPhaseApertureNode")
__all__.append("LinacEnergyApertureNode")


__all__.append("RF_Cavity")
__all__.append("Sequence")

__all__.append("LinacStructureTree")
__all__.append("LinacStructureSeq")
__all__.append("LinacStuctureNode")

__all__.append("BaseRF_Gap")
__all__.append("AxisFieldRF_Gap")
__all__.append("RF_AxisFieldsStore")

__all__.append("AxisField_and_Quad_RF_Gap")
__all__.append("OverlappingQuadsNode")

__all__.append("GetGlobalQuadGradient")
__all__.append("GetGlobalQuadGradientDerivative")
__all__.append("GetGlobalRF_AxisField")
__all__.append("getNodeForNameFromWholeLattice")
__all__.append("getNodePosDictForWholeLattice")
__all__.append("getAllNodesInLattice")
__all__.append("getAllMagnetsInLattice")

__all__.append("LinacTrMatrixGenNode")
__all__.append("LinacTrMatricesController")

__all__.append("LinacBPM")
