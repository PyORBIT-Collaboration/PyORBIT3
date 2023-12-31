## \namespace orbit::collimation
## \brief The classes and functions for collimations
##
## Classes:
##  - TeapotCollimatorNode - collimation node for the TEAPOT lattices
##
## addTeapotCollimatorNode - function to add one collimator node to the lattice

from .TeapotCollimatorNode import TeapotCollimatorNode
from .collimationLatticeModifications import addTeapotCollimatorNode

__all__ = []
__all__.append("TeapotCollimatorNode")
__all__.append("addTeapotCollimatorNode")
