## \namespace orbit::sc1d
## \brief The classes and functions for foils
##
## Classes:
##   sc1DNode -longitudinal space charge node for the TEAPOT lattices
##
## Functions:
##   addLongitudinalSpaceChargeNode- function to add one longitudinal
##   space charge node to the lattice

from pyorbit.space_charge.sc1d.sc1DNode import SC1D_AccNode
from pyorbit.space_charge.sc1d.sc1DNode import FreqDep_SC1D_AccNode
from pyorbit.space_charge.sc1d.sc1DNode import BetFreqDep_SC1D_AccNode

from pyorbit.space_charge.sc1d.scLatticeModifications import addLongitudinalSpaceChargeNode
from pyorbit.space_charge.sc1d.scLatticeModifications import addLongitudinalSpaceChargeNodeAsChild

__all__ = []
__all__.append("sc1DNode")
__all__.append("addLongitudinalSpaceChargeNode")
__all__.append("addLongitudinalSpaceChargeNodeAsChild")
