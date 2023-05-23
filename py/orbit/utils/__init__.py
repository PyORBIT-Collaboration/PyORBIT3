## \namespace orbit::utils
## \brief Utility classes.
##
## Classes:
## - multiDimDoubleArray - Method. Generates multi-dimensional array with doubles.
## - multiDimIntArray    - Method. Generates multi-dimensional array with integers.
## - orbitFinalize    - Method. Finalizes ORBIT script execution.
## - NamedObject      - Class. Represents an object with a name.
## - TypedObject      - Class. Represents an object with a type.
## - ParamsDictObject - Class. Represents an object that has a parameters dictionary.

from .multiDimArray import multiDimDoubleArray
from .multiDimArray import multiDimIntArray
from .orbitFinalize import orbitFinalize
from .NamedObject import NamedObject
from .TypedObject import TypedObject
from .ParamsDictObject import ParamsDictObject

from .phaseOperations import phaseNearTargetPhase, phaseNearTargetPhaseDeg
from .consts import speed_of_light

__all__ = []
__all__.append("multiDimDoubleArray")
__all__.append("multiDimIntArray")
__all__.append("orbitFinalize")
__all__.append("NamedObject")
__all__.append("TypedObject")
__all__.append("ParamsDictObject")
__all__.append("phaseNearTargetPhase")
__all__.append("phaseNearTargetPhaseDeg")
__all__.append("speed_of_light")
