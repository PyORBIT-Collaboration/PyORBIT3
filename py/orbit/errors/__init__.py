## \namespace orbit::errors
## \brief The classes and functions for errors
##
## Classes:
##   ErrorNode - error node for TEAPOT lattices
##
## Functions:
##   addErrorNode - function to add one error
##                  node to the lattice

from pyorbit.errors.ErrorNode import coorddisplacement
from pyorbit.errors.ErrorNode import longdisplacement
from pyorbit.errors.ErrorNode import straightrotationxy
from pyorbit.errors.ErrorNode import straightrotationxsi
from pyorbit.errors.ErrorNode import straightrotationxsf
from pyorbit.errors.ErrorNode import straightrotationysi
from pyorbit.errors.ErrorNode import straightrotationysf
from pyorbit.errors.ErrorNode import bendfieldi
from pyorbit.errors.ErrorNode import bendfieldf
from pyorbit.errors.ErrorNode import benddisplacementxi
from pyorbit.errors.ErrorNode import benddisplacementxf
from pyorbit.errors.ErrorNode import benddisplacementyi
from pyorbit.errors.ErrorNode import benddisplacementyf
from pyorbit.errors.ErrorNode import benddisplacementli
from pyorbit.errors.ErrorNode import benddisplacementlf
from pyorbit.errors.ErrorNode import rotationi
from pyorbit.errors.ErrorNode import rotationf
from pyorbit.errors.ErrorNode import dipolekicker
from pyorbit.errors.ErrorNode import dipolekickerosc
from pyorbit.errors.ErrorNode import quadkicker
from pyorbit.errors.ErrorNode import quadkickerosc
from pyorbit.errors.ErrorNode import AddErrorNode
from pyorbit.errors.ErrorNode import AddErrorSet

from pyorbit.errors.ErrorLatticeModifications import addErrorNode
from pyorbit.errors.ErrorLatticeModifications import addErrorNodeAsChild
from pyorbit.errors.ErrorLatticeModifications import addErrorNodeAsChild_I
from pyorbit.errors.ErrorLatticeModifications import addErrorNodeAsChild_F

__all__ = []
__all__.append("")
__all__.append("addErrorNode")
__all__.append("addErrorNodeAsChild")
__all__.append("addErrorNodeAsChild_I")
__all__.append("addErrorNodeAsChild_F")
