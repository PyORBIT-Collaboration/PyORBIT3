## \namespace orbit::errors
## \brief The classes and functions for errors
##
## Classes:
##   ErrorNode - error node for TEAPOT lattices
##
## Functions:
##   addErrorNode - function to add one error
##                  node to the lattice

__all__ = []
try:
    from .orbit_correction import orbit
    from .orbit_correction import correction
    __all__.append("orbit")
    __all__.append("correction")
except ModuleNotFoundError:
    pass

