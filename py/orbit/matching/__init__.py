## \namespace orbit::injection
## \brief These classes are for displacement bumps
##
## Classes:
## - simpleBump  - Class for simple coordinate transverse bump.

## - addTeapotBumpNode - Adds a teapot bump node to a teapot lattice
## - TeapotSimpleBumpNode - Creates a teapot instance of a simple bump nodes

__all__ = []

# guard against missing numpy/scipy
try:
    from .matching import Optics, EnvelopeSolver
    __all__.append("Optics")
    __all__.append("EnvelopeSolver")
except ModuleNotFoundError:
    pass

