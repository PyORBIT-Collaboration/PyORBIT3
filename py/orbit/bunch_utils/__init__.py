## \namespace orbit::bunch_utils
## \brief These classes are for bunch utilities
##
## Classes:
## - ParticleIdNumber  - Class for adding unique id numbers to particle in a bunch
#

from .particleidnumber import ParticleIdNumber

# This guards against missing numpy.
# Should be imporved with some meaningful (and MPI friendly?) warning printed out.
try:
    from .serialize import collect_bunch, save_bunch, load_bunch
    from .serialize import BunchDict
    from .serialize import FileHandler, NumPyHandler
except:
    pass

__all__ = []
# __all__.append("addParticleIdNumbers")  # doesn't exist
__all__.append("ParticleIdNumber")
__all__.append("collect_bunch")
__all__.append("save_bunch")
__all__.append("load_bunch")
__all__.append("BunchDict")
__all__.append("FileHandler")
__all__.append("NumPyHandler")

