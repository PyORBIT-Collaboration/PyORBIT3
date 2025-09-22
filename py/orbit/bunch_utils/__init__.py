## \namespace orbit::bunch_utils
## \brief These classes are for bunch utilities
##
## Classes:
## - ParticleIdNumber  - Class for adding unique id numbers to particle in a bunch
#

from orbit.bunch_utils.particleidnumber import ParticleIdNumber

# This guards against missing numpy.
# Should be imporved with some meaningful (and MPI friendly?) warning printed out.
try:
    from orbit.bunch_utils.serialize import collect_bunch, save_bunch, load_bunch
except:
    pass

__all__ = []
__all__.append("addParticleIdNumbers")
