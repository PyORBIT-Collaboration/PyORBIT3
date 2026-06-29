## \namespace orbit::bunch_generators
## \brief The classes for bunch generation according to different models.
##
## Classes:


from .TeapotApertureNode import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from .ApertureLatticeModifications import addTeapotApertureNode
from .ApertureLatticeRangeModifications import addCircleApertureSet, addEllipseApertureSet, addRectangleApertureSet

from orbit.core.aperture import Aperture

# from TeapotApertureShapeNode import CircleApertureNode

__all__ = []
__all__.append("Aperture")
# __all__.append("checkBunch")  # non-static method of the Aperture Class
