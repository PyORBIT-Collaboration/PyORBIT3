"""
This module defines x and y kicker classes for TEAPOT lattice
"""

import os
import math

# import the auxiliary classes
from ..utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from ..lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from ..teapot import DriftTEAPOT

# import injection class
from . import waveforms
from . import XKicker, YKicker


class TeapotXKickerNode(DriftTEAPOT):
    """
    The kicker node class for TEAPOT lattice
    """

    def __init__(self, bunch, kx, waveform, name="kicker"):
        """
        Constructor. Creates the Kicker TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.kicker = XKicker(bunch, kx, waveform)
        self.setType("XKicker")
        self.setLength(0.0)

    def track(self, paramsDict):
        """
        The kicker-teapot class implementation of the
        AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.kicker.kick()


class TeapotYKickerNode(DriftTEAPOT):
    """
    The kicker node class for TEAPOT lattice
    """

    def __init__(self, bunch, ky, waveform, name="kicker"):
        """
        Constructor. Creates the Kicker TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.kicker = YKicker(bunch, ky, waveform)
        self.setType("YKicker")
        self.setLength(0.0)

    def track(self, paramsDict):
        """
        The kicker-teapot class implementation of the
        AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        self.kicker.kick()
