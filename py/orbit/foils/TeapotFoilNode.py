"""
This module is a foil node class for TEAPOT lattice
"""

import os
import math

# import the auxiliary classes
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccNode, AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

# import Foil class
from orbit.core.foil import Foil


class TeapotFoilNode(DriftTEAPOT):
    """
    The foil node class for TEAPOT lattice
    """

    def __init__(self, xmin, xmax, ymin, ymax, thick, name="foil no name"):
        """
        Constructor. Creates the Foil TEAPOT element.
        """
        DriftTEAPOT.__init__(self, name)
        self.foil = Foil(xmin, xmax, ymin, ymax, thick)
        self.setType("foil teapot")
        self.setLength(0.0)
        # Scattering routine: 0 = full (default), 1 = simple, 2 = geometrical.
        self.scatterChoice = 0

    def track(self, paramsDict):
        """
        The foil-teapot class implementation of the AccNodeBunchTracker class track(probe) method.
        """
        length = self.getLength(self.getActivePartIndex())
        bunch = paramsDict["bunch"]
        if self.scatterChoice == 0:
            self.foil.traverseFoilFullScatter(bunch, paramsDict["lostbunch"])
        elif self.scatterChoice == 1:
            self.foil.traverseFoilSimpleScatter(bunch)
        else:
            self.foil.traverseFoilGeometricalScatter(bunch)


    def setScatterChoice(self, choice):
        self.scatterChoice = choice

    def getnHits(self):
        #returns the particle foil hits
        return self.foil.getnHits()
