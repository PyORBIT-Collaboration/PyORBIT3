import math
import random

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import SpaceChargeCalcUnifEllipse
from orbit.space_charge.sc2p5d import SCUnifEllipse2D_AccNode
from orbit.space_charge.sc2p5d import setSCUnifEllipse2DAccNodes
from orbit.space_charge.sc3d import SCUnifEllipse_AccNode
from orbit.space_charge.sc3d import setSCUnifEllipseAccNodes


def test_unif_ellipse_node_track_bunch():
    bunch = Bunch()
    bunch.mass(0.938)
    bunch.getSyncParticle().kinEnergy(0.001)

    for i in range(10):
        x = random.gauss(0.0, 0.001)
        y = random.gauss(0.0, 0.001)
        z = random.gauss(0.0, 0.001)
        bunch.addParticle(x, 0.0, y, 0.0, z, 0.0)

    params_dict = {"bunch": bunch}

    calculator = SpaceChargeCalcUnifEllipse()
    node = SCUnifEllipse2D_AccNode(calculator)
    node.track(params_dict)


def test_unif_ellipse_node_setters_are_exported():
    assert callable(setSCUnifEllipseAccNodes)
    assert callable(setSCUnifEllipse2DAccNodes)
