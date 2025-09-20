from typing import Any

from ..lattice import AccActionsContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..lattice import AccNodeBunchTracker
from ..utils import orbitFinalize

from orbit.core.envelope import KVEnvelopeTracker
from orbit.core.envelope import DanilovEnvelopeTracker


class EnvelopeTrackerNode(AccNodeBunchTracker):
    def __init__(self, name: str = None, kick_length: float = 0.0) -> None:
        super().__init__(name=name)
        self.setType("EnvelopeTracker")
        self.setLength(0.0)

        self.kick_length = kick_length
        self.active = True
        self.tracker = None

    def setActive(self, setting: bool) -> None:
        self.active = setting

    def setKickLength(self, kick_length: float) -> None:
        self.kick_length = kick_length

    def isRFGap(self) -> bool:
        # In case this node is used in linac tracking
        return False

    def trackDesign(self, params_dict):
        # In case this node is used in linac tracking
        pass

    def track(self, params_dict: dict) -> None:
        if not self.active:
            return
        bunch = params_dict["bunch"]
        self.tracker.trackBunch(bunch, self.kick_length)


class KVEnvelopeTrackerNode(EnvelopeTrackerNode):
    def __init__(self, eps_x: float, eps_y: float, perveance: float, **kwargs) -> None:
        super().__init__(**kwargs)
        self.tracker = KVEnvelopeTracker(perveance, 4.0 * eps_x, 4.0 * eps_y)

    def setPerveance(self, perveance: float) -> None:
        self.tracker.setPerveance(perveance)

    def setEmittances(self, eps_x: float, eps_y: float) -> None:
        self.tracker.setEmittanceX(4.0 * eps_x)
        self.tracker.setEmittanceY(4.0 * eps_y)
        

class DanilovEnvelopeTrackerNode(EnvelopeTrackerNode):
    def __init__(self, perveance: float, **kwargs) -> None:
        super().__init__(**kwargs)
        self.tracker = DanilovEnvelopeTracker(perveance)

    def setPerveance(self, perveance: float) -> None:
        self.tracker.setPerveance(perveance)

        