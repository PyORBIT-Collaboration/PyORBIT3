from typing import Any

from ..lattice import AccActionsContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..lattice import AccNodeBunchTracker
from ..utils import orbitFinalize

from orbit.core.envelope import EnvSolverDanilov20
from orbit.core.envelope import EnvSolverDanilov22


class EnvSolverNode(AccNodeBunchTracker):
    def __init__(
        self,
        solver: Any,
        name: str = "",
        kick_length: float = 0.0,
    ) -> None:
        super().__init__(name=name)
        self.setType("EnvSolver")
        self.setLength(0.0)

        self.solver = None
        self.kick_length = kick_length
        self.active = True

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
        self.solver.trackBunch(bunch, self.kick_length)


class Danilov22EnvSolverNode(EnvSolverNode):
    def __init__(self, perveance: float, **kwargs) -> None:
        super().__init__(**kwargs)
        self.solver = EnvSolverDanilov22(perveance)

    def setPerveance(self, perveance: float) -> None:
        self.solver.setPerveance(perveance)


class Danilov20EnvSolverNode(EnvSolverNode):
    def __init__(self, eps_x: float, eps_y: float, perveance: float, **kwargs) -> None:
        super().__init__(**kwargs)
        self.solver = EnvSolverDanilov20(perveance, eps_x, eps_y)

    def setPerveance(self, perveance: float) -> None:
        self.solver.setPerveance(perveance)

    def setEmittances(self, eps_x: float, eps_y: float) -> None:
        self.solver.setEmittanceX(eps_x)
        self.solver.setEmittanceY(eps_y)