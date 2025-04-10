from ..lattice import AccActionsContainer
from ..lattice import AccLattice
from ..lattice import AccNode
from ..lattice import AccNodeBunchTracker
from ..utils import orbitFinalize

from orbit.ext.danilov_envelope import DanilovEnvelopeSolver20
from orbit.ext.danilov_envelope import DanilovEnvelopeSolver22


class DanilovEnvelopeSolverNode(AccNodeBunchTracker):
    def __init__(
        self,
        solver: DanilovEnvelopeSolver20 | DanilovEnvelopeSolver22,
        name: str = None,
        kick_length: float = 0.0,
        perveance: float = 0.0,
    ) -> None:
        super().__init__(name=name)
        self.setType("DanilovEnvSolver")
        self.setLength(0.0)

        self.solver = solver
        self.kick_length = kick_length
        self.perveance = perveance
        self.active = True

    def set_active(self, setting: bool) -> None:
        self.active = setting

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

    def set_perveance(self, perveance: float) -> None:
        self.solver.setPerveance(perveance)

    def set_kick_length(self, kick_length: float) -> None:
        self.kick_length = kick_length


class DanilovEnvelopeSolverNode20(DanilovEnvelopeSolverNode):
    def __init__(self, eps_x: float, eps_y: float, **kwargs) -> None:
        super().__init__(solver=None, **kwargs)
        self.eps_x = eps_x
        self.eps_y = eps_y
        self.solver = DanilovEnvelopeSolver20(self.perveance, self.eps_x, self.eps_y)

    def set_emittances(self, eps_x: float, eps_y: float) -> None:
        self.eps_x = eps_x
        self.eps_y = eps_y
        self.solver.setEmittanceX(eps_x)
        self.solver.setEmittanceY(eps_y)


class DanilovEnvelopeSolverNode22(DanilovEnvelopeSolverNode):
    def __init__(self, **kwargs) -> None:
        super().__init__(solver=None, **kwargs)
        self.solver = DanilovEnvelopeSolver22(self.perveance)
