## \namespace orbit::utils::fitting
##
## Classes:
## - PolynomialFit - fitting a Function or SplineCH instances with the plynomial


from .PolynomialFit import PolynomialFit

from .SimplexSearch import SimplexSearchAlgorithm
from .RandomSearch import RandomSearchAlgorithm
from .GoldenSectionSearch1D import GoldenSectionSearchAlgorithm
from .BisectionSearch1D import BisectionSearchAlgorithm

from .Solver_lib import Solver
from .Solver_lib import TrialPoint
from .Solver_lib import SolveStopper
from .Solver_lib import SolveStopperFactory
from .Solver_lib import ScoreboardActionListener
from .Solver_lib import VariableProxy
from .Solver_lib import Scorer
from .Solver_lib import SearchAgorithm

__all__ = []
__all__.append("PolynomialFit")
__all__.append("SimplexSearchAlgorithm")
__all__.append("RandomSearchAlgorithm")
__all__.append("GoldenSectionSearchAlgorithm")
__all__.append("BisectionSearchAlgorithm")
__all__.append("Solver")
__all__.append("TrialPoint")
__all__.append("SolveStopper")
__all__.append("SolveStopperFactory")
__all__.append("ScoreboardActionListener")
__all__.append("VariableProxy")
__all__.append("Scorer")
__all__.append("SearchAgorithm")
