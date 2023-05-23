## \namespace orbit::utils::fitting
##
## Classes:
## - PolynomialFit - fitting a Function or SplineCH instances with the plynomial


from orbit.utils.fitting.PolynomialFit import PolynomialFit

from orbit.utils.fitting.SimplexSearch import SimplexSearchAlgorithm
from orbit.utils.fitting.RandomSearch import RandomSearchAlgorithm
from orbit.utils.fitting.GoldenSectionSearch1D import GoldenSectionSearchAlgorithm
from orbit.utils.fitting.BisectionSearch1D import BisectionSearchAlgorithm

from orbit.utils.fitting.Solver_lib import Solver
from orbit.utils.fitting.Solver_lib import TrialPoint
from orbit.utils.fitting.Solver_lib import SolveStopper
from orbit.utils.fitting.Solver_lib import SolveStopperFactory
from orbit.utils.fitting.Solver_lib import ScoreboardActionListener
from orbit.utils.fitting.Solver_lib import VariableProxy
from orbit.utils.fitting.Solver_lib import Scorer
from orbit.utils.fitting.Solver_lib import SearchAgorithm

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
