# -----------------------------------------------------------------------
# -----Test of the Search Algorithms
#  1D Search:
#      BisectionSearchAlgorithm
#      GoldenSectionSearchAlgorithm
# -------------------------------------------
#  General Search (1d,2D, ...):
#      Simplex
#      Random search
# -----------------------------------------------------------------------

import os
import math
import sys
import time
import orbit.core
import pytest

from orbit.utils.fitting.Solver_lib import Solver, Scorer, SolveStopperFactory, VariableProxy, TrialPoint

from orbit.utils.fitting.BisectionSearch1D import BisectionSearchAlgorithm
from orbit.utils.fitting.GoldenSectionSearch1D import GoldenSectionSearchAlgorithm
from orbit.utils.fitting.SimplexSearch import SimplexSearchAlgorithm
from orbit.utils.fitting.RandomSearch import RandomSearchAlgorithm

from orbit.utils.fitting.BisectionSearch1D import BisectionSearchAlgorithm
from orbit.utils.fitting.GoldenSectionSearch1D import GoldenSectionSearchAlgorithm
from orbit.utils.fitting.SimplexSearch import SimplexSearchAlgorithm
from orbit.utils.fitting.RandomSearch import RandomSearchAlgorithm

# ------------------------------------------------------
#  Functions for minimization
# ------------------------------------------------------


class SquareValueScorer1D(Scorer):
    """
    The implementation of the abstract Score class
    for search algorithms testing.
    """

    def __init__(self):
        self.min_pos = 2.0
        self.base = 10.0

    def getScore(self, trialPoint):
        x = trialPoint.getVariableProxyArr()[0].getValue()
        y = (x - self.min_pos) ** 2 + self.base
        return y

    def getAnswer(self):
        """
        Returns true position of the minimum and base value.
        """
        return (self.min_pos, self.base)


class SquareValueScorer3D(Scorer):
    """
    The implementation of the abstract Score class
    for search algorithms testing.
    """

    def __init__(self):
        self.min_pos_arr = [2.0, 3.0, 4.0]
        self.base = 5.0

    def getScore(self, trialPoint):
        x0 = trialPoint.getVariableProxyArr()[0].getValue()
        x1 = trialPoint.getVariableProxyArr()[1].getValue()
        x2 = trialPoint.getVariableProxyArr()[2].getValue()
        x_arr = self.min_pos_arr
        y = (x0 - x_arr[0]) ** 2 + (x1 - x_arr[1]) ** 2 + (x2 - x_arr[2]) ** 2 + self.base
        return y

    def getAnswer(self):
        """
        Returns true position of the minimum and base value.
        """
        return (self.min_pos_arr, self.base)


# -------------------------------------------------------
# Convenience function for testing
# -------------------------------------------------------


def Test_FitTest(scorer, searchAlgorithm, variableProxy_arr, maxIter=10):
    """
    Test for 1D search algorithms.
    """

    # ---- max number of iteration is maxIter
    solverStopper = SolveStopperFactory.maxIterationStopper(maxIter)

    solver = Solver()
    solver.setAlgorithm(searchAlgorithm)
    solver.setStopper(solverStopper)

    trialPoint = TrialPoint()
    for variableProxy in variableProxy_arr:
        trialPoint.addVariableProxy(variableProxy)

    solver.solve(scorer, trialPoint)

    output_string = ""

    # ---- the fitting process ended, now about results
    search_alg_name = searchAlgorithm.getName()
    output_string += "==============Fitting results for algorithm: " + search_alg_name + "\n"
    solver.getScoreboard().printScoreBoard()
    (pos_min_arr, value_min) = scorer.getAnswer()
    output_string += "===== best score ========== exact answer (pos. min, min value)=" + str((pos_min_arr, value_min)) + "\n"
    bestScore = format(solver.getScoreboard().getBestScore(), ".12g")
    output_string += "best score=" + str(bestScore) + " iteration=" + str(solver.getScoreboard().getIteration()) + "\n"
    trialPoint = solver.getScoreboard().getBestTrialPoint()
    output_string += trialPoint.textDesciption() + "\n"
    success_limit = 0.01
    return output_string

    if abs(bestScore - value_min) > success_limit:
        print("Search failed! Stop!")
        sys.exit(1)


# --------------------------------------------
#    Main part of the tests
# --------------------------------------------

scorer = SquareValueScorer1D()
searchAlgorithm = BisectionSearchAlgorithm()

# ---- for Bisection serach we set (name,initial value, step)
# ---- and limits for serach. The step parameter will not be used.
variableProxy = VariableProxy("(x-2)^2+10", 0.0, 0.3333)
variableProxy.setLowerLimit(-5.0)
variableProxy.setUpperLimit(+10.0)

bisection_output = Test_FitTest(
    scorer,
    searchAlgorithm,
    [
        variableProxy,
    ],
)
print(bisection_output)
# -----------------------------------

scorer = SquareValueScorer1D()
searchAlgorithm = GoldenSectionSearchAlgorithm()

# ---- for Golden Section serach we set (name,initial value, step)
# ---- and limits for serach. The step parameter will not be used.
variableProxy = VariableProxy("(x-2)^2+10", 0.0, 0.3333)
variableProxy.setLowerLimit(-5.0)
variableProxy.setUpperLimit(+10.0)

goldensection_output = Test_FitTest(
    scorer,
    searchAlgorithm,
    [
        variableProxy,
    ],
)
print(goldensection_output)

# -----------------------------------

scorer = SquareValueScorer1D()
searchAlgorithm = SimplexSearchAlgorithm()
# searchAlgorithm = RandomSearchAlgorithm()

# ---- for Simplex serach we set (name,initial value, initial search step)
# ---- The *math.sqrt(2.) introduced to avoid exact solutions
variableProxy = VariableProxy("(x-2)^2+10", 0.0, 0.5 * math.sqrt(2.0))

# ---- Number of iterations for SimplexSearchAlgorithm is 20
# ---- but for RandomSearchAlgorithm it should be 50
solvex_output_1D = Test_FitTest(
    scorer,
    searchAlgorithm,
    [
        variableProxy,
    ],
    20,
)
print(solvex_output_1D)
# -----------------------------------

scorer = SquareValueScorer3D()
searchAlgorithm = SimplexSearchAlgorithm()
# searchAlgorithm = RandomSearchAlgorithm()

# ---- for Simplex serach we set (name,initial value, initial search step)
# ---- The *math.sqrt(2.) introduced to avoid exact solutions
variableProxy0 = VariableProxy("one", 1.0, 0.2 * math.sqrt(2.0))
variableProxy1 = VariableProxy("two", 2.0, 0.2 * math.sqrt(2.0))
# ---- test that we can switch off a varibale from the minimization process
# variableProxy1.setUseInSolver(False)
variableProxy2 = VariableProxy("three", 3.0, 0.2 * math.sqrt(2.0))
variableProxy_arr = []
variableProxy_arr.append(variableProxy0)
variableProxy_arr.append(variableProxy1)
variableProxy_arr.append(variableProxy2)

# ---- Number of iterations for SimplexSearchAlgorithm is 50
# ---- but for RandomSearchAlgorithm it should be 200
solvex_output_3D = Test_FitTest(scorer, searchAlgorithm, variableProxy_arr, 50)
print(solvex_output_3D)


def test_bisection_search():
    expected_bisection_output = """==============Fitting results for algorithm: Bisection Search
===== best score ========== exact answer (pos. min, min value)=(2.0, 10.0)
best score=10.0007476807 iteration=10
======== TrialPoint ===========
 Name                       Value          Step       Use      Limit_Min       Limit_Max  
               (x-2)^2+10         1.972656          0.3333    1                -5              10  
"""
    assert bisection_output == expected_bisection_output


def test_goldensection_solve():
    expected_goldensection_output = """==============Fitting results for algorithm: Golden Section Search
===== best score ========== exact answer (pos. min, min value)=(2.0, 10.0)
best score=10.0067304508 iteration=10
======== TrialPoint ===========
 Name                       Value          Step       Use      Limit_Min       Limit_Max  
               (x-2)^2+10         2.082039          0.3333    1                -5              10  
"""
    assert goldensection_output == expected_goldensection_output


def test_simplex_solve_1D():
    expected_solvex_output_1D = """==============Fitting results for algorithm: Simplex Search
===== best score ========== exact answer (pos. min, min value)=(2.0, 10.0)
best score=10.0001268367 iteration=20
======== TrialPoint ===========
 Name                       Value          Step       Use      Limit_Min       Limit_Max  
               (x-2)^2+10         1.988738      0.08838835    1    -1.797693e+300   1.797693e+300  
"""
    assert solvex_output_1D == expected_solvex_output_1D


def test_simplex_solve_3D():
    expected_solvex_output_3D = """==============Fitting results for algorithm: Simplex Search
===== best score ========== exact answer (pos. min, min value)=([2.0, 3.0, 4.0], 5.0)
best score=5.00004003151 iteration=50
======== TrialPoint ===========
 Name                       Value          Step       Use      Limit_Min       Limit_Max  
                      one         1.999807     0.009020704    1    -1.797693e+300   1.797693e+300  
                      two         3.001916      0.01349468    1    -1.797693e+300   1.797693e+300  
                    three         4.006027      0.01065753    1    -1.797693e+300   1.797693e+300  
"""
    assert solvex_output_3D == expected_solvex_output_3D
