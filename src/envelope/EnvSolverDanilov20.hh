#ifndef ENV_SOLVER_DANILOV_20_H
#define ENV_SOLVER_DANILOV_20_H

#include "Bunch.hh"
#include "CppPyWrapper.hh"

using namespace std;

/** Envelope solver for the {2, 0} Danilov distribution (upright KV distribution.) 
 The ellipse in the x-y plane can be parameterized as:
     x = cx * cos(psi),
     y = cy * sin(psi),
 where 0 <= psi <= 2pi. The first particle in the bunch is used to track the
 envelope parameters {a, b}, which are used to apply space charge kicks to the
 other particles in the bunch.
 */
class EnvSolverDanilov20 : public OrbitUtils::CppPyWrapper {
public:
  EnvSolverDanilov20(double perveance, double emittanceX, double emittanceY);
  void trackBunch(Bunch *bunch, double length);
  void setPerveance(double perveance);
  void setEmittanceX(double emittanceX);
  void setEmittanceY(double emittanceY);
  double getPerveance();
  double getEmittanceX();
  double getEmittanceY();

private:
  double Q;     // beam perveance
  double epsX; // (4 * sqrt(<xx><x'x'> - <xx'><xx'>))
  double epsY; // (4 * sqrt(<yy><y'y'> - <yy'><yy'>))
};

#endif