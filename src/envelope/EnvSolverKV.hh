#ifndef ENV_SOLVER_KV_H
#define ENV_SOLVER_KV_H

#include "Bunch.hh"
#include "CppPyWrapper.hh"

using namespace std;

/** Envelope solver for upright KV distribution in uncoupled lattice.

The first particle in the bunch is used to track the
envelope parameters. The envelope parameters are use to apply space charge
kicks to the other particles in the bunch.

This class does not yet handle boudary conditions or dispersion.
 */
class EnvSolverKV : public OrbitUtils::CppPyWrapper {
public:
  EnvSolverKV(double perveance, double emittanceX, double emittanceY);
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