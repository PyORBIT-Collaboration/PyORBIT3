#ifndef DANILOV_ENVELOPE_SOLVER_20_H
#define DANILOV_ENVELOPE_SOLVER_20_H

#include "Bunch.hh"
#include "CppPyWrapper.hh"

using namespace std;

/** Envelope solver for the {2, 0} Danilov distribution (KV distribution).

 The ellipse in the x-y plane can be parameterized as:
     x = cx * cos(psi),
     y = cy * sin(psi),
 where 0 <= psi <= 2pi. The first particle in the bunch is used to track the
 envelope parameters {a, b}, which are used to apply space charge kicks to the
 other particles in the bunch.
 */
class DanilovEnvelopeSolver20 : public OrbitUtils::CppPyWrapper {
public:
  DanilovEnvelopeSolver20(double perveance, double eps_x, double eps_y);
  void trackBunch(Bunch *bunch, double length);
  void setPerveance(double perveance);
  void setEmittanceX(double eps_x);
  void setEmittanceY(double eps_y);
  double getPerveance();
  double getEmittanceX();
  double getEmittanceY();

private:
  double _perveance;     // beam perveance
  double _eps_x; // (4 * sqrt(<xx><x'x'> - <xx'><xx'>))
  double _eps_y; // (4 * sqrt(<yy><y'y'> - <yy'><yy'>))
};

#endif