#ifndef DANILOV_ENVELOPE_TRACKER_H
#define DANILOV_ENVELOPE_TRACKER_H

#include "Bunch.hh"
#include "CppPyWrapper.hh"

using namespace std;

/** Envelope solver for the Danilov distribution (KV distribution with
 zero emittance in one plane).

 The ellipse in the x-y plane can be parameterized as [1]:
     x = a * cos(psi) + b * sin(psi),
     y = e * cos(psi) + f * sin(psi),
 where 0 <= psi <= 2pi. The first two particles in the bunch are used to track
 the envelope parameters {a, b, e, f}, which are used to apply space charge
 kicks to the other particles in the bunch.

 References
 ----------
 [1] V. Danilov, S. Cousineau, S. Henderson, and J. Holmes, "Self-consistent time
     dependent two dimensional and three dimensional space charge distributions with
     linear force", PPRAB 6, 74–85 (2003).
*/
class DanilovEnvelopeTracker : public OrbitUtils::CppPyWrapper {
public:
  DanilovEnvelopeTracker(double perveanceQ);
  void trackBunch(Bunch *bunch, double length);
  void setPerveance(double perveance);
  double getPerveance();

private:
  double Q;  // beam perveance
};

#endif