#ifndef KV_ENVELOPE_TRACKER_H
#define KV_ENVELOPE_TRACKER_H

#include "Bunch.hh"
#include "CppPyWrapper.hh"

using namespace std;

/** Envelope solver for upright KV distribution in uncoupled lattice.

The first particle in the bunch is used to track the envelope parameters.
The envelope parameters are use to apply space charge kicks to the other
particles in the bunch.

This class does not yet handle boudary conditions or dispersion.
 */
class KVEnvelopeTracker : public OrbitUtils::CppPyWrapper {
  public:
    KVEnvelopeTracker(double perveance, double emittance_x, double emittance_y);
    void trackBunch(Bunch *bunch, double length);
    void setPerveance(double perveance);
    void setEmittanceX(double emittance);
    void setEmittanceY(double emittance);
    double getPerveance();
    double getEmittanceX();
    double getEmittanceY();

  private:
    double Q;     // beam perveance
    double eps_x; // (4 * sqrt(<xx><x'x'> - <xx'><xx'>))
    double eps_y; // (4 * sqrt(<yy><y'y'> - <yy'><yy'>))
};

#endif