// This class represents a simple RF gap.
// For this RF gap we know the E0TL parameter only.

#ifndef BASE_RF_GAP_H
#define BASE_RF_GAP_H

//MPI Function Wrappers
#include "mpi/orbit_mpi.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "orbit/Bunch.hh"

//pyORBIT utils


/**
  This class represents a 2D rectangular grid.
*/

class BaseRfGap
{
public:

  /** Constructor for Base RF gap*/
  BaseRfGap();

  /** Destructor */
  virtual ~BaseRfGap();

  /** Tracks the Bunch trough the RF gap. */
  static void trackBunch(Bunch* bunch, double frequency, double E0TL, double phase);

};

#endif
