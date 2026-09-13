/**
This class represents a simplified RF gap. It acts on the coordinates like a transport matrix.
There are no nonlinear effects. It should be analog of RF Gap of XAL online model or Trace3D.
For this RF gap we know the E0TL, frequency, and phase only.
*/

#ifndef MATRIX_RF_GAP_H
#define MATRIX_RF_GAP_H

//MPI Function Wrappers
#include "mpi/orbit_mpi.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "orbit/Bunch.hh"

//pyORBIT utils


/**
  This class represents a RF gap as transport matrix. No nonlinear effects.
*/

class MatrixRfGap
{
public:

	/** Constructor for Base RF gap*/
  MatrixRfGap();

  /** Destructor */
  virtual ~MatrixRfGap();

	/** Tracks the Bunch trough the RF gap. */
	static void trackBunch(Bunch* bunch, double frequency, double E0TL, double phase);

};

#endif
