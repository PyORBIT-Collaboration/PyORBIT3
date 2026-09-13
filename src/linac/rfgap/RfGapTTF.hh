/**
   This class represents a Parmila type RF gap. It acts on the coordinates
   of the particle by using the transit time factors. The model includes
   non-linearity in transverse direction.
*/

#ifndef TTF_RF_GAP_H
#define TTF_RF_GAP_H

//MPI Function Wrappers
#include "mpi/orbit_mpi.hh"


#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "orbit/Bunch.hh"

//pyORBIT utils
#include "utils/polynomial/OU_Polynomial.hh"


/**
  This class represents a RF gap as a Parmila type gap.
*/

class RfGapTTF
{
public:

	/** Constructor for Parmila's type RF gap with TTF */
  RfGapTTF();

  /** Destructor */
  virtual ~RfGapTTF();

	/** Tracks the Bunch through the RF gap. */
	static void trackBunch(Bunch* bunch, double frequency, double E0L, double phase,
	                       OrbitUtils::Polynomial* Tttf,
	                       OrbitUtils::Polynomial* Sttf,
	                       OrbitUtils::Polynomial* Tpttf,
	                       OrbitUtils::Polynomial* Spttf);

};

#endif
