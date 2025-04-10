//Calculate the space charge effect of the bunch in 2D slice by slice

#ifndef SC_SPACECHARGE_CALC_2P5D_H
#define SC_SPACECHARGE_CALC_2P5D_H

//MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "Bunch.hh"

//pyORBIT utils
#include "CppPyWrapper.hh"
#include "BunchExtremaCalculator.hh"

#include "Grid3D.hh"
#include "PoissonSolverFFT2D.hh"
#include "BaseBoundary2D.hh"

using namespace std;

class SpaceChargeCalcSliceBySlice2D: public OrbitUtils::CppPyWrapper
{
public:

	/** Constructor */
	SpaceChargeCalcSliceBySlice2D(int xSize, int ySize, int zSize, double xy_ratio_in);

	SpaceChargeCalcSliceBySlice2D(int xSize, int ySize, int zSize);

	/** Destructor */
	virtual ~SpaceChargeCalcSliceBySlice2D();

	/** Calculates space charge and applies the transverse
	SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary);

	/** Returns the 3D rho grid with a transverse density distribution. */
	Grid3D* getRhoGrid();

	/** Returns the 3D phi grid with a transverse potential. */
	Grid3D* getPhiGrid();
	
	/** 
	    Sets if we will use the longitudinal dynamics for bunch tracking.
	    Remember if there is no boundary this settings is irrelevant,
	    we cannot perform longitudinal electric field from the Space Charge 
	    Solver. useLongTracking = 1 - will use, 0 - no. 
	*/
	void longTracking(int useLongTracking);
	int getLongitudinalTracking();

private:
	/** Analyses the bunch and does bining. */
	void bunchAnalysis(Bunch* bunch, double& totalMacrosize, BaseBoundary2D* boundary);
	int useLongTracking;	


protected:
	PoissonSolverFFT2D* poissonSolver;
	Grid3D* rhoGrid3D;
	Grid3D* phiGrid3D;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;

	double xy_ratio;
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
