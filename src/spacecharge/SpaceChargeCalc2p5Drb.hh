/**
 This class calculates the space charge kicks for bunch using 2.5D approach in transverse direction
 and Rick Baartman's approach (RB) to the longitudinal kicks. It is not sutable for short bunches, but can be used in
 rings.
*/

#ifndef SC_SPACECHARGE_CALC_2P5D_RB_H
#define SC_SPACECHARGE_CALC_2P5D_RB_H

//MPI Function Wrappers
#include "mpi/orbit_mpi.hh"

#include <cstdlib>
#include <cmath>
#include <memory>

//ORBIT bunch
#include "orbit/Bunch.hh"

//pyORBIT utils
#include "utils/bunch/BunchExtremaCalculator.hh"

#include "spacecharge/Grid1D.hh"
#include "spacecharge/Grid2D.hh"
#include "spacecharge/PoissonSolverFFT2D.hh"


class SpaceChargeCalc2p5Drb
{
public:

	/** Constructor with the "x to y ratio" parameter. */
	SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize, double xy_ratio_in);

	/** Constructor with "ratio" parameter equals 1. */
	SpaceChargeCalc2p5Drb(int xSize, int ySize, int zSize);

	/** Destructor */
	virtual ~SpaceChargeCalc2p5Drb();

	/** Calculates space charge and applies the transverse and
	    longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length, double pipe_radius);

	/** Returns the 2D rho grid with a transverse density distribution. **/
	Grid2D* getRhoGrid();

	/** Returns the 2D phi grid with a transverse potential. **/
	Grid2D* getPhiGrid();

	/** Returns the 1D grid with a longitudinal density. **/
	Grid1D* getLongGrid();

	/** Returns the 1D grid with a derivative of the longitudinal density. **/
	Grid1D* getLongDerivativeGrid();

	/** Sets the number of smoothing points to calculate the derivative of the longitudinal density. */
	void setLongAveragingPointsN(int n_points);

	/** Returns the number of smoothing points to calculate the derivative of the longitudinal density. */
	int getLongAveragingPointsN();

private:

	/** Analyses the bunch and does bining. */
 void bunchAnalysis(Bunch* bunch, double& totalMacrosize, double& x_c, double& y_c, double& a_bunch);

 /** Calculates the derivative of the longitudinal density by using Quadratic Curve Fitting */
 void calculateLongDerivative();

protected:
	PoissonSolverFFT2D* poissonSolver;
	std::unique_ptr<Grid2D> rhoGrid;
	std::unique_ptr<Grid2D> phiGrid;
	std::unique_ptr<Grid1D> zGrid;
	std::unique_ptr<Grid1D> zDerivGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;

	double xy_ratio;
	int n_long_avg;

	//auxiliary 5x2 array for Quadratic Curve Fitting
	double** S_arr;
};
//end of SC_SPACECHARGE_CALC_2P5D_RB_H
#endif
