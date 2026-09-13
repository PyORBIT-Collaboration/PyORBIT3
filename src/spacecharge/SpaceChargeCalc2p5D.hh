//Calculate the space charge effect of the bunch in the 2.5D

#ifndef SC_SPACECHARGE_CALC_2P5D_H
#define SC_SPACECHARGE_CALC_2P5D_H

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
#include "spacecharge/BaseBoundary2D.hh"


class SpaceChargeCalc2p5D
{
public:

	/** Constructor */
	SpaceChargeCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio_in);

	SpaceChargeCalc2p5D(int xSize, int ySize, int zSize);

	/** Destructor */
	virtual ~SpaceChargeCalc2p5D();

	/** Calculates space charge and applies the transverse and
	longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary);

	/** Returns the 2D rho grid with a transverse density distribution. **/
	Grid2D* getRhoGrid();

	/** Returns the 2D phi grid with a transverse potential. **/
	Grid2D* getPhiGrid();

	/** Returns the 1D grid with a longitudinal density. **/
	Grid1D* getLongGrid();

private:
	/** Analyses the bunch and does bining. */
 void bunchAnalysis(Bunch* bunch, double& totalMacrosize, BaseBoundary2D* boundary);

protected:
	PoissonSolverFFT2D* poissonSolver;
	std::unique_ptr<Grid2D> rhoGrid;
	std::unique_ptr<Grid2D> phiGrid;
	std::unique_ptr<Grid1D> zGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;

	double xy_ratio;
};
//end of SC_SPACECHARGE_CALC_2P5D_H
#endif
