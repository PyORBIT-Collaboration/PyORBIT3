//Calculate the space charge effect of the bunch in the 2.5D

#ifndef SC_SPACEFORCECHARGE_CALC_2P5D_H
#define SC_SPACEFORCECHARGE_CALC_2P5D_H

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
#include "spacecharge/ForceSolverFFT2D.hh"
#include "spacecharge/BaseBoundary2D.hh"


class SpaceChargeForceCalc2p5D
{
public:

	/** Constructor */
	SpaceChargeForceCalc2p5D(int xSize, int ySize, int zSize, double xy_ratio_in);

	SpaceChargeForceCalc2p5D(int xSize, int ySize, int zSize);

	/** Destructor */
	virtual ~SpaceChargeForceCalc2p5D();

	/** Calculates space charge and applies the transverse and
	longitudinal SC kicks to the macro-particles in the bunch. */
	void trackBunch(Bunch* bunch, double length);

	/** Returns the 2D rho grid with a transverse density distribution. **/
	Grid2D* getRhoGrid();

	/** Returns the 2D horizontal force grid with a transverse force. **/
	Grid2D* getForceGridX();

	/** Returns the 2D vertical force grid with a transverse force. **/
	Grid2D* getForceGridY();

	/** Returns the 1D grid with a longitudinal density. **/
	Grid1D* getLongGrid();

private:
	/** Analyses the bunch and does bining. */
 void bunchAnalysis(Bunch* bunch, double& totalMacrosize);

protected:
	ForceSolverFFT2D* forceSolver;
	std::unique_ptr<Grid2D> rhoGrid;
	std::unique_ptr<Grid2D> phiGrid;
	std::unique_ptr<Grid2D> forceGridX;
	std::unique_ptr<Grid2D> forceGridY;
	std::unique_ptr<Grid1D> zGrid;
	OrbitUtils::BunchExtremaCalculator* bunchExtremaCalc;

};
//end of SC_SPACEFORCECHARGE_CALC_2P5D_H
#endif
