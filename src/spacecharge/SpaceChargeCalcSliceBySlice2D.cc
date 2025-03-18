/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalcSliceBySlice2D.cc
//
//   04/06/2014 H. Bartosik
//
// DESCRIPTION
//   Calculate the space charge effect of the bunch in 2D slice by slice
//
/////////////////////////////////////////////////////////////////////////////

#include "Grid3D.hh"
#include "PoissonSolverFFT2D.hh"
#include "SpaceChargeCalcSliceBySlice2D.hh"
#include "BufferStore.hh"

#include <iostream>
#include <cmath>
#include <cfloat>

using namespace OrbitUtils;

SpaceChargeCalcSliceBySlice2D::SpaceChargeCalcSliceBySlice2D(int xSize, int ySize, int zSize, double xy_ratio_in): CppPyWrapper(NULL)
{
	xy_ratio = xy_ratio_in;
	useLongTracking = 0;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -xy_ratio, xy_ratio, -1.0, 1.0);
	rhoGrid3D = new Grid3D(xSize, ySize, zSize);
	phiGrid3D = new Grid3D(xSize, ySize, zSize);
	bunchExtremaCalc = new BunchExtremaCalculator();
}

SpaceChargeCalcSliceBySlice2D::SpaceChargeCalcSliceBySlice2D(int xSize, int ySize, int zSize): CppPyWrapper(NULL)
{
	xy_ratio = 1.0;
	useLongTracking = 0;
	poissonSolver = new PoissonSolverFFT2D(xSize, ySize, -xy_ratio, xy_ratio, -1.0, 1.0);
	rhoGrid3D = new Grid3D(xSize, ySize, zSize);
	phiGrid3D = new Grid3D(xSize, ySize, zSize);
	bunchExtremaCalc = new BunchExtremaCalculator();
}

SpaceChargeCalcSliceBySlice2D::~SpaceChargeCalcSliceBySlice2D(){
	delete poissonSolver;
	if(rhoGrid3D->getPyWrapper() != NULL){
		Py_DECREF(rhoGrid3D->getPyWrapper());
	} else {
		delete rhoGrid3D;
	}
	if(phiGrid3D->getPyWrapper() != NULL){
		Py_DECREF(phiGrid3D->getPyWrapper());
	} else {
		delete phiGrid3D;
	}
	delete bunchExtremaCalc;
}

void SpaceChargeCalcSliceBySlice2D::longTracking(int useLongTracking){
	this->useLongTracking = useLongTracking;
}

int SpaceChargeCalcSliceBySlice2D::getLongitudinalTracking()
{
	return useLongTracking;
}

Grid3D* SpaceChargeCalcSliceBySlice2D::getRhoGrid(){
	return rhoGrid3D;
}

Grid3D* SpaceChargeCalcSliceBySlice2D::getPhiGrid(){
	return phiGrid3D;
}


void SpaceChargeCalcSliceBySlice2D::trackBunch(Bunch* bunch, double length, BaseBoundary2D* boundary){

	int nPartsGlobal = bunch->getSizeGlobal();
	if(nPartsGlobal < 2) return;

	double totalMacrosize = 0.;
	this->bunchAnalysis(bunch, totalMacrosize, boundary);
	//std::cerr<<"totalMacrosize="<<totalMacrosize;

	int rank = 0;
	int size = 0;
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size);
	int nZ = rhoGrid3D->getSizeZ();

	for(int iz = 0; iz < nZ; iz++){
		if(rank == iz%size){
			// calculate the potential on the temporary 2D phiGrid
			poissonSolver->findPotential(rhoGrid3D->getGrid2D(iz), phiGrid3D->getGrid2D(iz));
			if(boundary != NULL){
				//update potential with boundary condition
				boundary->addBoundaryPotential(rhoGrid3D->getGrid2D(iz), phiGrid3D->getGrid2D(iz));
				//std::cerr<<"Boundary ADDED."<<std::endl;
			}
		}
		else{
			phiGrid3D->getGrid2D(iz)->setZero();
		}
	}

	// synchronize the 3D phiGrid
	phiGrid3D->synchronizeMPI(bunch->getMPI_Comm_Local());
	
	SyncPart* syncPart = bunch->getSyncPart();
	
	double beta = syncPart->getBeta();
	double gamma = syncPart->getGamma();

	double trans_factor = length*bunch->getClassicalRadius()/(pow(beta,2)*pow(gamma,2));
	double long_factor = length*bunch->getClassicalRadius()*bunch->getMass();

	double x,y,z,ex,ey,ez;
	
	// defines if we use longitudinal electric field in particles tracking
	int long_tracking = useLongTracking;
	if(boundary == NULL){
		long_tracking = 0;
	}

	for (int i = 0, n = bunch->getSize(); i < n; i++){
		x = bunch->x(i);
		y = bunch->y(i);
		z = bunch->z(i)*gamma;

		if(boundary == NULL || (boundary != NULL && boundary->isInside(x,y) == BaseBoundary2D::IS_INSIDE)){
			phiGrid3D->calcGradient(x,ex,y,ey,z,ez);

			bunch->xp(i) -= ex * trans_factor;
			bunch->yp(i) -= ey * trans_factor;
			
			if(long_tracking == 1){
				bunch->dE(i) -= ez * long_factor;
			}
		}
	}
}

void SpaceChargeCalcSliceBySlice2D::bunchAnalysis(Bunch* bunch, double& totalMacrosize, BaseBoundary2D* boundary){

	double xMin, xMax, yMin, yMax, zMin, zMax;

	if(boundary == NULL){

		bunchExtremaCalc->getExtremaXYZ(bunch, xMin, xMax, yMin, yMax, zMin, zMax);

		//check if the beam size is not zero
		if( xMin >=  xMax || yMin >=  yMax || zMin >=  zMax){
			int rank = 0;
			ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if(rank == 0){
				std::cerr << "SpaceChargeCalcSliceBySlice2D::bunchAnalysis(bunch,...)" << std::endl
				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
				<< "x min ="<< xMin <<" max="<< xMax << std::endl
				<< "y min ="<< yMin <<" max="<< yMax << std::endl
				<< "z min ="<< zMin <<" max="<< zMax << std::endl
				<< "Stop."<< std::endl;
			}
			ORBIT_MPI_Finalize();
		}

		/** Here we will define the x/y ratio to use it in the FFT solver.
		    The same ratio means that we do not have to recalculate Green
				functions for FFT Possion solver. So, if the existing ratio changed
				less for 25% we will fix the space charge region which is faster
				than the Green functions recalculations in poissonSolver->setGridXY(...)
				method.
		*/
		double xy_ratio_beam = (xMax - xMin)/(yMax - yMin);
		if(fabs(xy_ratio_beam-xy_ratio) < 0.25*xy_ratio){
			double width, center;
			if(xy_ratio_beam > xy_ratio){
				center = (yMax + yMin)/2.0;
				width = ((yMax - yMin)*(xy_ratio_beam/xy_ratio))/2.0;
				yMin = center - width;
				yMax = center + width;
			} else {
				center = (xMax + xMin)/2.0;
				width = ((xMax - xMin)/(xy_ratio_beam/xy_ratio))/2.0;
				xMin = center - width;
				xMax = center + width;
			}
			xy_ratio = (xMax - xMin)/(yMax - yMin);
		}
		else{
			xy_ratio = xy_ratio_beam;
			poissonSolver->setGridXY(xMin,xMax,yMin,yMax);
			//std::cerr << "debug v0 grid changed r="<<xy_ratio<< std::endl;
		}
	}
	else{
		xMax = boundary->getMaxX();
		xMin = boundary->getMinX();
		yMax = boundary->getMaxY();
		yMin = boundary->getMinY();

		xy_ratio = (xMax - xMin)/(yMax - yMin);

		bunchExtremaCalc->getExtremaZ(bunch, zMin, zMax);

		//check if the beam size is not zero
		if(zMin >=  zMax){
			int rank = 0;
			ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if(rank == 0){
				std::cerr << "SpaceChargeCalcSliceBySlice2D::bunchAnalysis(bunch,...)" << std::endl
				<< "The bunch min and max sizes are wrong! Cannot calculate space charge!" << std::endl
				<< "z min ="<< zMin <<" max="<< zMax << std::endl
				<< "Stop."<< std::endl;
			}
			ORBIT_MPI_Finalize();
		}
	}

	//set Grids' limits
	rhoGrid3D->setGridX(xMin,xMax);
	rhoGrid3D->setGridY(yMin,yMax);
	rhoGrid3D->setGridZ(zMin,zMax);
	phiGrid3D->setGridX(xMin,xMax);
	phiGrid3D->setGridY(yMin,yMax);
	phiGrid3D->setGridZ(zMin,zMax);


	//this one just for case boundary != null, and will work only once
	double solver_xMin = poissonSolver->getMinX();
	double solver_xMax = poissonSolver->getMaxX();
	double solver_yMin = poissonSolver->getMinY();
	double solver_yMax = poissonSolver->getMaxY();
	double shape_diff_limit = 0.00000001;
	if(fabs((solver_xMax-solver_xMin)/(solver_yMax-solver_yMin)- xy_ratio) > shape_diff_limit){
		poissonSolver->setGridXY(xMin,xMax,yMin,yMax);
		//std::cerr << "debug v1 grid changed r="<<xy_ratio<< std::endl;
	}

	//bin Bunch to the Grid
	rhoGrid3D->setZero();
	rhoGrid3D->binBunchSlice2D(bunch);
	// Before we solve 2D Poisson problem for 2D potentials 
	// we have to redefine the rho values because for 2D
	// problem we will treat them as linear density not like
	// point-like charges.
	// Also we have a factor 2 missing from from our Poisson 2D solver
	// because we use lambda/r_transv instead of 2*lambda/r_transv.
	// We also have to take into account
	// the Lorentz length transformation from moving bunch length
	// to the synch. particle rest system of coordinates.
	// N.B. - after this the sum of all elements will not be equal to
	// the total macro-size of the bunch
	double gamma = bunch->getSyncPart()->getGamma();
	zMin *= gamma;
	zMax *= gamma;
	rhoGrid3D->setGridZ(zMin,zMax);
	phiGrid3D->setGridZ(zMin,zMax);
	rhoGrid3D->multiply(2./rhoGrid3D->getStepZ());
	//-----------------------------
	rhoGrid3D->synchronizeMPI(bunch->getMPI_Comm_Local());
}
