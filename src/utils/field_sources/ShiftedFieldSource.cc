//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   ShiftedFieldSource.cc
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    12/23/2020
//
// DESCRIPTION
//    A class is an implementation of BaseFiledSource class.
//    This class describe the instance of the Field Source (FS) class which is
//    shifted and rotated with respect of the center (0,0,0) original
//    coordinates.
//    The coordTransform Matrix is (4x4) matrix object that translates the
//    spacial coordinates from external coordinates to the inner coordinates
//    of this FS. The last raw in coordTransform Matrix is [0,0,0,1], and
//    The last column is a [a,b,c,1] where [-a,-b,-c] coordinates of this FS
//    origin in the external coordinates system.
//    The vectors (Ex,Ey,Ez) and (Bx,By,Bz) are transformed back to
//    the external coordinate system by coordTransformBack (3x3) matrix.
//
///////////////////////////////////////////////////////////////////////////
#include "mpi/orbit_mpi.hh"

#include "utils/field_sources/ShiftedFieldSource.hh"
#include "utils/matrix/MatrixOperations.hh"

using namespace OrbitUtils;

/** Constructor. */
ShiftedFieldSource::ShiftedFieldSource(): BaseFieldSource()
{

	ownedCoordTransformM4x4.reset(new Matrix(4,4));
	coordTransformM4x4 = ownedCoordTransformM4x4.get();
	coordTransformM3x3.reset(new Matrix(3,3));
	coordTransformM4x4->unit();
	coordTransformM3x3->unit();

	coordVectExt = new PhaseVector(4);
	coordVectInn = new PhaseVector(4);
	coordVectExt->zero();
	coordVectInn->zero();
	coordVectExt->getArray()[3] = 1.0;
	coordVectInn->getArray()[3] = 1.0;

	fieldVectExt = new PhaseVector(3);
	fieldVectInn = new PhaseVector(3);
	fieldVectExt->zero();
	fieldVectInn->zero();

}

/** Destructor */
ShiftedFieldSource::~ShiftedFieldSource()
{
	delete coordVectExt;
	delete coordVectInn;

	delete fieldVectExt;
	delete fieldVectInn;
}

/** Returns components of the electric and magnetic fields. */
void ShiftedFieldSource::getElectricMagneticField(
	  double x, double y, double z, double t,
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z){

		double* extArr = coordVectExt->getArray();
		extArr[0] = x;
		extArr[1] = y;
		extArr[2] = z;

		MatrixOperations::mult(coordTransformM4x4,coordVectExt,coordVectInn);

		double* innArr = coordVectInn->getArray();
		x = innArr[0];
		y = innArr[1];
		z = innArr[2];

		this->getInnerElectricMagneticField(x,y,z,t,E_x,E_y,E_z,H_x,H_y,H_z);

		//Electric field transformatiom

		innArr = fieldVectInn->getArray();
		innArr[0] = E_x;
		innArr[1] = E_y;
		innArr[2] = E_z;

		MatrixOperations::mult(coordTransformM3x3.get(),fieldVectInn,fieldVectExt);

		extArr = fieldVectExt->getArray();
		E_x = extArr[0];
		E_y = extArr[1];
		E_z = extArr[2];

		//Magnetic field transformatiom

		innArr = fieldVectInn->getArray();
		innArr[0] = H_x;
		innArr[1] = H_y;
		innArr[2] = H_z;

		MatrixOperations::mult(coordTransformM3x3.get(),fieldVectInn,fieldVectExt);

		extArr = fieldVectExt->getArray();
		H_x = extArr[0];
		H_y = extArr[1];
		H_z = extArr[2];

}

/** This method should be implemented in the sub-classes */
void ShiftedFieldSource::getInnerElectricMagneticField(
	  double x, double y, double z, double t,
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z){
	E_x=0.;	E_y=0.;	E_z=0.;
	H_x=0.;	H_y=0.;	H_z=0.;
}

/**
   Returns coordinates transformation matrix 4x4.
   It includes rotation 3x3 matrix and origin shift.
*/
Matrix* ShiftedFieldSource::getCoordsTransformMatrix()
{
	return coordTransformM4x4;
}

/**
   Sets coordinates transformation matrix 4x4.
   It includes rotation 3x3 matrix and origin shift.
*/
void ShiftedFieldSource::setCoordsTransformMatrix(Matrix* coordTransformM4x4_In)
{
	ownedCoordTransformM4x4.reset();
	coordTransformM4x4 = coordTransformM4x4_In;

	for(int ind_x = 0; ind_x < 3; ind_x++){
		for(int ind_y = 0; ind_y < 3; ind_y++){
			coordTransformM3x3->getArray()[ind_x][ind_y] = coordTransformM4x4->getArray()[ind_x][ind_y];
		}
	}
	int res = MatrixOperations::invert(coordTransformM3x3->getArray(),3);
	if(res == 0){
		ORBIT_MPI_Finalize("ShiftedFieldSource::setCoordsTransformMatrix - cannot invert coords. transformation matrix! Stop.");
	}
}
