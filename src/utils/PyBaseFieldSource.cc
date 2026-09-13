//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyBaseFieldSource.cc
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source.
//    It should be sub-classed on Python level and implements
//    getElectricMagneticField(x,y,z,t) method (returns (Ex,Ey,Ez,Bx,By,Bz)).
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise!
//               It is slow, because it delegates the field calculations
//               to the Python level. It is only for prototyping and
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#include "utils/PyBaseFieldSource.hh"

#include "mpi/orbit_mpi.hh"
#include <iostream>

using namespace OrbitUtils;

PyBaseFieldSource::PyBaseFieldSource(PyObject* py_wrapperIn)
{
	py_wrapper = py_wrapperIn;
}

PyBaseFieldSource::~PyBaseFieldSource()
{
}

void PyBaseFieldSource::getElectricMagneticField(
	double x, double y, double z, double t,
	double& fe_x, double& fe_y, double& fe_z,
	double& fm_x, double& fm_y, double& fm_z)
{
	PyObject* ef_tuple = PyObject_CallMethod(py_wrapper,const_cast<char*>("getElectricMagneticField"),const_cast<char*>("dddd"),x,y,z,t);
    if(!PyArg_ParseTuple(	ef_tuple,"dddddd:electric_magnetic_field",&fe_x,&fe_y,&fe_z,&fm_x,&fm_y,&fm_z)){
      ORBIT_MPI_Finalize("PyBaseFieldSource - getElectricMagneticField(x,y,z,t0 method does not work!");
    }
		//the references should be decreased because they were created as "new reference"
		Py_DECREF(ef_tuple);
}
