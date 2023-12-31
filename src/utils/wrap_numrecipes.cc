#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"

#include <iostream>
#include <string>

#include "bessel.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_numrecipes{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	static PyObject* numrecipes_bessj0(PyObject* self, PyObject* args)
	{
		double x;
		if(!PyArg_ParseTuple(	args,"d:bessj0",&x)){
			error("bessj0(x) - parameter is needed.");
		}
		return Py_BuildValue("d",bessj0(x));
	}

	static PyObject* numrecipes_bessj1(PyObject* self, PyObject* args)
	{
		double x;
		if(!PyArg_ParseTuple(	args,"d:bessj1",&x)){
			error("bessj1(x) - parameter is needed.");
		}
		return Py_BuildValue("d",bessj1(x));
	}

	static PyObject* numrecipes_bessj(PyObject* self, PyObject* args)
	{
		double x;
		int n;
		if(!PyArg_ParseTuple(	args,"id:bessj_n_x",&n,&x)){
			error("bessj(n,x) - parameters are needed.");
		}
		return Py_BuildValue("d",bessj(n,x));
	}


	static PyObject* numrecipes_bessi0(PyObject* self, PyObject* args)
	{
		double x;
		if(!PyArg_ParseTuple(	args,"d:bessi0",&x)){
			error("bessi0(x) - parameter is needed.");
		}
		return Py_BuildValue("d",bessi0(x));
	}

	static PyObject* numrecipes_bessi1(PyObject* self, PyObject* args)
	{
		double x;
		if(!PyArg_ParseTuple(	args,"d:bessi1",&x)){
			error("bessi1(x) - parameter is needed.");
		}
		return Py_BuildValue("d",bessi1(x));
	}

	static PyObject* numrecipes_bessi(PyObject* self, PyObject* args)
	{
		double x;
		int n;
		if(!PyArg_ParseTuple(	args,"id:bessi_n_x",&n,&x)){
			error("bessi(n,x) - parameters are needed.");
		}
		return Py_BuildValue("d",bessi(n,x));
	}

	static PyMethodDef NumrecipesModuleMethods[] = {
		{"bessi0", numrecipes_bessi0 , METH_VARARGS, "bessi0(x) - I0(x) First kind modified Bessel Function."},
		{"bessi1", numrecipes_bessi1 , METH_VARARGS, "bessi1(x) - I1(x) First kind modified Bessel Function."},
		{"bessi",  numrecipes_bessi  , METH_VARARGS, "bessi(n,x) - I(n,x) First kind modified Bessel Function."},
		{"bessj0", numrecipes_bessj0 , METH_VARARGS, "bessj0(x) -  J0(x) First kind Bessel Function."},
		{"bessj1", numrecipes_bessj1 , METH_VARARGS, "bessj1(x) -  J1(x) First kind Bessel Function."},
		{"bessj",  numrecipes_bessj  , METH_VARARGS, "bessj(n,x) - J(n,x) First kind Bessel Function."},
		{NULL, NULL, 0, NULL}        /* Sentinel */
	};


	static struct PyModuleDef cModPyDem =
	{
		PyModuleDef_HEAD_INIT,
		"num_recipes", "Bessel functions implementations from Numerical Recipes.",
		-1,
		NumrecipesModuleMethods
	};

	//--------------------------------------------------
	//Initialization functions of the numrecipes module
	//--------------------------------------------------
  void initNumrecipes(PyObject* module){
    //create numrecipes module
    PyObject* module_nr = PyModule_Create(&cModPyDem);
		Py_INCREF(module_nr);
		PyModule_AddObject(module,const_cast<char*>("num_recipes"), module_nr);
	}

#ifdef __cplusplus
}
#endif


}
