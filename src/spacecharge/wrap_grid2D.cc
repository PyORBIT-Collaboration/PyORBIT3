#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_grid2D.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "Grid2D.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python Grid2D class definition
	//---------------------------------------------------------

	//constructor for python class wrapping Grid2D instance
	//It never will be called directly
	static PyObject* Grid2D_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The Grid2D new has been called!"<<std::endl;
		return (PyObject *) self;
	}

  //initializator for python  Grid2D class
  //this is implementation of the __init__ method
  static int Grid2D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
   int binX, binY;
	 double xMin = -1.0, yMin = -1.0 , xMax = +1.0 , yMax = +1.0;
	 if(!PyArg_ParseTuple(args,"ii|dddd:__init__",&binX,&binY,&xMin,&xMax,&yMin,&yMax)){
				ORBIT_MPI_Finalize("PyGrid2D - Grid2D(nX,nY[,xMin,xMax,yMin,yMax]) - constructor needs parameters.");
		}
		self->cpp_obj = new Grid2D(binX,binY,xMin,xMax,yMin,yMax);
		((Grid2D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//std::cerr<<"The Grid2D __init__ has been called!"<<std::endl;
		return 0;
  }

  //setZero()
  static PyObject* Grid2D_setZero(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		cpp_Grid2D->setZero();
		Py_INCREF(Py_None);
		return Py_None;
	}

	//getValue(double x, double y)
  static PyObject* Grid2D_getValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double x,y;
		if(!PyArg_ParseTuple(args,"dd:getValue",&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - getValue(x,y) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid2D->getValue(x,y));
	}

	//getValueBilinear(double x, double y)
  static PyObject* Grid2D_getValueBilinear(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double x,y;
		if(!PyArg_ParseTuple(args,"dd:getValue",&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - getValueBilinear(x,y) - parameters are needed.");
		}
		double val;
		cpp_Grid2D->interpolateBilinear(x,y,val);
		return Py_BuildValue("d",val);
	}

	//setValue(double value, int ix, int iy)
  static PyObject* Grid2D_setValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double val;
		int ix,iy;
		if(!PyArg_ParseTuple(args,"dii:setValue",&val,&ix,&iy)){
			ORBIT_MPI_Finalize("PyGrid2D - setValue(val,ix,iy) - parameters are needed.");
		}
		cpp_Grid2D->setValue(val,ix,iy);
		Py_INCREF(Py_None);
		return Py_None;
	}

	//getValueOnGrid(int ix, int iy)
  static PyObject* Grid2D_getValueOnGrid(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		int ix,iy;
		if(!PyArg_ParseTuple(args,"ii:getValueOnGrid",&ix,&iy)){
			ORBIT_MPI_Finalize("PyGrid2D - getValueOnGrid(ix,iy) - parameters are needed.");
		}
		return Py_BuildValue("d",cpp_Grid2D->getValueOnGrid(ix,iy));
	}

	//setGridX(double min, double max, int n)
  static PyObject* Grid2D_setGridX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double min,max;
		if(!PyArg_ParseTuple(args,"dd:setGridX",&min,&max)){
			ORBIT_MPI_Finalize("PyGrid2D - setGridX(min,max) - parameters are needed.");
		}
		cpp_Grid2D->setGridX(min,max);
		Py_INCREF(Py_None);
		return Py_None;
	}

	//setGridY(double min, double max, int n)
  static PyObject* Grid2D_setGridY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double min,max;
		if(!PyArg_ParseTuple(args,"dd:setGridY",&min,&max)){
			ORBIT_MPI_Finalize("PyGrid2D - setGridY(min,max) - parameters are needed.");
		}
		cpp_Grid2D->setGridY(min,max);
		Py_INCREF(Py_None);
		return Py_None;
	}

	//getGridX(ix)
  static PyObject* Grid2D_getGridX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(args,"i:getGridX",&ind) || ind < 0 || ind >= cpp_Grid2D->getSizeX()){
			ORBIT_MPI_Finalize("PyGrid2D - getGridX(ix) - parameter is needed. [0 - sizeX[");
		}
		return Py_BuildValue("d",cpp_Grid2D->getGridX(ind));
	}

	//getGridY(iy)
  static PyObject* Grid2D_getGridY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		int ind = -1;
		if(!PyArg_ParseTuple(args,"i:getGridY",&ind) || ind < 0 || ind >= cpp_Grid2D->getSizeY()){
			ORBIT_MPI_Finalize("PyGrid2D - getGridY(iy) - parameter is needed. [0 - sizeY[");
		}
		return Py_BuildValue("d",cpp_Grid2D->getGridY(ind));
	}

	//getSizeX()
  static PyObject* Grid2D_getSizeX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		return Py_BuildValue("i",cpp_Grid2D->getSizeX());
	}

	//getSizeY()
  static PyObject* Grid2D_getSizeY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		return Py_BuildValue("i",cpp_Grid2D->getSizeY());
	}

  //It will synchronize through the MPI communicator
  static PyObject* Grid2D_synchronizeMPI(PyObject *self, PyObject *args){
     pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		int nVars = PyTuple_Size(args);
		if(nVars == 0){
			cpp_Grid2D->synchronizeMPI(NULL);
		}
		else {
			PyObject* py_mpi_comm_type = wrap_orbit_mpi_comm::getMPI_CommType("MPI_Comm");
			PyObject* pyMPIComm = PyTuple_GetItem(args,0);
			if((!PyObject_IsInstance(pyMPIComm,py_mpi_comm_type))){
				ORBIT_MPI_Finalize("Grid2D.synchronizeMPI(MPI_Comm) - input parameter is not MPI_Comm");
			}
			cpp_Grid2D->synchronizeMPI((pyORBIT_MPI_Comm*) pyMPIComm);
		}
	 	Py_INCREF(Py_None);
		return Py_None;
  }


	//getMinX()
  static PyObject* Grid2D_getMinX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		return Py_BuildValue("d",cpp_Grid2D->getMinX());
	}

	//getMaxX()
  static PyObject* Grid2D_getMaxX(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		return Py_BuildValue("d",cpp_Grid2D->getMaxX());
	}

	//getMinY()
  static PyObject* Grid2D_getMinY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		return Py_BuildValue("d",cpp_Grid2D->getMinY());
	}

	//getMaxY()
  static PyObject* Grid2D_getMaxY(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		return Py_BuildValue("d",cpp_Grid2D->getMaxY());
	}

	//getSum()
  static PyObject* Grid2D_getSum(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		return Py_BuildValue("d",cpp_Grid2D->getSum());
	}

	//multiply()
  static PyObject* Grid2D_multiply(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double coeff;
		if(!PyArg_ParseTuple(args,"d:multiply",&coeff)){
			ORBIT_MPI_Finalize("PyGrid2D - multiply(coeff) - parameters are needed.");
		}
		cpp_Grid2D->multiply(coeff);
		Py_INCREF(Py_None);
    	return Py_None;
	}

	//isInside(x,y)
  static PyObject* Grid2D_isInside(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double x,y;
		if(!PyArg_ParseTuple(args,"dd:isInside",&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - isInside(x,y) - parameters are needed.");
		}
		return Py_BuildValue("i",cpp_Grid2D->isInside(x,y));
	}

	//binBunch(Bunch* bunch, [ind0,ind1]), by default ind0 = 0, ind1 = 2 (XY) plane
  static PyObject* Grid2D_binBunch(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		PyObject* pyBunch;
		int ind0 = -1;
		int ind1 = -1;
		if(!PyArg_ParseTuple(args,"O|ii:binBunch",&pyBunch,&ind0,&ind1)){
			ORBIT_MPI_Finalize("PyGrid2D - binBunch(Bunch* bunch, [ind0,ind1]) - parameter are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyGrid2D - binBunch(Bunch* bunch, [ind0,ind1]) - method needs a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		if(ind0 < 0 || ind1 < 0 || ind0 > 5 || ind1 > 5){
		  cpp_Grid2D->binBunch(cpp_bunch);
		} else {
			cpp_Grid2D->binBunch(cpp_bunch,ind0,ind1);
		}
		Py_INCREF(Py_None);
    return Py_None;
	}

	//binBunchBilinear(Bunch* bunch, [ind0,ind1]), by default ind0 = 0, ind1 = 2 (XY) plane
  static PyObject* Grid2D_binBunchBilinear(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		PyObject* pyBunch;
		int ind0 = -1;
		int ind1 = -1;
		if(!PyArg_ParseTuple(args,"O|ii:binBunchBilinear",&pyBunch,&ind0,&ind1)){
			ORBIT_MPI_Finalize("PyGrid2D - binBunchBilinear(Bunch* bunch, [ind0,ind1]) - parameter are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PyGrid2D - binBunchBilinear(Bunch* bunch, [ind0,ind1]) - method needs a Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		if(ind0 < 0 || ind1 < 0 || ind0 > 5 || ind1 > 5){
		  cpp_Grid2D->binBunchBilinear(cpp_bunch);
		} else {
			cpp_Grid2D->binBunchBilinear(cpp_bunch,ind0,ind1);
		}
		Py_INCREF(Py_None);
    return Py_None;
	}

	//binValue(double value, double x, double y)
  static PyObject* Grid2D_binValue(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double val,x,y;
		if(!PyArg_ParseTuple(args,"ddd:binValue",&val,&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - binValue(val,x,y) - parameters are needed.");
		}
		cpp_Grid2D->binValue(val,x,y);
		Py_INCREF(Py_None);
    return Py_None;
	}

	//binValueBilinear(double value, double x, double y)
  static PyObject* Grid2D_binValueBilinear(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double val,x,y;
		if(!PyArg_ParseTuple(args,"ddd:binValueBilinear",&val,&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - binValueBilinear(val,x,y) - parameters are needed.");
		}
		cpp_Grid2D->binValueBilinear(val,x,y);
		Py_INCREF(Py_None);
    return Py_None;
	}

	//calcGradient(double x, double y)
  static PyObject* Grid2D_calcGradient(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double x,y;
		double ex, ey;
		if(!PyArg_ParseTuple(args,"dd:calcGradient",&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - calcGradient(x,y) - parameters are needed.");
		}
		cpp_Grid2D->calcGradient(x,y,ex,ey);
		return Py_BuildValue("(dd)",ex,ey);
	}

	//calcGradientBilinear(double x, double y)
  static PyObject* Grid2D_calcGradientBilinear(PyObject *self, PyObject *args){
    pyORBIT_Object* pyGrid2D = (pyORBIT_Object*) self;
		Grid2D* cpp_Grid2D = (Grid2D*) pyGrid2D->cpp_obj;
		double x,y;
		double ex, ey;
		if(!PyArg_ParseTuple(args,"dd:calcGradientBilinear",&x,&y)){
			ORBIT_MPI_Finalize("PyGrid2D - calcGradientBilinear(x,y) - parameters are needed.");
		}
		cpp_Grid2D->calcGradientBilinear(x,y,ex,ey);
		return Py_BuildValue("(dd)",ex,ey);
	}


  //-----------------------------------------------------
  //destructor for python Grid2D class (__del__ method).
  //-----------------------------------------------------
  static void Grid2D_del(pyORBIT_Object* self){
		//std::cerr<<"The Grid2D __del__ has been called!"<<std::endl;
		Grid2D* cpp_Grid2D = (Grid2D*) self->cpp_obj;
		delete cpp_Grid2D;
		self->ob_base.ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python Grid2D wrapper class
	// they will be vailable from python level
  static PyMethodDef Grid2DClassMethods[] = {
		{ "setZero",              Grid2D_setZero,              METH_VARARGS,"sets all points on the grid to zero"},
		{ "getValue",             Grid2D_getValue,             METH_VARARGS,"returns value for (x,y) point calculated 9-points weighting scheme"},
		{ "getValueBilinear",     Grid2D_getValueBilinear,     METH_VARARGS,"returns value for (x,y) point using bi-linear interpolation"},
		{ "setValue",             Grid2D_setValue,             METH_VARARGS,"sets value for (ix,iy) point - (val,ix,iy)"},
		{ "getValueOnGrid",       Grid2D_getValueOnGrid,       METH_VARARGS,"returns value on the grid for (ind_x,ind_y) point"},
		{ "setGridX",             Grid2D_setGridX,             METH_VARARGS,"sets the X grid with min,max"},
		{ "setGridY",             Grid2D_setGridY,             METH_VARARGS,"sets the Y grid with min,max"},
		{ "getGridX",             Grid2D_getGridX,             METH_VARARGS,"returns the x-grid point with index ind"},
		{ "getGridY",             Grid2D_getGridY,             METH_VARARGS,"returns the x-grid point with index ind"},
		{ "getSizeX",             Grid2D_getSizeX,             METH_VARARGS,"returns the size of grid in X dir."},
		{ "getSizeY",             Grid2D_getSizeY,             METH_VARARGS,"returns the size of grid in Y dir."},
		{ "getMinX",              Grid2D_getMinX,              METH_VARARGS,"returns the min grid point in X dir."},
		{ "getMaxX",              Grid2D_getMaxX,              METH_VARARGS,"returns the max grid point in X dir."},
		{ "getMinY",              Grid2D_getMinY,              METH_VARARGS,"returns the min grid point in Y dir."},
		{ "getMaxY",              Grid2D_getMaxY,              METH_VARARGS,"returns the max grid point in Y dir."},
		{ "getSum",               Grid2D_getSum,               METH_VARARGS,"returns the sum of all grid points."},
		{ "multiply",             Grid2D_multiply,             METH_VARARGS,"multiplies all elements of Grid2D by constant coefficient"},
		{ "isInside",             Grid2D_isInside,             METH_VARARGS,"returns 1 or 0 if (x,y) inside grid or not"},
		{ "binValue",             Grid2D_binValue,             METH_VARARGS,"bins the value into the 2D mesh"},
		{ "binValueBilinear",     Grid2D_binValueBilinear,     METH_VARARGS,"bins the value into the 2D mesh bi-linearly"},
		{ "binBunch",             Grid2D_binBunch,             METH_VARARGS,"bins the Bunch instance into the 2D mesh (XY plane by default)"},
		{ "binBunchBilinear",     Grid2D_binBunchBilinear,     METH_VARARGS,"bins the Bunch instance into the 2D mesh bi-linearly (XY plane by default)"},
		{ "calcGradient",         Grid2D_calcGradient,         METH_VARARGS,"returns gradient as (gx,gy) for point (x,y) calculated by 9-points weighting scheme"},
		{ "calcGradientBilinear", Grid2D_calcGradientBilinear, METH_VARARGS,"returns gradient as (gx,gy) for point (x,y) calculated bi-linerly"},
		{ "synchronizeMPI",       Grid2D_synchronizeMPI,       METH_VARARGS,"synchronize through the MPI communicator"},
    {NULL}
  };

	// defenition of the memebers of the python Grid2D wrapper class
	// they will be vailable from python level
	static PyMemberDef Grid2DClassMembers [] = {
		{NULL}
	};

	//new python Grid2D wrapper type definition
	static PyTypeObject pyORBIT_Grid2D_Type = {
		PyVarObject_HEAD_INIT(NULL, 0)
		"Grid2D", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Grid2D_del , /*tp_dealloc*/
		0, /*tp_print*/
		0, /*tp_getattr*/
		0, /*tp_setattr*/
		0, /*tp_compare*/
		0, /*tp_repr*/
		0, /*tp_as_number*/
		0, /*tp_as_sequence*/
		0, /*tp_as_mapping*/
		0, /*tp_hash */
		0, /*tp_call*/
		0, /*tp_str*/
		0, /*tp_getattro*/
		0, /*tp_setattro*/
		0, /*tp_as_buffer*/
		Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
		"The Grid2D python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		Grid2DClassMethods, /* tp_methods */
		Grid2DClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Grid2D_init, /* tp_init */
		0, /* tp_alloc */
		Grid2D_new, /* tp_new */
	};

	//--------------------------------------------------
	//Initialization function of the pyGrid2D class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initGrid2D(PyObject* module){
		if (PyType_Ready(&pyORBIT_Grid2D_Type) < 0) return;
		Py_INCREF(&pyORBIT_Grid2D_Type);
		PyModule_AddObject(module, "Grid2D", (PyObject *)&pyORBIT_Grid2D_Type);
		//std::cout<<"debug Grid2D added! "<<std::endl;
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
