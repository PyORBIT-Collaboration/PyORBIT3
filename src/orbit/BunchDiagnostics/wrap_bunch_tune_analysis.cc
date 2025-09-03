#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch_tune_analysis.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "BunchTuneAnalysis.hh"

namespace wrap_bunch_tune_analysis{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

/**
	Constructor for python class wrapping c++ BunchTuneAnalysis instance.
	It never will be called directly.
*/
static PyObject* BunchTuneAnalysis_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
	pyORBIT_Object* self;
	self = (pyORBIT_Object *) type->tp_alloc(type, 0);
	self->cpp_obj = NULL;
	return (PyObject *) self;
}

/** This is implementation of the __init__ method */
static int BunchTuneAnalysis_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	self->cpp_obj =  new BunchTuneAnalysis();
	((BunchTuneAnalysis*) self->cpp_obj)->setPyWrapper((PyObject*) self);
return 0;
}

/** Performs the Tune analysis of the bunch */
static PyObject* BunchTuneAnalysis_analyzeBunch(PyObject *self, PyObject *args){
	BunchTuneAnalysis* cpp_BunchTuneAnalysis = (BunchTuneAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	PyObject* pyBunch;
	if(!PyArg_ParseTuple(args,"O:analyzeBunch",&pyBunch)){
		ORBIT_MPI_Finalize("BunchTuneAnalysis - analyzeBunch(Bunch* bunch) - parameter are needed.");
	}
	PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
	if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
		ORBIT_MPI_Finalize("BunchTuneAnalysis - analyzeBunch(Bunch* bunch) - method needs a Bunch.");
	}
	Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
	cpp_BunchTuneAnalysis->analyzeBunch(cpp_bunch);
	Py_INCREF(Py_None);
	return Py_None;
}

/** Performs the Tune analysis of the bunch */
static PyObject* BunchTuneAnalysis_assignTwiss(PyObject *self, PyObject *args){
	BunchTuneAnalysis* cpp_BunchTuneAnalysis = (BunchTuneAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	double betax;
	double alphax;
	double etax;
	double etapx;
	double betay;
	double alphay;
	if(!PyArg_ParseTuple(args,"dddddd:assignTwiss",&betax, &alphax,&etax,&etapx,&betay,&alphay)){
		ORBIT_MPI_Finalize("BunchTuneAnalysis - getTwiss(double betax, double alphax, double etax, double etapx, double betay, double alphay) - parameter are needed.");
	}
	cpp_BunchTuneAnalysis->assignTwiss(betax, alphax, etax, etapx, betay, alphay);
	Py_INCREF(Py_None);
	return Py_None;
}



//--------------------------------------------------------------
//destructor for python BunchTuneAnalysis class (__del__ method).
//---------------------------------------------------------------
static void BunchTuneAnalysis_del(pyORBIT_Object* self){
	//std::cerr<<"The BunchTuneAnalysis __del__ has been called!"<<std::endl;
	delete ((BunchTuneAnalysis*)self->cpp_obj);
	self->ob_base.ob_type->tp_free((PyObject*)self);
}

// defenition of the methods of the python BunchTuneAnalysis wrapper class
// they will be vailable from python level
static PyMethodDef BunchTuneAnalysisClassMethods[] = {
	{ "analyzeBunch", BunchTuneAnalysis_analyzeBunch, METH_VARARGS,"Performs the Tune analysis of the bunch."},
	{ "assignTwiss", BunchTuneAnalysis_assignTwiss, METH_VARARGS,"Assigns Twiss at location of tune calculator."},
	{NULL}
};

// defenition of the memebers of the python BunchTwissAnalysis wrapper class
// they will be vailable from python level
static PyMemberDef BunchTuneAnalysisClassMembers [] = {
	{NULL}
};

//new python BunchTwissAnalysis wrapper type definition
static PyTypeObject pyORBIT_BunchTuneAnalysis_Type = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"BunchTuneAnalysis", /*tp_name*/
	sizeof(pyORBIT_Object), /*tp_basicsize*/
	0, /*tp_itemsize*/
	(destructor) BunchTuneAnalysis_del , /*tp_dealloc*/
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
	"The BunchTuneAnalysis python wrapper", /* tp_doc */
	0, /* tp_traverse */
	0, /* tp_clear */
	0, /* tp_richcompare */
	0, /* tp_weaklistoffset */
	0, /* tp_iter */
	0, /* tp_iternext */
	BunchTuneAnalysisClassMethods, /* tp_methods */
	BunchTuneAnalysisClassMembers, /* tp_members */
	0, /* tp_getset */
	0, /* tp_base */
	0, /* tp_dict */
	0, /* tp_descr_get */
	0, /* tp_descr_set */
	0, /* tp_dictoffset */
	(initproc) BunchTuneAnalysis_init, /* tp_init */
	0, /* tp_alloc */
	BunchTuneAnalysis_new, /* tp_new */
};


//--------------------------------------------------
//Initialization BunchTwissAnalysis of the pyBunchTwissAnalysis class
//--------------------------------------------------
void initbunchtuneanalysis(PyObject* module){
	if (PyType_Ready(&pyORBIT_BunchTuneAnalysis_Type) < 0) return;
	Py_INCREF(&pyORBIT_BunchTuneAnalysis_Type);
	PyModule_AddObject(module, "BunchTuneAnalysis", (PyObject *)&pyORBIT_BunchTuneAnalysis_Type);
}





//-----------------------------------------------------------------------------------
// BunchTwissAnalysis4D
//-----------------------------------------------------------------------------------


/**
Constructor for python class wrapping C++ BunchTuneAnalysis4D instance.
It never will be called directly.
*/
static PyObject* BunchTuneAnalysis4D_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
	pyORBIT_Object* self;
	self = (pyORBIT_Object *) type->tp_alloc(type, 0);
	self->cpp_obj = NULL;
	return (PyObject *) self;
}

/** Implementation of the __init__ method */
static int BunchTuneAnalysis4D_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	self->cpp_obj =  new BunchTuneAnalysis4D();
	((BunchTuneAnalysis4D*) self->cpp_obj)->setPyWrapper((PyObject*) self);
	return 0;
}

/** Analyzes the bunch. */
static PyObject* BunchTuneAnalysis4D_analyzeBunch(PyObject *self, PyObject *args){
	BunchTuneAnalysis4D* cpp_BunchTuneAnalysis4D = (BunchTuneAnalysis4D*)((pyORBIT_Object*) self)->cpp_obj;
	PyObject* pyBunch;
	if(!PyArg_ParseTuple(args,"O:analyzeBunch",&pyBunch)){
		ORBIT_MPI_Finalize("BunchTuneAnalysis4D - analyzeBunch(Bunch* bunch) - parameter are needed.");
	}
	PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
	if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
		ORBIT_MPI_Finalize("BunchTuneAnalysis4D - analyzeBunch(Bunch* bunch) - method needs a Bunch.");
	}
	Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
	cpp_BunchTuneAnalysis4D->analyzeBunch(cpp_bunch);
	Py_INCREF(Py_None);
	return Py_None;
}

/** Sets normalization matrix. */
// [TO DO] How to parse 2D Python list?
static PyObject* BunchTuneAnalysis4D_setMatrix(PyObject *self, PyObject *args){
	BunchTuneAnalysis4D* cpp_BunchTuneAnalysis4D = (BunchTuneAnalysis4D*)((pyORBIT_Object*) self)->cpp_obj;
	double dummy;
	if(!PyArg_ParseTuple(args,"d:setMatrix",&dummy)){
		ORBIT_MPI_Finalize("BunchTuneAnalysis4D - setMatrix(double dummy) - parameters are needed.");
	}
	cpp_BunchTuneAnalysis4D->setMatrix(dummy);
	Py_INCREF(Py_None);
	return Py_None;
}

// Destructor for python BunchTuneAnalysis4D class (__del__ method).
static void BunchTuneAnalysis4D_del(pyORBIT_Object* self){
	//std::cerr<<"The BunchTuneAnalysis4D __del__ has been called!"<<std::endl;
	delete ((BunchTuneAnalysis4D*)self->cpp_obj);
	self->ob_base.ob_type->tp_free((PyObject*)self);
}

// Definition of the methods of the python BunchTuneAnalysis4D wrapper class.
// They will be available from python level.
static PyMethodDef BunchTuneAnalysis4DClassMethods[] = {
	{ "analyzeBunch", BunchTuneAnalysis4D_analyzeBunch, METH_VARARGS,"Analyzes the bunch."},
	{ "setMatrix", BunchTuneAnalysis4D_setMatrix, METH_VARARGS,"Sets normalization matrix."},
	{NULL}
};

// Definition of the memebers of the python BunchTuneAnalysis4D wrapper class.
// They will be vailable from python level.
static PyMemberDef BunchTuneAnalysis4DClassMembers [] = {
	{NULL}
};

// New python BunchTuneAnalysis4D wrapper type definition.
static PyTypeObject pyORBIT_BunchTuneAnalysis4D_Type = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"BunchTuneAnalysis4D", /*tp_name*/
	sizeof(pyORBIT_Object), /*tp_basicsize*/
	0, /*tp_itemsize*/
	(destructor) BunchTuneAnalysis4D_del , /*tp_dealloc*/
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
	"The BunchTuneAnalysis4D python wrapper", /* tp_doc */
	0, /* tp_traverse */
	0, /* tp_clear */
	0, /* tp_richcompare */
	0, /* tp_weaklistoffset */
	0, /* tp_iter */
	0, /* tp_iternext */
	BunchTuneAnalysis4DClassMethods, /* tp_methods */
	BunchTuneAnalysis4DClassMembers, /* tp_members */
	0, /* tp_getset */
	0, /* tp_base */
	0, /* tp_dict */
	0, /* tp_descr_get */
	0, /* tp_descr_set */
	0, /* tp_dictoffset */
	(initproc) BunchTuneAnalysis4D_init, /* tp_init */
	0, /* tp_alloc */
	BunchTuneAnalysis4D_new, /* tp_new */
};

// Initialization of the pyBunchTuneAnalysis4D class.
void initbunchtuneanalysis4d(PyObject* module){
	if (PyType_Ready(&pyORBIT_BunchTuneAnalysis4D_Type) < 0) return;
	Py_INCREF(&pyORBIT_BunchTuneAnalysis4D_Type);
	PyModule_AddObject(module, "BunchTuneAnalysis4D", (PyObject *)&pyORBIT_BunchTuneAnalysis4D_Type);
}





#ifdef __cplusplus
}
#endif


}
