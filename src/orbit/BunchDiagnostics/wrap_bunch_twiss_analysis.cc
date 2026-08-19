#include "orbit/BunchDiagnostics/wrap_bunch_twiss_analysis.hh"

#include "main/pyORBIT_Object.hh"
#include "mpi/orbit_mpi.hh"
#include "orbit/BunchDiagnostics/BunchTwissAnalysis.hh"
#include "orbit/wrap_bunch.hh"

namespace wrap_bunch_twiss_analysis
{

void error(const char* msg)
{
  ORBIT_MPI_Finalize(msg);
}

#ifdef __cplusplus
extern "C" {
#endif

/**
    Constructor for python class wrapping c++ BunchTwissAnalysis instance.
    It never will be called directly.
*/
static PyObject* BunchTwissAnalysis_new(PyTypeObject* type, PyObject* Py_UNUSED(args), PyObject* Py_UNUSED(kwds))
{
  pyORBIT_Object* self;
  self = (pyORBIT_Object*)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject*)self;
}

/** This is implementation of the __init__ method */
static int BunchTwissAnalysis_init(pyORBIT_Object* self, PyObject* Py_UNUSED(args), PyObject* Py_UNUSED(kwds))
{
  self->cpp_obj = new BunchTwissAnalysis();
  ((BunchTwissAnalysis*)self->cpp_obj)->setPyWrapper((PyObject*)self);
  return 0;
}

/** Performs the Twiss analysis of the bunch */
static PyObject* BunchTwissAnalysis_analyzeBunch(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(arg, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("BunchTwissAnalysis - analyzeBunch(Bunch* bunch) - method needs a Bunch.");
  }

  Bunch* cpp_bunch = (Bunch*)((pyORBIT_Object*)arg)->cpp_obj;
  cpp_BunchTwissAnalysis->analyzeBunch(cpp_bunch);

  Py_RETURN_NONE;
}

/** Returns the XY moments of the beam up to a prescribed order */
static PyObject* BunchTwissAnalysis_computeBunchMoments(PyObject* self, PyObject* args)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;
  PyObject* pyBunch;
  int order = 2;
  int normalize = 0;
  int emitnormterm = 0;
  int dispterm = 0;
  if (!PyArg_ParseTuple(
        args,
        "O|ippp:computeBunchMoments",
        &pyBunch,
        &order,
        &normalize,
        &emitnormterm,
        &dispterm
      )) {
    ORBIT_MPI_Finalize(
      "BunchTwissAnalysis - computeBunchMoments(Bunch* bunch, int order, int dispterm, int "
      "emitnormterm) - parameters are needed."
    );
  }
  PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize(
      "BunchTwissAnalysis - computeBunchMoments(Bunch* bunch, int order, int dispterm, int "
      "emitnormterm) - method needs a Bunch."
    );
  }

  Bunch* cpp_bunch = (Bunch*)((pyORBIT_Object*)pyBunch)->cpp_obj;
  cpp_BunchTwissAnalysis->computeBunchMoments(cpp_bunch, order, normalize,  emitnormterm, dispterm);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* BunchTwissAnalysis_getCovariance(PyObject* self, PyObject* args)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;
  int i, j;
  if (!PyArg_ParseTuple(args, "ii:getCovariance", &i, &j)) {
    error("pyBunchTwissAnalysis.getCovariance(i, j) - parameters are needed");
  }
  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getCovariance(i, j));
}

static PyObject* BunchTwissAnalysis_getCorrelation(PyObject* self, PyObject* args)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;
  int i, j;
  if (!PyArg_ParseTuple(args, "ii:getCorrelation", &i, &j)) {
    error("pyBunchTwissAnalysis.getCorrelation(i, j) - parameters are needed");
  }
  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getCorrelation(i, j));
}

/** It will return the (i,j) XY moment of the beam */
static PyObject* BunchTwissAnalysis_getBunchMoment(PyObject* self, PyObject* args)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;
  int i, j;
  if (!PyArg_ParseTuple(args, "ii:getBunchMoment", &i, &j)) {
    error("pyBunchTwissAnalysis.getBunchMoment(i,j) - parameters are needed");
  }
  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getBunchMoment(i, j));
}

/** Returns the average value for coordinate with index ic */
static PyObject* BunchTwissAnalysis_getAverage(PyObject* self, PyObject* args)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;
  int ic;
  if (!PyArg_ParseTuple(args, "i:", &ic)) {
    error("pyBunchTwissAnalysis.getAverage(ic) - parameter is needed");
  }
  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getAverage(ic));
}

/** Returns the total number of analysed macroparticles */
static PyObject* BunchTwissAnalysis_getGlobalCount(PyObject* self, PyObject* Py_UNUSED(ignored))
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;
  return Py_BuildValue("i", cpp_BunchTwissAnalysis->getGlobalCount());
}

/** Returns the total macrosize */
static PyObject* BunchTwissAnalysis_getGlobalMacrosize(PyObject* self, PyObject* Py_UNUSED(ignored))
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;
  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getGlobalMacrosize());
}

//------------------------------------------------------------
// Twiss functions
//------------------------------------------------------------
/** It returns the emittance for index 0,1,2 - x,y,z planes*/
static PyObject* BunchTwissAnalysis_getEmittance(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getEmittance", &ic)) {
    error("pyBunchTwissAnalysis.getEmittance(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getEmittance(ic));
}

/** It returns the normalized emittance for index 0,1 - x,y planes*/
static PyObject* BunchTwissAnalysis_getEmittanceNormalized(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getEmittanceNormalized", &ic)) {
    error("pyBunchTwissAnalysis.getEmittanceNormalized(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getEmittanceNormalized(ic));
}

/** It returns the Twiss alpha for index 0,1,2 - x,y,z planes*/
static PyObject* BunchTwissAnalysis_getAlpha(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getAlpha", &ic)) {
    error("pyBunchTwissAnalysis.getAlpha(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getAlpha(ic));
}

/** It returns the Twiss beta for index 0,1,2 - x,y,z planes*/
static PyObject* BunchTwissAnalysis_getBeta(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_ParseTuple(arg, "i:getBeta", &ic)) {
    error("pyBunchTwissAnalysis.getBeta(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getBeta(ic));
}

/** It returns the Twiss gamma for index 0,1,2 - x,y,z planes*/
static PyObject* BunchTwissAnalysis_getGamma(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getGamma", &ic)) {
    error("pyBunchTwissAnalysis.getGamma(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getGamma(ic));
}

/** It returns Twiss dispersion function for index 0,1 - x,y planes*/
static PyObject* BunchTwissAnalysis_getDispersion(PyObject* self, PyObject* arg)
{

  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getDispersion", &ic)) {
    error("pyBunchTwissAnalysis.getDispersion(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getDispersion(ic));
}

/** It returns Twiss dispersion_prime function for index 0,1 - x,y planes*/
static PyObject* BunchTwissAnalysis_getDispersionDerivative(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getDispersionDerivative", &ic)) {
    error("pyBunchTwissAnalysis.getDispersionDerivative(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getDispersionDerivative(ic));
}

/** It returns the Twiss array (alpha,beta,gamma,emittance) for index 0,1,2 - x,y,z planes*/
static PyObject* BunchTwissAnalysis_getTwiss(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getTwiss", &ic)) {
    error("pyBunchTwissAnalysis.getTwiss(ic) - parameter is needed");
  }

  return PyTuple_Pack(4,
    PyFloat_FromDouble(cpp_BunchTwissAnalysis->getAlpha(ic)),
    PyFloat_FromDouble(cpp_BunchTwissAnalysis->getBeta(ic)),
    PyFloat_FromDouble(cpp_BunchTwissAnalysis->getGamma(ic)),
    PyFloat_FromDouble(cpp_BunchTwissAnalysis->getEmittance(ic))
  );
}

/** It returns the effective emittance for index 0,1 - x,y planes*/
static PyObject* BunchTwissAnalysis_getEffectiveEmittance(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getEffectiveEmittance", &ic)) {
    error("pyBunchTwissAnalysis.getEffectiveEmittance(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getEffectiveEmittance(ic));
}

/** It returns the effective Twiss alpha for index 0,1 - x,y planes*/
static PyObject* BunchTwissAnalysis_getEffectiveAlpha(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getEffectiveAlpha", &ic)) {
    error("pyBunchTwissAnalysis.getEffectiveAlpha(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getEffectiveAlpha(ic));
}

/** It returns the effective Twiss beta for index 0,1 - x,y planes*/
static PyObject* BunchTwissAnalysis_getEffectiveBeta(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getEffectiveBeta", &ic)) {
    error("pyBunchTwissAnalysis.getEffectiveBeta(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getEffectiveBeta(ic));
}

/** It returns the effective Twiss gamma for index 0,1 - x,y planes*/
static PyObject* BunchTwissAnalysis_getEffectiveGamma(PyObject* self, PyObject* arg)
{
  BunchTwissAnalysis* cpp_BunchTwissAnalysis =
    (BunchTwissAnalysis*)((pyORBIT_Object*)self)->cpp_obj;

  int ic;

  if (!PyArg_Parse(arg, "i:getEffectiveGamma", &ic)) {
    error("pyBunchTwissAnalysis.getEffectiveGamma(ic) - parameter is needed");
  }

  return PyFloat_FromDouble(cpp_BunchTwissAnalysis->getEffectiveGamma(ic));
}

//-----------------------------------------------------
// destructor for python BunchTwissAnalysis class (__del__ method).
//-----------------------------------------------------
static void BunchTwissAnalysis_del(pyORBIT_Object* self)
{
  // std::cerr<<"The BunchTwissAnalysis __del__ has been called!"<<std::endl;
  delete ((BunchTwissAnalysis*)self->cpp_obj);
  self->ob_base.ob_type->tp_free((PyObject*)self);
}

// defenition of the methods of the python BunchTwissAnalysis wrapper class
// they will be vailable from python level
static PyMethodDef BunchTwissAnalysisClassMethods[] =
  {{"analyzeBunch",
    BunchTwissAnalysis_analyzeBunch,
    METH_O,
    "Performs the Twiss analysis of the bunch."},
   {"computeBunchMoments",
    BunchTwissAnalysis_computeBunchMoments,
    METH_VARARGS,
    "Returns the XY moments of the beam up to a prescribed order"},
   {"getCovariance",
    BunchTwissAnalysis_getCovariance,
    METH_VARARGS,
    "Returns the centered covariance <u*v> - <u>*<v> for coordinates with "
    "indices (i, j), where u is coordinate i and v is coordinate j. "
    "Coordinate indices 0-5: 0=x, 1=x', 2=y, 3=y', 4=z, 5=dE. "
    "Returns 0.0 if i or j is out of range."},
   {"getCorrelation",
    BunchTwissAnalysis_getCorrelation,
    METH_VARARGS,
    "Returns the Pearson correlation coefficient for coordinates with "
    "indices (i, j). Coordinate indices 0-5: 0=x, 1=x', 2=y, 3=y', 4=z, 5=dE. "
    "Returns 0.0 if i or j is out of range, or if either variance is non-positive."},
   {"getBunchMoment",
    BunchTwissAnalysis_getBunchMoment,
    METH_VARARGS,
    "Returns the (i,j) xy moment of the beam"},
   {"getAverage",
    BunchTwissAnalysis_getAverage,
    METH_VARARGS,
    "Returns the average value for coordinate with index ic"},
   {"getGlobalCount",
    BunchTwissAnalysis_getGlobalCount,
    METH_NOARGS,
    "Returns the total number of analysed macroparticles"},
   {"getGlobalMacrosize",
    BunchTwissAnalysis_getGlobalMacrosize,
    METH_NOARGS,
    "Returns the total macrosize"},
   {"getEmittance",
    BunchTwissAnalysis_getEmittance,
    METH_O,
    "Returns the emittance for index 0,1,2 - x,y,z planes"},
   {"getEmittanceNormalized",
    BunchTwissAnalysis_getEmittanceNormalized,
    METH_O,
    "Returns the normalized emittance for index 0,1 - x,y planes"},
   {"getAlpha",
    BunchTwissAnalysis_getAlpha,
    METH_O,
    "Returns Twiss alpha for index 0,1,2 - x,y,z planes"},
   {"getBeta",
    BunchTwissAnalysis_getBeta,
    METH_O,
    "Returns Twiss beta for index 0,1,2 - x,y,z planes"},
   {"getGamma",
    BunchTwissAnalysis_getGamma,
    METH_O,
    "Returns Twiss gamma for index 0,1,2 - x,y,z planes"},
   {"getTwiss",
    BunchTwissAnalysis_getTwiss,
    METH_O,
    "Returns Twiss tuple (alpha,beta,gamma,emitt) for index 0,1,2 - x,y,z planes"},
   {"getDispersion",
    BunchTwissAnalysis_getDispersion,
    METH_O,
    "Returns Twiss dispersion function for index 0,1 - x,y planes"},
   {"getDispersionDerivative",
    BunchTwissAnalysis_getDispersionDerivative,
    METH_O,
    "Returns Twiss dispersion' function for index 0,1 - x,y planes"},
   {"getEffectiveEmittance",
    BunchTwissAnalysis_getEffectiveEmittance,
    METH_O,
    "Returns the effective emittance for index 0,1 - x,y planes"},
   {"getEffectiveAlpha",
    BunchTwissAnalysis_getEffectiveAlpha,
    METH_O,
    "Returns effective Twiss alpha for index 0,1 - x,y planes"},
   {"getEffectiveBeta",
    BunchTwissAnalysis_getEffectiveBeta,
    METH_O,
    "Returns effective Twiss beta for index 0,1 - x,y planes"},
   {"getEffectiveGamma",
    BunchTwissAnalysis_getEffectiveGamma,
    METH_O,
    "Returns effective Twiss gamma for index 0,1 - x,y planes"},
   {NULL, NULL, 0, NULL}};

// defenition of the memebers of the python BunchTwissAnalysis wrapper class
// they will be vailable from python level
static PyMemberDef BunchTwissAnalysisClassMembers[] = {{NULL}};

// new python BunchTwissAnalysis wrapper type definition
static PyTypeObject pyORBIT_BunchTwissAnalysis_Type = {
  PyVarObject_HEAD_INIT(NULL, 0) "BunchTwissAnalysis", /*tp_name*/
  sizeof(pyORBIT_Object),                              /*tp_basicsize*/
  0,                                                   /*tp_itemsize*/
  (destructor)BunchTwissAnalysis_del,                  /*tp_dealloc*/
  0,                                                   /*tp_print*/
  0,                                                   /*tp_getattr*/
  0,                                                   /*tp_setattr*/
  0,                                                   /*tp_compare*/
  0,                                                   /*tp_repr*/
  0,                                                   /*tp_as_number*/
  0,                                                   /*tp_as_sequence*/
  0,                                                   /*tp_as_mapping*/
  0,                                                   /*tp_hash */
  0,                                                   /*tp_call*/
  0,                                                   /*tp_str*/
  0,                                                   /*tp_getattro*/
  0,                                                   /*tp_setattro*/
  0,                                                   /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,            /*tp_flags*/
  "The BunchTwissAnalysis python wrapper",             /* tp_doc */
  0,                                                   /* tp_traverse */
  0,                                                   /* tp_clear */
  0,                                                   /* tp_richcompare */
  0,                                                   /* tp_weaklistoffset */
  0,                                                   /* tp_iter */
  0,                                                   /* tp_iternext */
  BunchTwissAnalysisClassMethods,                      /* tp_methods */
  BunchTwissAnalysisClassMembers,                      /* tp_members */
  0,                                                   /* tp_getset */
  0,                                                   /* tp_base */
  0,                                                   /* tp_dict */
  0,                                                   /* tp_descr_get */
  0,                                                   /* tp_descr_set */
  0,                                                   /* tp_dictoffset */
  (initproc)BunchTwissAnalysis_init,                   /* tp_init */
  0,                                                   /* tp_alloc */
  BunchTwissAnalysis_new,                              /* tp_new */
};

//--------------------------------------------------
// Initialization BunchTwissAnalysis of the pyBunchTwissAnalysis class
//--------------------------------------------------
void initbunchtwissanalysis(PyObject* module)
{
  if (PyType_Ready(&pyORBIT_BunchTwissAnalysis_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_BunchTwissAnalysis_Type);
  PyModule_AddObject(module, "BunchTwissAnalysis", (PyObject*)&pyORBIT_BunchTwissAnalysis_Type);
}

#ifdef __cplusplus
}
#endif

} // namespace wrap_bunch_twiss_analysis
