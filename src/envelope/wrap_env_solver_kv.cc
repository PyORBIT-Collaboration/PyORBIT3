#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "EnvSolverKV.hh"
#include "wrap_bunch.hh"
#include "wrap_envelope.hh"
#include "wrap_env_solver_kv.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping EnvSolverKV instance.
// It never will be called directly.
static PyObject *EnvSolverKV_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  pyORBIT_Object *self;
  self = (pyORBIT_Object *)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject *)self;
}

// Initialization of Python EnvSolverKV class.
// This is implementation of the __init__ method.
static int EnvSolverKV_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
  double perveance = 0.0;
  double eps_x = 1.0;
  double eps_y = 1.0;
  self->cpp_obj = new EnvSolverKV(perveance, eps_x, eps_y);
  ((EnvSolverKV *)self->cpp_obj)->setPyWrapper((PyObject *)self);
  return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *EnvSolverKV_trackBunch(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverKV = (pyORBIT_Object *)self;
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)pyEnvSolverKV->cpp_obj;
  PyObject *pyBunch;
  double length;
  if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
    ORBIT_MPI_Finalize("PyEnvSolverKV - trackBunch(Bunch* bunch, double length) - parameters are needed.");
  }
  PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("PyEnvSolverKV - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
  }
  Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
  cpp_EnvSolverKV->trackBunch(cpp_bunch, length);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *EnvSolverKV_setPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverKV = (pyORBIT_Object *)self;
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)pyEnvSolverKV->cpp_obj;
  double perveance;
  if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
    ORBIT_MPI_Finalize("PyEnvSolverKV - setPerveance(double perveance) - parameters are needed.");
  }
  cpp_EnvSolverKV->setPerveance(perveance);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceX(double emittance_x)
static PyObject *EnvSolverKV_setEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverKV = (pyORBIT_Object *)self;
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)pyEnvSolverKV->cpp_obj;
  double emittance_x;
  if (!PyArg_ParseTuple(args, "d:setEmittanceX", &emittance_x)) {
    ORBIT_MPI_Finalize("PyEnvSolverKV - setEmittanceX(double emittance_x) - parameters are needed.");
  }
  cpp_EnvSolverKV->setEmittanceX(emittance_x);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceY(double emittance_y)
static PyObject *EnvSolverKV_setEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverKV = (pyORBIT_Object *)self;
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)pyEnvSolverKV->cpp_obj;
  double emittance_y;
  if (!PyArg_ParseTuple(args, "d:setEmittanceY", &emittance_y)) {
    ORBIT_MPI_Finalize("PyEnvSolverKV - setEmittanceY(double emittance_y) - parameters are needed.");
  }
  cpp_EnvSolverKV->setEmittanceY(emittance_y);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: getPerveance()
static PyObject *EnvSolverKV_getPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverKV = (pyORBIT_Object *)self;
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)pyEnvSolverKV->cpp_obj;
  double perveance = cpp_EnvSolverKV->getPerveance();
  return Py_BuildValue("d", perveance);
}

// Method: getEmittanceX()
static PyObject *EnvSolverKV_getEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverKV = (pyORBIT_Object *)self;
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)pyEnvSolverKV->cpp_obj;
  double eps_x = cpp_EnvSolverKV->getEmittanceX();
  return Py_BuildValue("d", eps_x);
}

// Method: getEmittanceY()
static PyObject *EnvSolverKV_getEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverKV = (pyORBIT_Object *)self;
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)pyEnvSolverKV->cpp_obj;
  double eps_y = cpp_EnvSolverKV->getEmittanceY();
  return Py_BuildValue("d", eps_y);
}

// Destructor for python EnvSolverKV class (__del__ method)
static void EnvSolverKV_del(pyORBIT_Object *self) {
  EnvSolverKV *cpp_EnvSolverKV = (EnvSolverKV *)self->cpp_obj;
  delete cpp_EnvSolverKV;
  self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python EnvSolverKV wrapper class methods.
// They will be available from the Python level.
static PyMethodDef EnvSolverKVClassMethods[] = {
  {"getEmittanceX", EnvSolverKV_getEmittanceX, METH_VARARGS, "Get emittance (x)."},
  {"getEmittanceY", EnvSolverKV_getEmittanceY, METH_VARARGS, "Get emittance (y)."},
  {"getPerveance", EnvSolverKV_getPerveance, METH_VARARGS, "Get space charge perveance."},
  {"setEmittanceX", EnvSolverKV_setEmittanceX, METH_VARARGS, "Set emittance (x)."},
  {"setEmittanceY", EnvSolverKV_setEmittanceY, METH_VARARGS, "Set emittance (y)."},
  {"setPerveance", EnvSolverKV_setPerveance, METH_VARARGS, "Set space charge perveance."},
  {"trackBunch", EnvSolverKV_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
  {NULL}
};

// Definition of Python EnvSolverKV wrapper class members.
// They will be available from the Python level.
static PyMemberDef EnvSolverKVClassMembers[] = {
    {NULL}};

// New python EnvSolverKV wrapper type definition.
static PyTypeObject pyORBIT_EnvSolverKV_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "EnvSolverKV", /*tp_name*/
    sizeof(pyORBIT_Object),                       /*tp_basicsize*/
    0,                                            /*tp_itemsize*/
    (destructor)EnvSolverKV_del,                  /*tp_dealloc*/
    0,                                            /*tp_print*/
    0,                                            /*tp_getattr*/
    0,                                            /*tp_setattr*/
    0,                                            /*tp_compare*/
    0,                                            /*tp_repr*/
    0,                                            /*tp_as_number*/
    0,                                            /*tp_as_sequence*/
    0,                                            /*tp_as_mapping*/
    0,                                            /*tp_hash */
    0,                                            /*tp_call*/
    0,                                            /*tp_str*/
    0,                                            /*tp_getattro*/
    0,                                            /*tp_setattro*/
    0,                                            /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     /*tp_flags*/
    "The EnvSolverKV python wrapper",             /* tp_doc */
    0,                                            /* tp_traverse */
    0,                                            /* tp_clear */
    0,                                            /* tp_richcompare */
    0,                                            /* tp_weaklistoffset */
    0,                                            /* tp_iter */
    0,                                            /* tp_iternext */
    EnvSolverKVClassMethods,                      /* tp_methods */
    EnvSolverKVClassMembers,                      /* tp_members */
    0,                                            /* tp_getset */
    0,                                            /* tp_base */
    0,                                            /* tp_dict */
    0,                                            /* tp_descr_get */
    0,                                            /* tp_descr_set */
    0,                                            /* tp_dictoffset */
    (initproc)EnvSolverKV_init,                   /* tp_init */
    0,                                            /* tp_alloc */
    EnvSolverKV_new,                              /* tp_new */
};

// Initialization function of the pyEnvSolverKV class
void initEnvSolverKV(PyObject *module) {
  if (PyType_Ready(&pyORBIT_EnvSolverKV_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_EnvSolverKV_Type);
  PyModule_AddObject(module, "EnvSolverKV", (PyObject *)&pyORBIT_EnvSolverKV_Type);
}

#ifdef __cplusplus
}
#endif

// end of namespace wrap_envelope
} // namespace wrap_envelope
