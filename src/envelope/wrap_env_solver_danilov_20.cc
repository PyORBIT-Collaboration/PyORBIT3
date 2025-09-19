#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "EnvSolverDanilov20.hh"
#include "wrap_bunch.hh"
#include "wrap_envelope.hh"
#include "wrap_env_solver_danilov_20.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping EnvSolverDanilov20 instance.
// It never will be called directly.
static PyObject *EnvSolverDanilov20_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  pyORBIT_Object *self;
  self = (pyORBIT_Object *)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject *)self;
}

// Initialization of Python EnvSolverDanilov20 class.
// This is implementation of the __init__ method.
static int EnvSolverDanilov20_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
  double perveance = 0.0;
  double eps_x = 1.0;
  double eps_y = 1.0;
  self->cpp_obj = new EnvSolverDanilov20(perveance, eps_x, eps_y);
  ((EnvSolverDanilov20 *)self->cpp_obj)->setPyWrapper((PyObject *)self);
  return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *EnvSolverDanilov20_trackBunch(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverDanilov20 = (pyORBIT_Object *)self;
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)pyEnvSolverDanilov20->cpp_obj;
  PyObject *pyBunch;
  double length;
  if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
    ORBIT_MPI_Finalize("PyEnvSolverDanilov20 - trackBunch(Bunch* bunch, double length) - parameters are needed.");
  }
  PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("PyEnvSolverDanilov20 - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
  }
  Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
  cpp_EnvSolverDanilov20->trackBunch(cpp_bunch, length);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *EnvSolverDanilov20_setPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverDanilov20 = (pyORBIT_Object *)self;
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)pyEnvSolverDanilov20->cpp_obj;
  double perveance;
  if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
    ORBIT_MPI_Finalize("PyEnvSolverDanilov20 - setPerveance(double perveance) - parameters are needed.");
  }
  cpp_EnvSolverDanilov20->setPerveance(perveance);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceX(double eps_x)
static PyObject *EnvSolverDanilov20_setEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverDanilov20 = (pyORBIT_Object *)self;
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)pyEnvSolverDanilov20->cpp_obj;
  double eps_x;
  if (!PyArg_ParseTuple(args, "d:setEmittanceX", &eps_x)) {
    ORBIT_MPI_Finalize("PyEnvSolverDanilov20 - setEmittanceX(double eps_x) - parameters are needed.");
  }
  cpp_EnvSolverDanilov20->setEmittanceX(eps_x);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceY(double eps_y)
static PyObject *EnvSolverDanilov20_setEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverDanilov20 = (pyORBIT_Object *)self;
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)pyEnvSolverDanilov20->cpp_obj;
  double eps_y;
  if (!PyArg_ParseTuple(args, "d:setEmittanceY", &eps_y)) {
    ORBIT_MPI_Finalize("PyEnvSolverDanilov20 - setEmittanceY(double eps_y) - parameters are needed.");
  }
  cpp_EnvSolverDanilov20->setEmittanceY(eps_y);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: getPerveance()
static PyObject *EnvSolverDanilov20_getPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverDanilov20 = (pyORBIT_Object *)self;
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)pyEnvSolverDanilov20->cpp_obj;
  double perveance = cpp_EnvSolverDanilov20->getPerveance();
  return Py_BuildValue("d", perveance);
}

// Method: getEmittanceX()
static PyObject *EnvSolverDanilov20_getEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverDanilov20 = (pyORBIT_Object *)self;
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)pyEnvSolverDanilov20->cpp_obj;
  double eps_x = cpp_EnvSolverDanilov20->getEmittanceX();
  return Py_BuildValue("d", eps_x);
}

// Method: getEmittanceY()
static PyObject *EnvSolverDanilov20_getEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyEnvSolverDanilov20 = (pyORBIT_Object *)self;
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)pyEnvSolverDanilov20->cpp_obj;
  double eps_y = cpp_EnvSolverDanilov20->getEmittanceY();
  return Py_BuildValue("d", eps_y);
}

// Destructor for python EnvSolverDanilov20 class (__del__ method)
static void EnvSolverDanilov20_del(pyORBIT_Object *self) {
  EnvSolverDanilov20 *cpp_EnvSolverDanilov20 = (EnvSolverDanilov20 *)self->cpp_obj;
  delete cpp_EnvSolverDanilov20;
  self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python EnvSolverDanilov20 wrapper class methods.
// They will be available from the Python level.
static PyMethodDef EnvSolverDanilov20ClassMethods[] = {
  {"getEmittanceX", EnvSolverDanilov20_getEmittanceX, METH_VARARGS, "Get emittance (x)."},
  {"getEmittanceY", EnvSolverDanilov20_getEmittanceY, METH_VARARGS, "Get emittance (y)."},
  {"getPerveance", EnvSolverDanilov20_getPerveance, METH_VARARGS, "Get space charge perveance."},
  {"setEmittanceX", EnvSolverDanilov20_setEmittanceX, METH_VARARGS, "Set emittance (x)."},
  {"setEmittanceY", EnvSolverDanilov20_setEmittanceY, METH_VARARGS, "Set emittance (y)."},
  {"setPerveance", EnvSolverDanilov20_setPerveance, METH_VARARGS, "Set space charge perveance."},
  {"trackBunch", EnvSolverDanilov20_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
  {NULL}
};

// Definition of Python EnvSolverDanilov20 wrapper class members.
// They will be available from the Python level.
static PyMemberDef EnvSolverDanilov20ClassMembers[] = {
    {NULL}};

// New python EnvSolverDanilov20 wrapper type definition.
static PyTypeObject pyORBIT_EnvSolverDanilov20_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "EnvSolverDanilov20", /*tp_name*/
    sizeof(pyORBIT_Object),                                   /*tp_basicsize*/
    0,                                                        /*tp_itemsize*/
    (destructor)EnvSolverDanilov20_del,                  /*tp_dealloc*/
    0,                                                        /*tp_print*/
    0,                                                        /*tp_getattr*/
    0,                                                        /*tp_setattr*/
    0,                                                        /*tp_compare*/
    0,                                                        /*tp_repr*/
    0,                                                        /*tp_as_number*/
    0,                                                        /*tp_as_sequence*/
    0,                                                        /*tp_as_mapping*/
    0,                                                        /*tp_hash */
    0,                                                        /*tp_call*/
    0,                                                        /*tp_str*/
    0,                                                        /*tp_getattro*/
    0,                                                        /*tp_setattro*/
    0,                                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,                 /*tp_flags*/
    "The EnvSolverDanilov20 python wrapper",             /* tp_doc */
    0,                                                        /* tp_traverse */
    0,                                                        /* tp_clear */
    0,                                                        /* tp_richcompare */
    0,                                                        /* tp_weaklistoffset */
    0,                                                        /* tp_iter */
    0,                                                        /* tp_iternext */
    EnvSolverDanilov20ClassMethods,                      /* tp_methods */
    EnvSolverDanilov20ClassMembers,                      /* tp_members */
    0,                                                        /* tp_getset */
    0,                                                        /* tp_base */
    0,                                                        /* tp_dict */
    0,                                                        /* tp_descr_get */
    0,                                                        /* tp_descr_set */
    0,                                                        /* tp_dictoffset */
    (initproc)EnvSolverDanilov20_init,                   /* tp_init */
    0,                                                        /* tp_alloc */
    EnvSolverDanilov20_new,                              /* tp_new */
};

// Initialization function of the pyEnvSolverDanilov20 class
void initEnvSolverDanilov20(PyObject *module) {
  if (PyType_Ready(&pyORBIT_EnvSolverDanilov20_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_EnvSolverDanilov20_Type);
  PyModule_AddObject(module, "EnvSolverDanilov20", (PyObject *)&pyORBIT_EnvSolverDanilov20_Type);
}

#ifdef __cplusplus
}
#endif

// end of namespace wrap_envelope
} // namespace wrap_envelope