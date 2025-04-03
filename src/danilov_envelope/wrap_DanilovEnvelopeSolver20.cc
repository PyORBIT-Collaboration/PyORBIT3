#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "DanilovEnvelopeSolver20.hh"
#include "wrap_DanilovEnvelopeSolver20.hh"
#include "wrap_bunch.hh"
#include "wrap_danilov_envelope.hh"

namespace wrap_danilov_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping DanilovEnvelopeSolver20 instance.
// It never will be called directly.
static PyObject *DanilovEnvelopeSolver20_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  pyORBIT_Object *self;
  self = (pyORBIT_Object *)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject *)self;
}

// Initialization of Python DanilovEnvelopeSolver20 class.
// This is implementation of the __init__ method.
static int DanilovEnvelopeSolver20_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
  double perveance = 0.0;
  double eps_x = 1.0;
  double eps_y = 1.0;
  self->cpp_obj = new DanilovEnvelopeSolver20(perveance, eps_x, eps_y);
  ((DanilovEnvelopeSolver20 *)self->cpp_obj)->setPyWrapper((PyObject *)self);
  return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *DanilovEnvelopeSolver20_trackBunch(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver20 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)pyDanilovEnvelopeSolver20->cpp_obj;
  PyObject *pyBunch;
  double length;
  if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver20 - trackBunch(Bunch* bunch, double length) - parameters are needed.");
  }
  PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver20 - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
  }
  Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
  cpp_DanilovEnvelopeSolver20->trackBunch(cpp_bunch, length);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *DanilovEnvelopeSolver20_setPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver20 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)pyDanilovEnvelopeSolver20->cpp_obj;
  double perveance;
  if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver20 - setPerveance(double perveance) - parameters are needed.");
  }
  cpp_DanilovEnvelopeSolver20->setPerveance(perveance);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceX(double eps_x)
static PyObject *DanilovEnvelopeSolver20_setEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver20 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)pyDanilovEnvelopeSolver20->cpp_obj;
  double eps_x;
  if (!PyArg_ParseTuple(args, "d:setEmittanceX", &eps_x)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver20 - setEmittanceX(double eps_x) - parameters are needed.");
  }
  cpp_DanilovEnvelopeSolver20->setEmittanceX(eps_x);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceY(double eps_y)
static PyObject *DanilovEnvelopeSolver20_setEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver20 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)pyDanilovEnvelopeSolver20->cpp_obj;
  double eps_y;
  if (!PyArg_ParseTuple(args, "d:setEmittanceY", &eps_y)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver20 - setEmittanceY(double eps_y) - parameters are needed.");
  }
  cpp_DanilovEnvelopeSolver20->setEmittanceY(eps_y);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: getPerveance()
static PyObject *DanilovEnvelopeSolver20_getPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver20 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)pyDanilovEnvelopeSolver20->cpp_obj;
  double perveance = cpp_DanilovEnvelopeSolver20->getPerveance();
  return Py_BuildValue("d", perveance);
}

// Method: getEmittanceX()
static PyObject *DanilovEnvelopeSolver20_getEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver20 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)pyDanilovEnvelopeSolver20->cpp_obj;
  double eps_x = cpp_DanilovEnvelopeSolver20->getEmittanceX();
  return Py_BuildValue("d", eps_x);
}

// Method: getEmittanceY()
static PyObject *DanilovEnvelopeSolver20_getEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver20 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)pyDanilovEnvelopeSolver20->cpp_obj;
  double eps_y = cpp_DanilovEnvelopeSolver20->getEmittanceY();
  return Py_BuildValue("d", eps_y);
}

// Destructor for python DanilovEnvelopeSolver20 class (__del__ method)
static void DanilovEnvelopeSolver20_del(pyORBIT_Object *self) {
  DanilovEnvelopeSolver20 *cpp_DanilovEnvelopeSolver20 = (DanilovEnvelopeSolver20 *)self->cpp_obj;
  delete cpp_DanilovEnvelopeSolver20;
  self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python DanilovEnvelopeSolver20 wrapper class methods.
// They will be available from the Python level.
static PyMethodDef DanilovEnvelopeSolver20ClassMethods[] = {
  {"getEmittanceX", DanilovEnvelopeSolver20_getEmittanceX, METH_VARARGS, "Get emittance (x)."},
  {"getEmittanceY", DanilovEnvelopeSolver20_getEmittanceY, METH_VARARGS, "Get emittance (y)."},
  {"getPerveance", DanilovEnvelopeSolver20_getPerveance, METH_VARARGS, "Get space charge perveance."},
  {"setEmittanceX", DanilovEnvelopeSolver20_setEmittanceX, METH_VARARGS, "Set emittance (x)."},
  {"setEmittanceY", DanilovEnvelopeSolver20_setEmittanceY, METH_VARARGS, "Set emittance (y)."},
  {"setPerveance", DanilovEnvelopeSolver20_setPerveance, METH_VARARGS, "Set space charge perveance."},
  {"trackBunch", DanilovEnvelopeSolver20_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
  {NULL}
};

// Definition of Python DanilovEnvelopeSolver20 wrapper class members.
// They will be available from the Python level.
static PyMemberDef DanilovEnvelopeSolver20ClassMembers[] = {
    {NULL}};

// New python DanilovEnvelopeSolver20 wrapper type definition.
static PyTypeObject pyORBIT_DanilovEnvelopeSolver20_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "DanilovEnvelopeSolver20", /*tp_name*/
    sizeof(pyORBIT_Object),                                   /*tp_basicsize*/
    0,                                                        /*tp_itemsize*/
    (destructor)DanilovEnvelopeSolver20_del,                  /*tp_dealloc*/
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
    "The DanilovEnvelopeSolver20 python wrapper",             /* tp_doc */
    0,                                                        /* tp_traverse */
    0,                                                        /* tp_clear */
    0,                                                        /* tp_richcompare */
    0,                                                        /* tp_weaklistoffset */
    0,                                                        /* tp_iter */
    0,                                                        /* tp_iternext */
    DanilovEnvelopeSolver20ClassMethods,                      /* tp_methods */
    DanilovEnvelopeSolver20ClassMembers,                      /* tp_members */
    0,                                                        /* tp_getset */
    0,                                                        /* tp_base */
    0,                                                        /* tp_dict */
    0,                                                        /* tp_descr_get */
    0,                                                        /* tp_descr_set */
    0,                                                        /* tp_dictoffset */
    (initproc)DanilovEnvelopeSolver20_init,                   /* tp_init */
    0,                                                        /* tp_alloc */
    DanilovEnvelopeSolver20_new,                              /* tp_new */
};

// Initialization function of the pyDanilovEnvelopeSolver20 class
void initDanilovEnvelopeSolver20(PyObject *module) {
  if (PyType_Ready(&pyORBIT_DanilovEnvelopeSolver20_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_DanilovEnvelopeSolver20_Type);
  PyModule_AddObject(module, "DanilovEnvelopeSolver20", (PyObject *)&pyORBIT_DanilovEnvelopeSolver20_Type);
}

#ifdef __cplusplus
}
#endif

// end of namespace wrap_danilov_envelope
} // namespace wrap_danilov_envelope