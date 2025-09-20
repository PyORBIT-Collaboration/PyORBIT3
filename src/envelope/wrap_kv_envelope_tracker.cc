#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "KVEnvelopeTracker.hh"
#include "wrap_bunch.hh"
#include "wrap_envelope.hh"
#include "wrap_kv_envelope_tracker.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping KVEnvelopeTracker instance.
// It never will be called directly.
static PyObject *KVEnvelopeTracker_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  pyORBIT_Object *self;
  self = (pyORBIT_Object *)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject *)self;
}

// Initialization of Python KVEnvelopeTracker class.
// This is implementation of the __init__ method.
static int KVEnvelopeTracker_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
  double perveance = 0.0;
  double eps_x = 1.0;
  double eps_y = 1.0;
  self->cpp_obj = new KVEnvelopeTracker(perveance, eps_x, eps_y);
  ((KVEnvelopeTracker *)self->cpp_obj)->setPyWrapper((PyObject *)self);
  return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *KVEnvelopeTracker_trackBunch(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyKVEnvelopeTracker = (pyORBIT_Object *)self;
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)pyKVEnvelopeTracker->cpp_obj;
  PyObject *pyBunch;
  double length;
  if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
    ORBIT_MPI_Finalize("PyKVEnvelopeTracker - trackBunch(Bunch* bunch, double length) - parameters are needed.");
  }
  PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("PyKVEnvelopeTracker - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
  }
  Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
  cpp_KVEnvelopeTracker->trackBunch(cpp_bunch, length);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *KVEnvelopeTracker_setPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyKVEnvelopeTracker = (pyORBIT_Object *)self;
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)pyKVEnvelopeTracker->cpp_obj;
  double perveance;
  if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
    ORBIT_MPI_Finalize("PyKVEnvelopeTracker - setPerveance(double perveance) - parameters are needed.");
  }
  cpp_KVEnvelopeTracker->setPerveance(perveance);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceX(double emittance_x)
static PyObject *KVEnvelopeTracker_setEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyKVEnvelopeTracker = (pyORBIT_Object *)self;
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)pyKVEnvelopeTracker->cpp_obj;
  double emittance_x;
  if (!PyArg_ParseTuple(args, "d:setEmittanceX", &emittance_x)) {
    ORBIT_MPI_Finalize("PyKVEnvelopeTracker - setEmittanceX(double emittance_x) - parameters are needed.");
  }
  cpp_KVEnvelopeTracker->setEmittanceX(emittance_x);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceY(double emittance_y)
static PyObject *KVEnvelopeTracker_setEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyKVEnvelopeTracker = (pyORBIT_Object *)self;
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)pyKVEnvelopeTracker->cpp_obj;
  double emittance_y;
  if (!PyArg_ParseTuple(args, "d:setEmittanceY", &emittance_y)) {
    ORBIT_MPI_Finalize("PyKVEnvelopeTracker - setEmittanceY(double emittance_y) - parameters are needed.");
  }
  cpp_KVEnvelopeTracker->setEmittanceY(emittance_y);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: getPerveance()
static PyObject *KVEnvelopeTracker_getPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyKVEnvelopeTracker = (pyORBIT_Object *)self;
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)pyKVEnvelopeTracker->cpp_obj;
  double perveance = cpp_KVEnvelopeTracker->getPerveance();
  return Py_BuildValue("d", perveance);
}

// Method: getEmittanceX()
static PyObject *KVEnvelopeTracker_getEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyKVEnvelopeTracker = (pyORBIT_Object *)self;
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)pyKVEnvelopeTracker->cpp_obj;
  double eps_x = cpp_KVEnvelopeTracker->getEmittanceX();
  return Py_BuildValue("d", eps_x);
}

// Method: getEmittanceY()
static PyObject *KVEnvelopeTracker_getEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyKVEnvelopeTracker = (pyORBIT_Object *)self;
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)pyKVEnvelopeTracker->cpp_obj;
  double eps_y = cpp_KVEnvelopeTracker->getEmittanceY();
  return Py_BuildValue("d", eps_y);
}

// Destructor for python KVEnvelopeTracker class (__del__ method)
static void KVEnvelopeTracker_del(pyORBIT_Object *self) {
  KVEnvelopeTracker *cpp_KVEnvelopeTracker = (KVEnvelopeTracker *)self->cpp_obj;
  delete cpp_KVEnvelopeTracker;
  self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python KVEnvelopeTracker wrapper class methods.
// They will be available from the Python level.
static PyMethodDef KVEnvelopeTrackerClassMethods[] = {
  {"getEmittanceX", KVEnvelopeTracker_getEmittanceX, METH_VARARGS, "Get emittance (x)."},
  {"getEmittanceY", KVEnvelopeTracker_getEmittanceY, METH_VARARGS, "Get emittance (y)."},
  {"getPerveance", KVEnvelopeTracker_getPerveance, METH_VARARGS, "Get space charge perveance."},
  {"setEmittanceX", KVEnvelopeTracker_setEmittanceX, METH_VARARGS, "Set emittance (x)."},
  {"setEmittanceY", KVEnvelopeTracker_setEmittanceY, METH_VARARGS, "Set emittance (y)."},
  {"setPerveance", KVEnvelopeTracker_setPerveance, METH_VARARGS, "Set space charge perveance."},
  {"trackBunch", KVEnvelopeTracker_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
  {NULL}
};

// Definition of Python KVEnvelopeTracker wrapper class members.
// They will be available from the Python level.
static PyMemberDef KVEnvelopeTrackerClassMembers[] = {
    {NULL}};

// New python KVEnvelopeTracker wrapper type definition.
static PyTypeObject pyORBIT_KVEnvelopeTracker_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "KVEnvelopeTracker", /*tp_name*/
    sizeof(pyORBIT_Object),                       /*tp_basicsize*/
    0,                                            /*tp_itemsize*/
    (destructor)KVEnvelopeTracker_del,                  /*tp_dealloc*/
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
    "The KVEnvelopeTracker python wrapper",             /* tp_doc */
    0,                                            /* tp_traverse */
    0,                                            /* tp_clear */
    0,                                            /* tp_richcompare */
    0,                                            /* tp_weaklistoffset */
    0,                                            /* tp_iter */
    0,                                            /* tp_iternext */
    KVEnvelopeTrackerClassMethods,                      /* tp_methods */
    KVEnvelopeTrackerClassMembers,                      /* tp_members */
    0,                                            /* tp_getset */
    0,                                            /* tp_base */
    0,                                            /* tp_dict */
    0,                                            /* tp_descr_get */
    0,                                            /* tp_descr_set */
    0,                                            /* tp_dictoffset */
    (initproc)KVEnvelopeTracker_init,                   /* tp_init */
    0,                                            /* tp_alloc */
    KVEnvelopeTracker_new,                              /* tp_new */
};

// Initialization function of the pyKVEnvelopeTracker class
void initKVEnvelopeTracker(PyObject *module) {
  if (PyType_Ready(&pyORBIT_KVEnvelopeTracker_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_KVEnvelopeTracker_Type);
  PyModule_AddObject(module, "KVEnvelopeTracker", (PyObject *)&pyORBIT_KVEnvelopeTracker_Type);
}

#ifdef __cplusplus
}
#endif

// end of namespace wrap_envelope
} // namespace wrap_envelope
