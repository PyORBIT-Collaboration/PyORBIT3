#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "Danilov20EnvelopeSolver.hh"
#include "wrap_danilov_20_envelope_solver.hh"
#include "wrap_bunch.hh"
#include "wrap_envelope.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping Danilov20EnvelopeSolver instance.
// It never will be called directly.
static PyObject *Danilov20EnvelopeSolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  pyORBIT_Object *self;
  self = (pyORBIT_Object *)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject *)self;
}

// Initialization of Python Danilov20EnvelopeSolver class.
// This is implementation of the __init__ method.
static int Danilov20EnvelopeSolver_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
  double perveance = 0.0;
  double eps_x = 1.0;
  double eps_y = 1.0;
  self->cpp_obj = new Danilov20EnvelopeSolver(perveance, eps_x, eps_y);
  ((Danilov20EnvelopeSolver *)self->cpp_obj)->setPyWrapper((PyObject *)self);
  return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *Danilov20EnvelopeSolver_trackBunch(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov20EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)pyDanilov20EnvelopeSolver->cpp_obj;
  PyObject *pyBunch;
  double length;
  if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
    ORBIT_MPI_Finalize("PyDanilov20EnvelopeSolver - trackBunch(Bunch* bunch, double length) - parameters are needed.");
  }
  PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("PyDanilov20EnvelopeSolver - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
  }
  Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
  cpp_Danilov20EnvelopeSolver->trackBunch(cpp_bunch, length);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *Danilov20EnvelopeSolver_setPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov20EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)pyDanilov20EnvelopeSolver->cpp_obj;
  double perveance;
  if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
    ORBIT_MPI_Finalize("PyDanilov20EnvelopeSolver - setPerveance(double perveance) - parameters are needed.");
  }
  cpp_Danilov20EnvelopeSolver->setPerveance(perveance);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceX(double eps_x)
static PyObject *Danilov20EnvelopeSolver_setEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov20EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)pyDanilov20EnvelopeSolver->cpp_obj;
  double emittanceX;
  if (!PyArg_ParseTuple(args, "d:setEmittanceX", &emittanceX)) {
    ORBIT_MPI_Finalize("PyDanilov20EnvelopeSolver - setEmittanceX(double emittanceX) - parameters are needed.");
  }
  cpp_Danilov20EnvelopeSolver->setEmittanceX(emittanceX);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setEmittanceY(double eps_y)
static PyObject *Danilov20EnvelopeSolver_setEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov20EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)pyDanilov20EnvelopeSolver->cpp_obj;
  double emittanceY;
  if (!PyArg_ParseTuple(args, "d:setEmittanceY", &emittanceY)) {
    ORBIT_MPI_Finalize("PyDanilov20EnvelopeSolver - setEmittanceY(double emittanceY) - parameters are needed.");
  }
  cpp_Danilov20EnvelopeSolver->setEmittanceY(emittanceY);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: getPerveance()
static PyObject *Danilov20EnvelopeSolver_getPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov20EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)pyDanilov20EnvelopeSolver->cpp_obj;
  double perveance = cpp_Danilov20EnvelopeSolver->getPerveance();
  return Py_BuildValue("d", perveance);
}

// Method: getEmittanceX()
static PyObject *Danilov20EnvelopeSolver_getEmittanceX(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov20EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)pyDanilov20EnvelopeSolver->cpp_obj;
  double emittanceX = cpp_Danilov20EnvelopeSolver->getEmittanceX();
  return Py_BuildValue("d", emittanceX);
}

// Method: getEmittanceY()
static PyObject *Danilov20EnvelopeSolver_getEmittanceY(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov20EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)pyDanilov20EnvelopeSolver->cpp_obj;
  double emittanceY = cpp_Danilov20EnvelopeSolver->getEmittanceY();
  return Py_BuildValue("d", emittanceY);
}

// Destructor for python Danilov20EnvelopeSolver class (__del__ method)
static void Danilov20EnvelopeSolver_del(pyORBIT_Object *self) {
  Danilov20EnvelopeSolver *cpp_Danilov20EnvelopeSolver = (Danilov20EnvelopeSolver *)self->cpp_obj;
  delete cpp_Danilov20EnvelopeSolver;
  self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python Danilov20EnvelopeSolver wrapper class methods.
// They will be available from the Python level.
static PyMethodDef Danilov20EnvelopeSolverClassMethods[] = {
  {"getEmittanceX", Danilov20EnvelopeSolver_getEmittanceX, METH_VARARGS, "Get emittance (x)."},
  {"getEmittanceY", Danilov20EnvelopeSolver_getEmittanceY, METH_VARARGS, "Get emittance (y)."},
  {"getPerveance", Danilov20EnvelopeSolver_getPerveance, METH_VARARGS, "Get space charge perveance."},
  {"setEmittanceX", Danilov20EnvelopeSolver_setEmittanceX, METH_VARARGS, "Set emittance (x)."},
  {"setEmittanceY", Danilov20EnvelopeSolver_setEmittanceY, METH_VARARGS, "Set emittance (y)."},
  {"setPerveance", Danilov20EnvelopeSolver_setPerveance, METH_VARARGS, "Set space charge perveance."},
  {"trackBunch", Danilov20EnvelopeSolver_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
  {NULL}
};

// Definition of Python Danilov20EnvelopeSolver wrapper class members.
// They will be available from the Python level.
static PyMemberDef Danilov20EnvelopeSolverClassMembers[] = {
    {NULL}};

// New python Danilov20EnvelopeSolver wrapper type definition.
static PyTypeObject pyORBIT_Danilov20EnvelopeSolver_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "Danilov20EnvelopeSolver", /*tp_name*/
    sizeof(pyORBIT_Object),                                   /*tp_basicsize*/
    0,                                                        /*tp_itemsize*/
    (destructor)Danilov20EnvelopeSolver_del,                  /*tp_dealloc*/
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
    "The Danilov20EnvelopeSolver python wrapper",             /* tp_doc */
    0,                                                        /* tp_traverse */
    0,                                                        /* tp_clear */
    0,                                                        /* tp_richcompare */
    0,                                                        /* tp_weaklistoffset */
    0,                                                        /* tp_iter */
    0,                                                        /* tp_iternext */
    Danilov20EnvelopeSolverClassMethods,                      /* tp_methods */
    Danilov20EnvelopeSolverClassMembers,                      /* tp_members */
    0,                                                        /* tp_getset */
    0,                                                        /* tp_base */
    0,                                                        /* tp_dict */
    0,                                                        /* tp_descr_get */
    0,                                                        /* tp_descr_set */
    0,                                                        /* tp_dictoffset */
    (initproc)Danilov20EnvelopeSolver_init,                   /* tp_init */
    0,                                                        /* tp_alloc */
    Danilov20EnvelopeSolver_new,                              /* tp_new */
};

// Initialization function of the pyDanilov20EnvelopeSolver class
void initDanilov20EnvelopeSolver(PyObject *module) {
  if (PyType_Ready(&pyORBIT_Danilov20EnvelopeSolver_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_Danilov20EnvelopeSolver_Type);
  PyModule_AddObject(module, "Danilov20EnvelopeSolver", (PyObject *)&pyORBIT_Danilov20EnvelopeSolver_Type);
}

#ifdef __cplusplus
}
#endif

// end of namespace wrap_danilov_envelope
} // namespace wrap_danilov_envelope
