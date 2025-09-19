#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "Danilov22EnvelopeSolver.hh"
#include "wrap_danilov_22_envelope_solver.hh"
#include "wrap_bunch.hh"
#include "wrap_envelope.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping Danilov22EnvelopeSolver instance.
// It never will be called directly.
static PyObject *Danilov22EnvelopeSolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  pyORBIT_Object *self;
  self = (pyORBIT_Object *)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject *)self;
}

// Initialization of Python Danilov22EnvelopeSolver class.
// This is implementation of the __init__ method.
static int Danilov22EnvelopeSolver_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
  double perveance = 0.0;
  self->cpp_obj = new Danilov22EnvelopeSolver(perveance);
  ((Danilov22EnvelopeSolver *)self->cpp_obj)->setPyWrapper((PyObject *)self);
  return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *Danilov22EnvelopeSolver_trackBunch(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov22EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov22EnvelopeSolver *cpp_Danilov22EnvelopeSolver = (Danilov22EnvelopeSolver *)pyDanilov22EnvelopeSolver->cpp_obj;
  PyObject *pyBunch;
  double length;
  if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
    ORBIT_MPI_Finalize("PyDanilov22EnvelopeSolver - trackBunch(Bunch* bunch, double length) - parameters are needed.");
  }
  PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("PyDanilov22EnvelopeSolver - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
  }
  Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
  cpp_Danilov22EnvelopeSolver->trackBunch(cpp_bunch, length);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *Danilov22EnvelopeSolver_setPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov22EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov22EnvelopeSolver *cpp_Danilov22EnvelopeSolver = (Danilov22EnvelopeSolver *)pyDanilov22EnvelopeSolver->cpp_obj;
  double perveance;
  if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
    ORBIT_MPI_Finalize("PyDanilov22EnvelopeSolver - setPerveance(double perveance) - parameters are needed.");
  }
  cpp_Danilov22EnvelopeSolver->setPerveance(perveance);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: getPerveance()
static PyObject *Danilov22EnvelopeSolver_getPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilov22EnvelopeSolver = (pyORBIT_Object *)self;
  Danilov22EnvelopeSolver *cpp_Danilov22EnvelopeSolver = (Danilov22EnvelopeSolver *)pyDanilov22EnvelopeSolver->cpp_obj;
  double perveance = cpp_Danilov22EnvelopeSolver->getPerveance();
  return Py_BuildValue("d", perveance);
}

// Destructor for python Danilov22EnvelopeSolver class (__del__ method)
static void Danilov22EnvelopeSolver_del(pyORBIT_Object *self) {
  Danilov22EnvelopeSolver *cpp_Danilov22EnvelopeSolver = (Danilov22EnvelopeSolver *)self->cpp_obj;
  delete cpp_Danilov22EnvelopeSolver;
  self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python Danilov22EnvelopeSolver wrapper class methods.
// They will be available from the Python level.
static PyMethodDef Danilov22EnvelopeSolverClassMethods[] = {
  {"trackBunch", Danilov22EnvelopeSolver_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
  {"setPerveance", Danilov22EnvelopeSolver_setPerveance, METH_VARARGS, "Set space charge perveance."},
  {"getPerveance", Danilov22EnvelopeSolver_getPerveance, METH_VARARGS, "Get space charge perveance."},
  {NULL}
};

// Definition of Python Danilov22EnvelopeSolver wrapper class members.
// They will be available from the Python level.
static PyMemberDef Danilov22EnvelopeSolverClassMembers[] = {
    {NULL}};

// New python Danilov22EnvelopeSolver wrapper type definition.
static PyTypeObject pyORBIT_Danilov22EnvelopeSolver_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "Danilov22EnvelopeSolver", /*tp_name*/
    sizeof(pyORBIT_Object),                                   /*tp_basicsize*/
    0,                                                        /*tp_itemsize*/
    (destructor)Danilov22EnvelopeSolver_del,                  /*tp_dealloc*/
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
    "The Danilov22EnvelopeSolver python wrapper",             /* tp_doc */
    0,                                                        /* tp_traverse */
    0,                                                        /* tp_clear */
    0,                                                        /* tp_richcompare */
    0,                                                        /* tp_weaklistoffset */
    0,                                                        /* tp_iter */
    0,                                                        /* tp_iternext */
    Danilov22EnvelopeSolverClassMethods,                      /* tp_methods */
    Danilov22EnvelopeSolverClassMembers,                      /* tp_members */
    0,                                                        /* tp_getset */
    0,                                                        /* tp_base */
    0,                                                        /* tp_dict */
    0,                                                        /* tp_descr_get */
    0,                                                        /* tp_descr_set */
    0,                                                        /* tp_dictoffset */
    (initproc)Danilov22EnvelopeSolver_init,                   /* tp_init */
    0,                                                        /* tp_alloc */
    Danilov22EnvelopeSolver_new,                              /* tp_new */
};

// Initialization function of the pyDanilov22EnvelopeSolver class
void initDanilov22EnvelopeSolver(PyObject *module) {
  if (PyType_Ready(&pyORBIT_Danilov22EnvelopeSolver_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_Danilov22EnvelopeSolver_Type);
  PyModule_AddObject(module, "Danilov22EnvelopeSolver", (PyObject *)&pyORBIT_Danilov22EnvelopeSolver_Type);
}

#ifdef __cplusplus
}
#endif

// end of namespace wrap_danilov_envelope
} // namespace wrap_danilov_envelope