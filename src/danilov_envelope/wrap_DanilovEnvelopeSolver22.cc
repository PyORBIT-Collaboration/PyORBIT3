#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "DanilovEnvelopeSolver22.hh"
#include "wrap_DanilovEnvelopeSolver22.hh"
#include "wrap_bunch.hh"
#include "wrap_danilov_envelope.hh"

namespace wrap_danilov_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping DanilovEnvelopeSolver22 instance.
// It never will be called directly.
static PyObject *DanilovEnvelopeSolver22_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  pyORBIT_Object *self;
  self = (pyORBIT_Object *)type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  return (PyObject *)self;
}

// Initialization of Python DanilovEnvelopeSolver22 class.
// This is implementation of the __init__ method.
static int DanilovEnvelopeSolver22_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
  double perveance = 0.0;
  self->cpp_obj = new DanilovEnvelopeSolver22(perveance);
  ((DanilovEnvelopeSolver22 *)self->cpp_obj)->setPyWrapper((PyObject *)self);
  return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *DanilovEnvelopeSolver22_trackBunch(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver22 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver22 *cpp_DanilovEnvelopeSolver22 = (DanilovEnvelopeSolver22 *)pyDanilovEnvelopeSolver22->cpp_obj;
  PyObject *pyBunch;
  double length;
  if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver22 - trackBunch(Bunch* bunch, double length) - parameters are needed.");
  }
  PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
  if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver22 - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
  }
  Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
  cpp_DanilovEnvelopeSolver22->trackBunch(cpp_bunch, length);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *DanilovEnvelopeSolver22_setPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver22 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver22 *cpp_DanilovEnvelopeSolver22 = (DanilovEnvelopeSolver22 *)pyDanilovEnvelopeSolver22->cpp_obj;
  double perveance;
  if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
    ORBIT_MPI_Finalize("PyDanilovEnvelopeSolver22 - setPerveance(double perveance) - parameters are needed.");
  }
  cpp_DanilovEnvelopeSolver22->setPerveance(perveance);
  Py_INCREF(Py_None);
  return Py_None;
}

// Method: getPerveance()
static PyObject *DanilovEnvelopeSolver22_getPerveance(PyObject *self, PyObject *args) {
  pyORBIT_Object *pyDanilovEnvelopeSolver22 = (pyORBIT_Object *)self;
  DanilovEnvelopeSolver22 *cpp_DanilovEnvelopeSolver22 = (DanilovEnvelopeSolver22 *)pyDanilovEnvelopeSolver22->cpp_obj;
  double perveance = cpp_DanilovEnvelopeSolver22->getPerveance();
  return Py_BuildValue("d", perveance);
}

// Destructor for python DanilovEnvelopeSolver22 class (__del__ method)
static void DanilovEnvelopeSolver22_del(pyORBIT_Object *self) {
  DanilovEnvelopeSolver22 *cpp_DanilovEnvelopeSolver22 = (DanilovEnvelopeSolver22 *)self->cpp_obj;
  delete cpp_DanilovEnvelopeSolver22;
  self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python DanilovEnvelopeSolver22 wrapper class methods.
// They will be available from the Python level.
static PyMethodDef DanilovEnvelopeSolver22ClassMethods[] = {
  {"trackBunch", DanilovEnvelopeSolver22_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
  {"setPerveance", DanilovEnvelopeSolver22_setPerveance, METH_VARARGS, "Set space charge perveance."},
  {"getPerveance", DanilovEnvelopeSolver22_getPerveance, METH_VARARGS, "Get space charge perveance."},
  {NULL}
};

// Definition of Python DanilovEnvelopeSolver22 wrapper class members.
// They will be available from the Python level.
static PyMemberDef DanilovEnvelopeSolver22ClassMembers[] = {
    {NULL}};

// New python DanilovEnvelopeSolver22 wrapper type definition.
static PyTypeObject pyORBIT_DanilovEnvelopeSolver22_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "DanilovEnvelopeSolver22", /*tp_name*/
    sizeof(pyORBIT_Object),                                   /*tp_basicsize*/
    0,                                                        /*tp_itemsize*/
    (destructor)DanilovEnvelopeSolver22_del,                  /*tp_dealloc*/
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
    "The DanilovEnvelopeSolver22 python wrapper",             /* tp_doc */
    0,                                                        /* tp_traverse */
    0,                                                        /* tp_clear */
    0,                                                        /* tp_richcompare */
    0,                                                        /* tp_weaklistoffset */
    0,                                                        /* tp_iter */
    0,                                                        /* tp_iternext */
    DanilovEnvelopeSolver22ClassMethods,                      /* tp_methods */
    DanilovEnvelopeSolver22ClassMembers,                      /* tp_members */
    0,                                                        /* tp_getset */
    0,                                                        /* tp_base */
    0,                                                        /* tp_dict */
    0,                                                        /* tp_descr_get */
    0,                                                        /* tp_descr_set */
    0,                                                        /* tp_dictoffset */
    (initproc)DanilovEnvelopeSolver22_init,                   /* tp_init */
    0,                                                        /* tp_alloc */
    DanilovEnvelopeSolver22_new,                              /* tp_new */
};

// Initialization function of the pyDanilovEnvelopeSolver22 class
void initDanilovEnvelopeSolver22(PyObject *module) {
  if (PyType_Ready(&pyORBIT_DanilovEnvelopeSolver22_Type) < 0)
    return;
  Py_INCREF(&pyORBIT_DanilovEnvelopeSolver22_Type);
  PyModule_AddObject(module, "DanilovEnvelopeSolver22", (PyObject *)&pyORBIT_DanilovEnvelopeSolver22_Type);
}

#ifdef __cplusplus
}
#endif

// end of namespace wrap_danilov_envelope
} // namespace wrap_danilov_envelope