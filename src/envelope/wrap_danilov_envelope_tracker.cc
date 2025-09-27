#include <iostream>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "DanilovEnvelopeTracker.hh"
#include "wrap_bunch.hh"
#include "wrap_danilov_envelope_tracker.hh"
#include "wrap_envelope.hh"

namespace wrap_envelope {

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for Python class wrapping DanilovEnvelopeTracker instance.
// It never will be called directly.
static PyObject *DanilovEnvelopeTracker_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    pyORBIT_Object *self;
    self = (pyORBIT_Object *)type->tp_alloc(type, 0);
    self->cpp_obj = NULL;
    return (PyObject *)self;
}

// Initialization of Python DanilovEnvelopeTracker class.
// This is implementation of the __init__ method.
static int DanilovEnvelopeTracker_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds) {
    double perveance = 0.0;
    self->cpp_obj = new DanilovEnvelopeTracker(perveance);
    ((DanilovEnvelopeTracker *)self->cpp_obj)->setPyWrapper((PyObject *)self);
    return 0;
}

// Method: trackBunch(Bunch* bunch, double length)
static PyObject *DanilovEnvelopeTracker_trackBunch(PyObject *self, PyObject *args) {
    pyORBIT_Object *pyDanilovEnvelopeTracker = (pyORBIT_Object *)self;
    DanilovEnvelopeTracker *cpp_DanilovEnvelopeTracker = (DanilovEnvelopeTracker *)pyDanilovEnvelopeTracker->cpp_obj;
    PyObject *pyBunch;
    double length;
    if (!PyArg_ParseTuple(args, "Od:trackBunch", &pyBunch, &length)) {
        ORBIT_MPI_Finalize("PyDanilovEnvelopeTracker - trackBunch(Bunch* bunch, double length) - parameters are needed.");
    }
    PyObject *pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
    if (!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type)) {
        ORBIT_MPI_Finalize("PyDanilovEnvelopeTracker - trackBunch(Bunch* bunch, double length) - first parameter should be Bunch.");
    }
    Bunch *cpp_bunch = (Bunch *)((pyORBIT_Object *)pyBunch)->cpp_obj;
    cpp_DanilovEnvelopeTracker->trackBunch(cpp_bunch, length);
    Py_INCREF(Py_None);
    return Py_None;
}

// Method: setPerveance(double perveance)
static PyObject *DanilovEnvelopeTracker_setPerveance(PyObject *self, PyObject *args) {
    pyORBIT_Object *pyDanilovEnvelopeTracker = (pyORBIT_Object *)self;
    DanilovEnvelopeTracker *cpp_DanilovEnvelopeTracker = (DanilovEnvelopeTracker *)pyDanilovEnvelopeTracker->cpp_obj;
    double perveance;
    if (!PyArg_ParseTuple(args, "d:setPerveance", &perveance)) {
        ORBIT_MPI_Finalize("PyDanilovEnvelopeTracker - setPerveance(double perveance) - parameters are needed.");
    }
    cpp_DanilovEnvelopeTracker->setPerveance(perveance);
    Py_INCREF(Py_None);
    return Py_None;
}

// Method: getPerveance()
static PyObject *DanilovEnvelopeTracker_getPerveance(PyObject *self, PyObject *args) {
    pyORBIT_Object *pyDanilovEnvelopeTracker = (pyORBIT_Object *)self;
    DanilovEnvelopeTracker *cpp_DanilovEnvelopeTracker = (DanilovEnvelopeTracker *)pyDanilovEnvelopeTracker->cpp_obj;
    double perveance = cpp_DanilovEnvelopeTracker->getPerveance();
    return Py_BuildValue("d", perveance);
}

// Destructor for python DanilovEnvelopeTracker class (__del__ method)
static void DanilovEnvelopeTracker_del(pyORBIT_Object *self) {
    DanilovEnvelopeTracker *cpp_DanilovEnvelopeTracker = (DanilovEnvelopeTracker *)self->cpp_obj;
    delete cpp_DanilovEnvelopeTracker;
    self->ob_base.ob_type->tp_free((PyObject *)self);
}

// Definition of Python DanilovEnvelopeTracker wrapper class methods.
// They will be available from the Python level.
static PyMethodDef DanilovEnvelopeTrackerClassMethods[] = {
    {"getPerveance", DanilovEnvelopeTracker_getPerveance, METH_VARARGS, "Get space charge perveance."},
    {"setPerveance", DanilovEnvelopeTracker_setPerveance, METH_VARARGS, "Set space charge perveance."},
    {"trackBunch", DanilovEnvelopeTracker_trackBunch, METH_VARARGS, "Apply space charge kick to beam envelope."},
    {NULL}
};

// Definition of Python DanilovEnvelopeTracker wrapper class members.
// They will be available from the Python level.
static PyMemberDef DanilovEnvelopeTrackerClassMembers[] = {{NULL}};

// New python DanilovEnvelopeTracker wrapper type definition.
static PyTypeObject pyORBIT_DanilovEnvelopeTracker_Type = {
    PyVarObject_HEAD_INIT(NULL, 0) "DanilovEnvelopeTracker", /*tp_name*/
    sizeof(pyORBIT_Object),                                  /*tp_basicsize*/
    0,                                                       /*tp_itemsize*/
    (destructor)DanilovEnvelopeTracker_del,                  /*tp_dealloc*/
    0,                                                       /*tp_print*/
    0,                                                       /*tp_getattr*/
    0,                                                       /*tp_setattr*/
    0,                                                       /*tp_compare*/
    0,                                                       /*tp_repr*/
    0,                                                       /*tp_as_number*/
    0,                                                       /*tp_as_sequence*/
    0,                                                       /*tp_as_mapping*/
    0,                                                       /*tp_hash */
    0,                                                       /*tp_call*/
    0,                                                       /*tp_str*/
    0,                                                       /*tp_getattro*/
    0,                                                       /*tp_setattro*/
    0,                                                       /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,                /*tp_flags*/
    "The DanilovEnvelopeTracker python wrapper",             /* tp_doc */
    0,                                                       /* tp_traverse */
    0,                                                       /* tp_clear */
    0,                                                       /* tp_richcompare */
    0,                                                       /* tp_weaklistoffset */
    0,                                                       /* tp_iter */
    0,                                                       /* tp_iternext */
    DanilovEnvelopeTrackerClassMethods,                      /* tp_methods */
    DanilovEnvelopeTrackerClassMembers,                      /* tp_members */
    0,                                                       /* tp_getset */
    0,                                                       /* tp_base */
    0,                                                       /* tp_dict */
    0,                                                       /* tp_descr_get */
    0,                                                       /* tp_descr_set */
    0,                                                       /* tp_dictoffset */
    (initproc)DanilovEnvelopeTracker_init,                   /* tp_init */
    0,                                                       /* tp_alloc */
    DanilovEnvelopeTracker_new,                              /* tp_new */
};

// Initialization function of the pyDanilovEnvelopeTracker class
void initDanilovEnvelopeTracker(PyObject *module) {
    if (PyType_Ready(&pyORBIT_DanilovEnvelopeTracker_Type) < 0) {
        return;
    }
    Py_INCREF(&pyORBIT_DanilovEnvelopeTracker_Type);
    PyModule_AddObject(module, "DanilovEnvelopeTracker", (PyObject *)&pyORBIT_DanilovEnvelopeTracker_Type);
}

#ifdef __cplusplus
}
#endif

} // namespace wrap_envelope
