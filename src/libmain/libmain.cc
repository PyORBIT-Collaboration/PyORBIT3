#include <Python.h>

#include "wrap_orbit_mpi.hh"
#include "wrap_bunch.hh"
#include "wrap_spacecharge.hh"
#include "wrap_linacmodule.hh"
#include "wrap_teapotbase.hh"
#include "wrap_utils.hh"
#include "wrap_aperture.hh"

// Define the module methods for the orbit module
static PyMethodDef OrbitMethods[] = {
    // Add methods here
    {nullptr, nullptr, 0, nullptr}
};


// Define the orbit module definition structure
static PyModuleDef OrbitModule = {
    PyModuleDef_HEAD_INIT,
    "orbit",
    "pORBIT modules implemented in C++ .",
    -1,
    OrbitMethods,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

// Define the module initialization function
// It is a top level placeholder, has no methods or anything
PyMODINIT_FUNC PyInit__orbit(void) {
    PyObject* orbit_module = PyModule_Create(&OrbitModule);
    if (!orbit_module) {
        return nullptr;
    }
    return orbit_module;
}


// Create wrappers with magic names
// that will be picked up by importlib

PyMODINIT_FUNC PyInit_orbit_mpi(void) {
    return wrap_orbit_mpi::initorbit_mpi();
}

PyMODINIT_FUNC PyInit_bunch(void) {
    return wrap_orbit_bunch::initbunch();
}

PyMODINIT_FUNC PyInit_spacecharge(void) {
    return initspacecharge();
}

PyMODINIT_FUNC PyInit_teapot_base(void) {
    return wrap_teapotbase::initteapotbase();
}

PyMODINIT_FUNC PyInit_linac(void) {
    return wrap_linac::initlinac();
}

PyMODINIT_FUNC PyInit_orbit_utils(void) {
    return wrap_orbit_utils::initutils();
}

PyMODINIT_FUNC PyInit_aperture(void) {
    return wrap_aperture::initaperture();
}

