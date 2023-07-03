#include <Python.h>

#include "wrap_orbit_mpi.hh"
#include "wrap_bunch.hh"
#include "wrap_spacecharge.hh"
#include "wrap_linacmodule.hh"
#include "wrap_teapotbase.hh"
#include "wrap_utils.hh"
#include "wrap_aperture.hh"
#include "wrap_foil.hh"
#include "wrap_collimator.hh"
#include "wrap_errorbase.hh"
#include "wrap_trackerrk4.hh"
#include "wrap_field_sources_module.hh"

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

PyMODINIT_FUNC PyInit_trackerrk4(void) {
    return inittrackerrk4();
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

PyMODINIT_FUNC PyInit_error_base(void) {
    return wrap_errorbase::initerrorbase();
}

PyMODINIT_FUNC PyInit_collimator(void) {
    return wrap_collimator::initcollimator();
}

PyMODINIT_FUNC PyInit_foil(void) {
    return wrap_foil::initfoil();
}

PyMODINIT_FUNC PyInit_field_sources(void) {
    return wrap_field_sources_module::initFieldSourcesModule();
}
