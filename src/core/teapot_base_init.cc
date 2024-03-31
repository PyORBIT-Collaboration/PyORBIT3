#include <Python.h>
#include "wrap_orbit_mpi.hh"
#include "wrap_teapotbase.hh"
#include "wrap_errorbase.hh"

PyMODINIT_FUNC PyInit_teapot_base(void) {
    return wrap_teapotbase::initteapotbase();
}