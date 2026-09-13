#include <Python.h>
#include "mpi/wrap_orbit_mpi.hh"
#include "teapot/wrap_teapotbase.hh"
#include "orbit/Errors/wrap_errorbase.hh"

PyMODINIT_FUNC PyInit_teapot_base(void) {
    return wrap_teapotbase::initteapotbase();
}