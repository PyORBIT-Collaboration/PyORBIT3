#include <Python.h>
#include "wrap_orbit_mpi.hh"
#include "wrap_danilov_envelope.hh"

PyMODINIT_FUNC PyInit_danilov_envelope(void) {
    return wrap_danilov_envelope::initdanilovenvelope();
}