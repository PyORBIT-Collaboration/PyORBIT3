#include "wrap_bunch.hh"
PyMODINIT_FUNC PyInit_bunch(void) {
    return wrap_orbit_bunch::initbunch();
}