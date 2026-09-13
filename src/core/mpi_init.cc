#include <Python.h>
#include "mpi/wrap_orbit_mpi.hh"

PyMODINIT_FUNC PyInit_orbit_mpi(void) {
    return wrap_orbit_mpi::initorbit_mpi();
}
