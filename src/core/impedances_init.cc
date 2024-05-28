#include "wrap_impedances.hh"
PyMODINIT_FUNC PyInit_impedances(void) {
    return wrap_impedances::initimpedances();
}