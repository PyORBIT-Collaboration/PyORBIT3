#include "wrap_fieldtracker.hh"
PyMODINIT_FUNC PyInit_fieldtracker(void) {
    return wrap_fieldtracker::initfieldtracker();
}