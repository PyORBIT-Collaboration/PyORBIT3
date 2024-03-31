#include "wrap_linacmodule.hh"
PyMODINIT_FUNC PyInit_linac(void) {
    return wrap_linac::initlinac();
}