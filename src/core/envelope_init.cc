#include "wrap_envelope.hh"
PyMODINIT_FUNC PyInit_envelope(void) {
    return wrap_envelope::initenvelope();
}