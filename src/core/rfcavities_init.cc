#include "wrap_rfcavities.hh"
PyMODINIT_FUNC PyInit_rfcavities(void) {
    return wrap_rfcavities::initrfcavities();
}