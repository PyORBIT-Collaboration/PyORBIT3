#include "wrap_foil.hh"
PyMODINIT_FUNC PyInit_foil(void) {
    return wrap_foil::initfoil();
}
