#include "wrap_aperture.hh"
PyMODINIT_FUNC PyInit_aperture(void) {
    return wrap_aperture::initaperture();
}