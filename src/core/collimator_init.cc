# include "orbit/MaterialInteractions/wrap_collimator.hh"
PyMODINIT_FUNC PyInit_collimator(void) {
    return wrap_collimator::initcollimator();
}
