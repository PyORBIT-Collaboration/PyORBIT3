#include "wrap_field_sources_module.hh"
PyMODINIT_FUNC PyInit_field_sources(void) {
    return wrap_field_sources_module::initFieldSourcesModule();
}