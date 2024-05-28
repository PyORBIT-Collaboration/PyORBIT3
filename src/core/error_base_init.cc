#include <Python.h>
# include "wrap_errorbase.hh"
PyMODINIT_FUNC PyInit_error_base(void) {
    return wrap_errorbase::initerrorbase();
}