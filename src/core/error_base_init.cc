#include <Python.h>
# include "orbit/Errors/wrap_errorbase.hh"
PyMODINIT_FUNC PyInit_error_base(void) {
    return wrap_errorbase::initerrorbase();
}