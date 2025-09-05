#ifndef WRAP_DANILOV_ENVELOPE_MODULE_HH_
#define WRAP_DANILOV_ENVELOPE_MODULE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace wrap_danilov_envelope {
  // The name of the function was changed to avoid collision with PyImport magic naming.
  PyMODINIT_FUNC initdanilovenvelope();
}

#ifdef __cplusplus
}
#endif

#endif /*WRAP_DANILOV_ENVELOPE_MODULE_HH_*/