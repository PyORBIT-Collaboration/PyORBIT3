#ifndef WRAP_ENVELOPE_MODULE_HH_
#define WRAP_ENVELOPE_MODULE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace wrap_envelope {
  // The name of the function was changed to avoid collision with PyImport magic naming.
  PyMODINIT_FUNC initenvelope();
}

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ENVELOPE_MODULE_HH_*/