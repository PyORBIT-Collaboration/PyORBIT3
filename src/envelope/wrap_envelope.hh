#ifndef WRAP_ENVELOPE_MODULE_HH_
#define WRAP_ENVELOPE_MODULE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace wrap_envelope {
PyMODINIT_FUNC initenvelope();
}

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ENVELOPE_MODULE_HH_*/
