#ifndef WRAP_KV_ENVELOPE_TRACKER_H
#define WRAP_KV_ENVELOPE_TRACKER_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace wrap_envelope {
void initKVEnvelopeSolver(PyObject *module);
}

#ifdef __cplusplus
}
#endif

#endif
