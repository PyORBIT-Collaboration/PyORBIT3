#ifndef WRAP_DANILOV_ENVELOPE_TRACKER_H
#define WRAP_DANILOV_ENVELOPE_TRACKER_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace wrap_envelope {
void initDanilovEnvelopeTracker(PyObject *module);
}

#ifdef __cplusplus
}
#endif

#endif
