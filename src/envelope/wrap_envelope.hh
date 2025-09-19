#ifndef WRAP_ENVELOPE_H
#define WRAP_ENVELOPE_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

	PyMODINIT_FUNC initenvelope(void);
	PyObject* getEnvelopeType(const char* name);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_ENVELOPE_H