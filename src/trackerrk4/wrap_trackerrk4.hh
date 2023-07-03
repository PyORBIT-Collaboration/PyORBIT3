#ifndef WRAP_RK4_TRACKER_H
#define WRAP_RK4_TRACKER_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

	PyMODINIT_FUNC inittrackerrk4(void);
	PyObject* getTrackerRK4Type(const char* name);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_RK4_TRACKER_H
