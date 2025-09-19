#ifndef WRAP_DANILOV_22_ENVELOPE_SOLVER_H
#define WRAP_DANILOV_22_ENVELOPE_SOLVER_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace wrap_envelope {
void initDanilov22EnvelopeSolver(PyObject *module);
}

#ifdef __cplusplus
}
#endif

#endif