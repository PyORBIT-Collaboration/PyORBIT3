#ifndef WRAP_IMPEDANCES_H
#define WRAP_IMPEDANCES_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_impedances
  {
    PyMODINIT_FUNC initimpedances();
    PyObject* getImpedanceType(char* name);
  }

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_IMPEDANCES_H
