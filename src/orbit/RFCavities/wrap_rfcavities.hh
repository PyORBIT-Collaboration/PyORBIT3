#ifndef WRAP_RFCAVITIES_H
#define WRAP_RFCAVITIES_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_rfcavities
  {
    PyMODINIT_FUNC initrfcavities();
    PyObject* getRFCavityType(char* name);
  }

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif // WRAP_RFCAVITIES_H

