#ifndef WRAP_ORBIT_UTILS_NUMRECIPES_HH_
#define WRAP_ORBIT_UTILS_NUMRECIPES_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_numrecipes{
    void initNumrecipes(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_UTILS_NUMRECIPES_HH_*/
