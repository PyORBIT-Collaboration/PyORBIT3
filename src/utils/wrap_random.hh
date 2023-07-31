#ifndef WRAP_ORBIT_UTILS_RANDOM_HH_
#define WRAP_ORBIT_UTILS_RANDOM_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_random{
    void initRandom(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_UTILS_RANDOM_HH_*/
