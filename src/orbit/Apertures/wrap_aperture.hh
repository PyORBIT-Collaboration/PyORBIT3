#ifndef WRAP_ORBIT_APERTURE_MODULE_HH_
#define WRAP_ORBIT_APERTURE_MODULE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_aperture{
    PyMODINIT_FUNC PyInit_aperture();
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_APERTURE_MODULE_HH_*/
