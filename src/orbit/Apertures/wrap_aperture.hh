#ifndef WRAP_ORBIT_APERTURE_MODULE_HH_
#define WRAP_ORBIT_APERTURE_MODULE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_aperture{
  /* The name of the function was changed to avoid collision with PyImport magic naming */
    PyMODINIT_FUNC initaperture();
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_ORBIT_APERTURE_MODULE_HH_*/
