#ifndef WRAP_UTILS_BUNCH_FUNCTIONS_H
#define WRAP_UTILS_BUNCH_FUNCTIONS_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_utils_bunch_functions{
    void initBunchUtilsFunctions(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
