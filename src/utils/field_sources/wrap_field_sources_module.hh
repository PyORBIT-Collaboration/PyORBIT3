#ifndef WRAP_FIELD_SOURCES_MODULE_H
#define WRAP_FIELD_SOURCES_MODULE_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_field_sources_module{
    PyMODINIT_FUNC initFieldSourcesModule();
  }

#ifdef __cplusplus
}
#endif

#endif
