#ifndef WRAP_UNIFORM_ELLIPSOID_FIELD_CALC_H
#define WRAP_UNIFORM_ELLIPSOID_FIELD_CALC_H

#include <Python.h>

class UniformEllipsoidFieldCalculator;

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initUniformEllipsoidFieldCalculator(PyObject* module);
    PyObject* wrapUniformEllipsoidFieldCalculator(
      UniformEllipsoidFieldCalculator* calculator,
      PyObject* owner
    );
  }

#ifdef __cplusplus
}
#endif

#endif
