#ifndef WRAP_SPACE_CHARGE_GRID_H
#define WRAP_SPACE_CHARGE_GRID_H

#include <Python.h>

typedef struct {
  PyObject_HEAD
  void* cpp_obj;
  int owns_cpp;
  PyObject* owner;
} pyORBIT_Grid;

#endif
