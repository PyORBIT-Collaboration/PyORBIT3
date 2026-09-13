#ifndef WRAP_SPACE_CHARGE_GRID_2D_H
#define WRAP_SPACE_CHARGE_GRID_2D_H

#include <Python.h>

class Grid2D;

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initGrid2D(PyObject* module);
    PyObject* wrapGrid2D(Grid2D* grid, PyObject* owner);
  }

#ifdef __cplusplus
}
#endif

#endif
