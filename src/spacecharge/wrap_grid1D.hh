#ifndef WRAP_SPACE_CHARGE_GRID_1D_H
#define WRAP_SPACE_CHARGE_GRID_1D_H

#include <Python.h>

class Grid1D;

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_spacecharge
  {
    void initGrid1D(PyObject* module);
    PyObject* wrapGrid1D(Grid1D* grid, PyObject* owner);
  }

#ifdef __cplusplus
}
#endif

#endif
