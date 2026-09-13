#ifndef WRAP_ORBIT_MPI_TYPES_H
#define WRAP_ORBIT_MPI_TYPES_H

#include <Python.h>

#include "mpi/orbit_mpi.hh"

typedef struct {
  PyObject_HEAD
  MPI_Comm comm;
} pyORBIT_MPI_Comm;

typedef struct {
  PyObject_HEAD
  MPI_Group group;
} pyORBIT_MPI_Group;

typedef struct {
  PyObject_HEAD
  MPI_Status status;
} pyORBIT_MPI_Status;

typedef struct {
  PyObject_HEAD
  MPI_Request request;
} pyORBIT_MPI_Request;

typedef struct {
  PyObject_HEAD
  MPI_Datatype datatype;
} pyORBIT_MPI_Datatype;

typedef struct {
  PyObject_HEAD
  MPI_Op op;
} pyORBIT_MPI_Op;

#endif
