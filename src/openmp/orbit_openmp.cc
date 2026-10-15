#include "orbit_openmp.hh"

#ifdef WITH_OPENMP
#include <omp.h>

int orbit_omp_get_max_threads()
{
  return omp_get_max_threads();
}

int orbit_omp_get_num_threads()
{
  return omp_get_num_threads();
}

int orbit_omp_get_thread_num()
{
  return omp_get_thread_num();
}

#else

int orbit_omp_get_max_threads()
{
  return 1;
}

int orbit_omp_get_thread_num()
{
  return 0;
}

int orbit_omp_get_num_threads()
{
  return 1;
}

#endif