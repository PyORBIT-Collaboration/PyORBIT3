#include "Python.h"
#include "orbit_mpi.hh"

#include <cstdlib>

#include <iostream>

//modules headers
#include "wrap_orbit_mpi.hh"
#include "wrap_bunch.hh"
#include "wrap_utils.hh"
#include "wrap_teapotbase.hh"
//#include "wrap_errorbase.hh"
//#include "wrap_trackerrk4.hh"
#include "wrap_spacecharge.hh"
#include "wrap_linacmodule.hh"
//#include "wrap_collimator.hh"
//#include "wrap_foil.hh"
//#include "wrap_rfcavities.hh"
//#include "wrap_aperture.hh"
//#include "wrap_fieldtracker.hh"
//#include "wrap_impedances.hh"

/**
 * The main function that will initialize the MPI and will
 * call the python interpreter: Py_Main(argc,argv).
 */

int main(int argc, char **argv)
{
  //  for(int i = 0; i < argc; i++)
  //  {
  //    std::cout << "before i = " << i << " arg = " << argv[i] << std::endl;
  //  }

  ORBIT_MPI_Init(&argc, &argv);
  
  wchar_t *program = Py_DecodeLocale(argv[0], NULL);

  //  for(int i = 0; i < argc; i++)
  //  {
  //    std::cout << "after i = " << i <<" arg = " << argv[i] << std::endl;
  //  }

  //  int rank = 0;
  //  int size = 0;
  //  ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //  ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size);
  //  std::cout << "rank = " << rank << " size = " << size << std::endl;

  
  
  /* Add a built-in module bunch, before Py_Initialize */
  if (PyImport_AppendInittab("bunch", wrap_orbit_bunch::PyInit_bunch) == -1) {
  	fprintf(stderr, "Error: could not extend in-built modules table. bunch\n");
  	exit(1);
  }
    
  /* Add a built-in module bunch, before Py_Initialize */
  if (PyImport_AppendInittab("orbit_mpi", wrap_orbit_mpi::initorbit_mpi) == -1) {
  	fprintf(stderr, "Error: could not extend in-built modules table. orbit_mpi\n");
  	exit(1);
  }
  
  /* Add a built-in module bunch, before Py_Initialize */
  if (PyImport_AppendInittab("spacecharge", initspacecharge) == -1) {
  	fprintf(stderr, "Error: could not extend in-built modules table. spacecharge\n");
  	exit(1);
  }

  /* Add a built-in module bunch, before Py_Initialize */
  if (PyImport_AppendInittab("teapot_base", wrap_teapotbase::initteapotbase) == -1) {
  	fprintf(stderr, "Error: could not extend in-built modules table. teapot_base\n");
  	exit(1);
  } 
  
  /* Add a built-in module bunch, before Py_Initialize */
  if (PyImport_AppendInittab("orbit_utils", wrap_orbit_utils::initutils) == -1) {
  	fprintf(stderr, "Error: could not extend in-built modules table. orbit_utils\n");
  	exit(1);
  } 
  
  /* Add a built-in module bunch, before Py_Initialize */
  if (PyImport_AppendInittab("linac", wrap_linac::initlinac) == -1) {
  	fprintf(stderr, "Error: could not extend in-built modules table. linac\n");
  	exit(1);
  } 

  // We need to initialize the extra ORBIT modules
  Py_Initialize();

  // ORBIT module initializations
  //wrap_errorbase::initerrorbase();
  //wrap_collimator::initcollimator();
  //wrap_aperture::initaperture();
  //wrap_foil::initfoil();
  //wrap_rfcavities::initrfcavities();
  //wrap_fieldtracker::initfieldtracker();
  //wrap_impedances::initimpedances();

  // Runge-Kutta tracker package

  //inittrackerrk4();

  // The python interpreter
  // It will call Py_Initialize() again, but there is no harm.

  
  wchar_t** wargv = new wchar_t*[argc];
  
  for(int i = 0; i < argc; i++)
  {
  	wargv[i] = Py_DecodeLocale(argv[i], NULL);
  	if(wargv[i] == NULL)
  	{
  		return EXIT_FAILURE;
  	}
  }
  
  //Py_SetProgramName(wargv[0]);

  //PySys_SetArgv(argc, wargv);

  Py_Main(argc, wargv);

  ORBIT_MPI_Finalize();
  
  for(int i = 0; i < argc; i++)
  {
  	PyMem_RawFree(wargv[i]);
  	wargv[i] = NULL;
  }    

  return 0;
}
