#ifndef RANDOM_H_
#define RANDOM_H_

//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    Random.hh
//
// AUTHOR
//    S. Cousineau
//
// CREATED
//    09/09/2011
//
// DESCRIPTION
//    This class for generating a random number between 0 and 1.
//
///////////////////////////////////////////////////////////////////////////
#include <complex>
#include <cmath>
#include <iostream>

/**
  This class is a methods to create random numbers
*/

namespace OrbitUtils
{
	class Random
	{
	public:
		/** Method to initialise the random seed */
		static void seed(int seed);

		/** The method calculates a random number between 0 and 1 */
		static double ran1();
	};
}

#endif
