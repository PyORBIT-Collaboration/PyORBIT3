//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   Random.cc
//
// AUTHOR
//   S. Cousineau
//   Changed by A. Shishlo 2022.02.14
//
// CREATED
//    09/09/2011
//
// DESCRIPTION
//    This class contains a standart c++ method for calculating
//    a random number between 0 and 1.
///////////////////////////////////////////////////////////////////////////

#include "Random.hh"
#include <random>

using namespace OrbitUtils;

// By default the random generator is initialized using a random seed (current time)
std::mt19937 mt(time(0));

void Random::seed(int seed){
	mt.seed(seed);
}

/** The method calculates a random number between 0 and 1. (0. and 1.0 excluded) */
double Random::ran1(){
	return ((double) mt() / (mt.max()));
}
