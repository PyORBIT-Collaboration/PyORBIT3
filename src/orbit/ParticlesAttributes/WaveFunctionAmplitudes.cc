//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//  WaveFunctionAmplitudes.cc
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    07/19/2005
//
// DESCRIPTION
//    A subclass of a ParticleAttributes class
//
//
///////////////////////////////////////////////////////////////////////////

#include "orbit/Bunch.hh"
#include "orbit/ParticlesAttributes/WaveFunctionAmplitudes.hh"

WaveFunctionAmplitudes::WaveFunctionAmplitudes(Bunch* bunch, int size_in):
  ParticleAttributes(bunch, size_in)
{
  cl_name_ = "Amplitudes";
  attrDescr = "Ampls";
}

WaveFunctionAmplitudes::~WaveFunctionAmplitudes()
{
}
