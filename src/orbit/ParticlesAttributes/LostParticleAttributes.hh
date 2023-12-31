/////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   LostParticleAttributes.hh
//
// AUTHOR
//    S. Cousineau
//
// CREATED
//    12/08/2011
//
// DESCRIPTION
//    A subclass of the particle attributes class.
//
///////////////////////////////////////////////////////////////////////////
#ifndef LOSTPARTICLEATTRIBUTES_HH_
#define LOSTPARTICLEATTRIBUTES_HH_

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include "ParticleAttributes.hh"

class LostParticleAttributes : public ParticleAttributes
{
public:

	/** This Attribute class contains additional properties of lost particles.
	*/
	LostParticleAttributes(Bunch* bunch);
	~LostParticleAttributes();
	double& getPosition(int particle_index);
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////


#endif
