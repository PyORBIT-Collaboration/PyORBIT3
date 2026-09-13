//The base class for BaseApertures. It defines the interface for BaseAperture
#ifndef BASE_APERTURE_H
#define BASE_APERTURE_H

#include "orbit/Bunch.hh"
#include "orbit/Apertures/BaseApertureShape.hh"


///////////////////////////////////////////////////////////////////////////
//
// NAME
//
//   BaseAperture
//
// AUTHOR:
//   Andrei Shishlo October 2022
//
///////////////////////////////////////////////////////////////////////////

/**
   BaseAperture class defines actions with macro particles in a main bunch
   and a lost bunch instances with respect of transverse coordinates of
   macro particles. It can keep macro particle in the main bunch or remove it
   after asking the BaseApertureShape class instance about the suitability
   of particle's coordinates.
*/

class BaseAperture
{
public:

	/** BaseAperture constructor */
	BaseAperture();

	/** BaseAperture decstructor */
	virtual ~BaseAperture();

	/** Returns aperture shape */
	BaseApertureShape* getApertureShape();

	/** Sets a borrowed aperture shape. The caller retains ownership. */
	void setApertureShape(BaseApertureShape* apertureShape);

	/** Routine for transfering particles through a aperture */
	void checkBunch(Bunch* bunch, Bunch* lostbunch);

	/** Returns the total number of partciles lost across all CPUs */
	int getNumberOfLost();

	/** Returns the aperture name */
	std::string getName();

	/** Sets the aperture name */
	void setName(std::string apertureNameIn);

	/** Sets the position of the node in the lattice */
	double getPosition();

	/** Returns the position of the node in the lattice */
	void setPosition(double position);

	/**
	Sets the aperture in an active ( 1 )/ not active ( 0 ) state
	*/
	void setOnOff(int isActive);

	/**
	Returns the aperture in an active ( 1 )/ not active ( 0 ) state
	*/
	int getOnOff();

protected:

	//name of the aperture
	std::string apertureName;

	//Counters
	int nLost_;

	//aperture position
	double pos_;

	//BaseApertureShape
	BaseApertureShape* apertureShape;

	//Info variable defining if the aperture is active or not.
	//isActive == 1 then it is active otherwise it is not.
	int isActive;

};

//end of BASE_APERTURE_H ifdef
#endif
