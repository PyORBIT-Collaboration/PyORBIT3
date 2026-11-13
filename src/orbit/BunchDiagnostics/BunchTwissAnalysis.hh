#ifndef BUNCH_TWISS_ANALYSIS_H
#define BUNCH_TWISS_ANALYSIS_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Bunch.hh"

#include <vector>


/**
  The BunchTwissAnalysis class calculates the average of 6D coordinates and they correlations.
  As results it returns the Twiss parameters for each plane.
*/

class BunchTwissAnalysis: public OrbitUtils::CppPyWrapper
{
	public:

		/** Constructor*/
		BunchTwissAnalysis();

		/** Destructor */
		virtual ~BunchTwissAnalysis();

		/** Performs the Twiss analysis of the bunch */
		void analyzeBunch(Bunch* bunch);

		/**
		 * @brief Returns the centered correlation between two components.
		 *
		 * \f$\langle (x-\langle x\rangle)(y-\langle y\rangle)\rangle
		 * = \langle x y\rangle - \langle x\rangle\langle y\rangle\f$
		 *
		 * @param ic Index of the first component (x).
		 * @param jc Index of the second component (y).
		 * @return The centered correlation value as a double.
		 */
		double getCorrelation(int ic, int jc) const;

		/** Returns the average value for coordinate with index ic */
		double getAverage(int ic) const;

		/** Returns the total number of analysed macroparticles */
		int getGlobalCount() const;

		/** Returns the total macrosize */
		double getGlobalMacrosize() const;

		/** Returns the emittance for index 0,1,2 - x,y,z planes. */
		double getEffectiveEmittance(int ic) const;

		/** Returns Twiss alpha for index 0,1,2 - x,y,z planes.*/
		double getEffectiveAlpha(int ic) const;

		/** Returns Twiss beta for index 0,1,2 - x,y,z planes.*/
		double getEffectiveBeta(int ic) const;

		/** Returns Twiss gamma for index 0,1,2 - x,y,z planes.*/
		double getEffectiveGamma(int ic) const;

		/** Returns the pure betatron emittance for index 0,1 - x,y planes. */
		double getEmittance(int ic) const;

		/** Returns the normalized pure betatron emittance for index 0,1 - x,y planes. */
		double getEmittanceNormalized(int ic) const;

		/** Returns Twiss alpha without dispersive part for index 0,1 - x,y planes.*/
		double getAlpha(int ic) const;

		/** Returns Twiss beta without dispersive part for index 0,1 - x,y planes.*/
		double getBeta(int ic) const;

		/** Returns Twiss gamma without dispersive part for index 0,1 - x,y planes.*/
		double getGamma(int ic) const;

		/** Returns Twiss dispersion function for index 0,1 - x,y planes.*/
		double getDispersion(int ic) const;

		/** Returns Twiss dispersion_prime function for index 0,1 - x,y planes.*/
		double getDispersionDerivative(int ic) const;

		/** Computes the XY moments of the bunch up to a prescribed order */
		void computeBunchMoments(Bunch* bunch, int order, int dispersionflag, int emitnormflag);

		/** Returns the XY moment of the beam */
		double getBunchMoment(int i, int j) const;

	private:

		/** Number of points accounted */
		int count;

		/** total macrosize accounted */
		double total_macrosize;

		/** array with average values for 6D coordinates */
		double avg_arr[6] = {};

		/** array with correlations between 6D coordinates.
		    It is a packed 1D array 6*6 = 36.
				To reach i,j = 0-5 coerrelation use [i+6*j]
		    It is excessive because <x*y> = <y*x> etc.*/
		double corr_arr[36] = {};

		/** Contiguous XY moments, indexed by j * (_order + 1) + i. */
		std::vector<double> momentXY;

		/** order for the moments **/
		int _order;

		/** bunch parameters */
		double bunch_momentum;
		double bunch_gamma;
		double bunch_beta;
	    double bunch_kinenergy;
		double bunch_mass;

};

#endif
//endif for BUNCH_TWISS_ANALYSIS_H
