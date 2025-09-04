#ifndef BUNCH_TUNE_ANALYSIS_H
#define BUNCH_TUNE_ANALYSIS_H

#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BunchTwissAnalysis.hh"

using namespace std;


/**
  The BunchTuneAnalysis class calculates the particle tunes
*/
class BunchTuneAnalysis: public OrbitUtils::CppPyWrapper
{
	public:

		/** Constructor*/
		BunchTuneAnalysis();

		/** Destructor */
		virtual ~BunchTuneAnalysis();

		/** Performs the Twiss analysis of the bunch */
		void analyzeBunch(Bunch* bunch);

		//** Assigns Twiss values at location of calculator */
		void assignTwiss(double bx, double ax, double dx, double dpx, double by, double ay);


	private:
		//** Twiss */
		double betax;
		double alphax;
		double etax;
		double etapx;
		double betay;
		double alphay;

};


/**
  Calculates tunes using 4D normalization.
*/
class BunchTuneAnalysis4D: public OrbitUtils::CppPyWrapper
{
	public:
		BunchTuneAnalysis4D();

		virtual ~BunchTuneAnalysis4D();

		void analyzeBunch(Bunch* bunch);

		void setMatrixElement(int i, int j, double value);

		double getMatrixElement(int i, int j);

	private:
		double matrix[4][4];
};


#endif
//endif for BUNCH_TUNE_ANALYSIS_H
