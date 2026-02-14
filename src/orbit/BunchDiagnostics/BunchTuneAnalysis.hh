#ifndef BUNCH_TUNE_ANALYSIS_H
#define BUNCH_TUNE_ANALYSIS_H

#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BunchTwissAnalysis.hh"

using namespace std;


/** Estimates particle tunes using average phase advance (APA) over one turn. */
class BunchTuneAnalysis: public OrbitUtils::CppPyWrapper
{
	public:
		/** Constructor*/
		BunchTuneAnalysis();

		/** Destructor */
		virtual ~BunchTuneAnalysis();

		/** Estimates tunes. */
		void analyzeBunch(Bunch* bunch);

		/** Sets element of normalization matrix. */
		void setNormMatrixElement(int i, int j, double value);

		/** Returns element of normalization matrix. */
		double getNormMatrixElement(int i, int j);

		/** Sets normalization matrix based on uncoupled Twiss parameters. */
		void assignTwiss(double betax, double alphax, double etax, double etapx, double betay, double alphay);

	private:
		// Normalization matrix V^{-1}
		double matrix[6][6];
};


#endif
//endif for BUNCH_TUNE_ANALYSIS_H
