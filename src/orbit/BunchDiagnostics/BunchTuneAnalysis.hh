#ifndef BUNCH_TUNE_ANALYSIS_H
#define BUNCH_TUNE_ANALYSIS_H

//pyORBIT utils

#include "orbit/Bunch.hh"
#include "orbit/BunchDiagnostics/BunchTwissAnalysis.hh"



/** Estimates particle tunes using average phase advance (APA) over one turn. */
class BunchTuneAnalysis
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
		void setNormMatrixFromTwiss(double betax, double alphax, double etax, double etapx, double betay, double alphay);

	private:
		// Normalization matrix V^{-1}
		double matrix[6][6];

		// Flag to erase stored phase/amplitude info.
		int erase;
};


#endif
//endif for BUNCH_TUNE_ANALYSIS_H
