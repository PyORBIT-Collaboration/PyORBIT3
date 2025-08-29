#ifndef BUNCH_TUNE_ANALYSIS_4D_H
#define BUNCH_TUNE_ANALYSIS_4D_H

#include "CppPyWrapper.hh"

#include "Bunch.hh"
#include "BunchTwissAnalysis.hh"

using namespace std;


class BunchTuneAnalysis4D: public OrbitUtils::CppPyWrapper
{
	public:
		BunchTuneAnalysis4D();

		virtual ~BunchTuneAnalysis4D();

		void analyzeBunch(Bunch* bunch);

		// void setMatrix(double _matrix[4][4]);
		void setMatrix(double dummy);

		double getTune(int ic);

	private:
		double matrix[4][4];
};


#endif
//endif for BUNCH_TUNE_ANALYSIS_4D_H
