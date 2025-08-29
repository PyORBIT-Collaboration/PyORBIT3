#include "BunchTuneAnalysis4D.hh"
#include "SyncPart.hh"
#include "OrbitConst.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>


BunchTuneAnalysis4D::BunchTuneAnalysis4D(): CppPyWrapper(NULL)
{
	double matrix[4][4] = {
		{1.0, 0.0, 0.0, 0.0},
		{0.0, 1.0, 0.0, 0.0},
		{0.0, 0.0, 1.0, 0.0},
		{0.0, 0.0, 0.0, 1.0}
	};
}


BunchTuneAnalysis4D::~BunchTuneAnalysis4D()
{
}


void BunchTuneAnalysis4D::setMatrix(double dummy) {	
	// for (int i = 0; i < 4; i++) {
	// 	for (int j = 0; j < 4; j++) {
	// 		matrix[i][j] = _matrix[i][j];
	// 	}
	// }
}


void BunchTuneAnalysis4D::analyzeBunch(Bunch* bunch){

	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double beta = syncPart->getBeta();
	double** part_coord_arr = bunch->coordArr();

	if(!bunch->hasParticleAttributes("ParticlePhaseAttributes")){
		cerr << "Adding particle phase information attribute.\n";
		std::map<std::string, double> tunemap;
		tunemap.insert(std::make_pair("phase1", 0));
		tunemap.insert(std::make_pair("phase2", 0));
		tunemap.insert(std::make_pair("tune1", 0));
		tunemap.insert(std::make_pair("tune2", 0));
		tunemap.insert(std::make_pair("action1", 0));
		tunemap.insert(std::make_pair("action2", 0));
		bunch->addParticleAttributes("ParticlePhaseAttributes", tunemap);
	}

	if (bunch->hasParticleAttributes("ParticlePhaseAttributes")) {
		for (int i=0; i < bunch->getSize(); i++) {
			// Get phase space coordinates
			double x = part_coord_arr[i][0];
			double xp = part_coord_arr[i][1];
			double y = part_coord_arr[i][2];
			double yp = part_coord_arr[i][3];
			double Etot = syncPart->getEnergy() + syncPart->getMass();
			double dpp = 1.0 / (beta * beta) * part_coord_arr[i][5] / Etot;

			// Normalize coordinates
            double xval  = matrix[0][0] * x + matrix[0][1] * xp + matrix[0][2] * y + matrix[0][3] * yp;
            double xpval = matrix[1][0] * x + matrix[1][1] * xp + matrix[1][2] * y + matrix[1][3] * yp;
            double yval  = matrix[2][0] * x + matrix[2][1] * xp + matrix[2][2] * y + matrix[2][3] * yp;
            double ypval = matrix[3][0] * x + matrix[3][1] * xp + matrix[3][2] * y + matrix[3][3] * yp;

			// Compute phases (mode 1)
			double angle = atan2(xpval, xval);
			if (angle < 0.0) {
				angle += (2.0 * OrbitConst::PI);
			}
			double xPhase = angle;
			double xPhaseOld = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 0);
			double xTune = (xPhaseOld - xPhase) / (2.0*OrbitConst::PI);
			if (xTune < 0.0) {
				xTune += 1.0;
			}
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 0) = xPhase;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 2) = xTune;

			angle = atan2(ypval, yval);
			if (angle < 0.0) {
				angle += (2.0 * OrbitConst::PI);
			}
			double yPhase = angle;
			double yPhaseOld = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 1);
			double yTune = (yPhaseOld - yPhase) / (2.0*OrbitConst::PI);
			if (yTune < 0.0) {
				yTune += 1.0;
			}
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 1) = yPhase;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 3) = yTune;

			double xAction = xval * xval + xpval * xpval;
			double yAction = yval * yval + ypval * ypval;
			
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 4) = xAction;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 5) = yAction;
		}
	}

}
