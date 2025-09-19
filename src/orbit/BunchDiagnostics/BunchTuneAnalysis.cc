#include "BunchTuneAnalysis.hh"
#include "SyncPart.hh"
#include "OrbitConst.hh"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>


BunchTuneAnalysis::BunchTuneAnalysis(): CppPyWrapper(NULL) {
    double matrix[6][6] = {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
    };
}


BunchTuneAnalysis::~BunchTuneAnalysis() {}

void BunchTuneAnalysis::setNormMatrixElement(int i, int j, double value) {
    matrix[i][j] = value;
}

double BunchTuneAnalysis::getNormMatrixElement(int i, int j) {
    return matrix[i][j];
}

void BunchTuneAnalysis::assignTwiss(double betax, double alphax, double etax, double etapx, double betay, double alphay) {
	// Set V_{-1} = I
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			matrix[i][j] = 0.0;
		}
	}
	for (int i = 0; i < 6; i++) {
		matrix[i][i] = 1.0;
	}

	// 2D normalization (x-x')
	matrix[0][0] = 1.0 / sqrt(betax);
	matrix[0][1] = 0.0;
	matrix[1][0] = alphax / sqrt(betax);
	matrix[1][1] = sqrt(betax);

	// 2D normalization (y-y')
	matrix[2][2] = 1.0 / sqrt(betay);
	matrix[2][3] = 0.0;
	matrix[3][2] = alphay / sqrt(betay);
	matrix[3][3] = sqrt(betay);

	// Dispersion (x-x')
	matrix[0][5] = -etax / sqrt(betax);
	matrix[1][5] = -etax * (alphax / sqrt(betax)) - etapx * sqrt(betax);

	// Dispersion (y-y')
	double etay = 0.0;
	double etapy = 0.0;
	matrix[2][5] = -etay / sqrt(betay);
	matrix[3][5] = -etay * (alphay / sqrt(betay)) - etapy * sqrt(betay);

}

void BunchTuneAnalysis::analyzeBunch(Bunch* bunch){

	bunch->compress();
	SyncPart* syncPart = bunch->getSyncPart();
	double beta = syncPart->getBeta();
	double** part_coord_arr = bunch->coordArr();

	if(!bunch->hasParticleAttributes("ParticlePhaseAttributes")){
		cerr<<"adding particle phase information attribute\n";
		std::map<std::string, double> tunemap;
		tunemap.insert(std::make_pair("phase_1", 0));
		tunemap.insert(std::make_pair("phase_2", 0));
		tunemap.insert(std::make_pair("tune_1", 0));
		tunemap.insert(std::make_pair("tune_2", 0));
		tunemap.insert(std::make_pair("action_1", 0));
		tunemap.insert(std::make_pair("action_2", 0));
		bunch->addParticleAttributes("ParticlePhaseAttributes", tunemap);
	}

	if (bunch->hasParticleAttributes("ParticlePhaseAttributes")){
		for (int i=0; i < bunch->getSize(); i++)
		{
			// Extract phase space coordinates
			double x = part_coord_arr[i][0];
			double xp = part_coord_arr[i][1];
			double y = part_coord_arr[i][2];
			double yp = part_coord_arr[i][3];
			double z = part_coord_arr[i][4];
			double Etot = syncPart->getEnergy() + syncPart->getMass();
			double dpp = 1.0 / (beta * beta) * part_coord_arr[i][5] / Etot;

			// Normalize phase space coordinates
            double xval  = matrix[0][0] * x + matrix[0][1] * xp + matrix[0][2] * y + matrix[0][3] * yp + matrix[0][4] * z + matrix[0][5] * dpp;
            double xpval = matrix[1][0] * x + matrix[1][1] * xp + matrix[1][2] * y + matrix[1][3] * yp + matrix[1][4] * z + matrix[1][5] * dpp;
            double yval  = matrix[2][0] * x + matrix[2][1] * xp + matrix[2][2] * y + matrix[2][3] * yp + matrix[2][4] * z + matrix[2][5] * dpp;
            double ypval = matrix[3][0] * x + matrix[3][1] * xp + matrix[3][2] * y + matrix[3][3] * yp + matrix[3][4] * z + matrix[3][5] * dpp;

			// Compute phase advance in normalized x-x'
			double angle = atan2(xpval, xval);
			if (angle < 0.0) {
				angle += (2.0 * OrbitConst::PI);
			}
			double xPhase = angle;
			double xPhaseOld = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 0);
			double xTune = (xPhaseOld - xPhase) / (2.0 * OrbitConst::PI);
			if (xTune < 0.0) { 
				xTune += 1.0;
			}
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 0) = xPhase;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 2) = xTune;

			// Compute phase advance in normalized y-y'
			angle = atan2(ypval, yval);
			if (angle < 0.0) {
				angle += (2.0 * OrbitConst::PI);
			}
			double yPhase = angle;
			double yPhaseOld = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 1);
			double yTune = (yPhaseOld - yPhase) / (2.0 * OrbitConst::PI);
			if (yTune < 0.0) {
				yTune += 1.0;
			}
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 1) = yPhase;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 3) = yTune;

			// Compute actions
			double xAction = (xval * xval + xpval * xpval) / 2.0;
			double yAction = (yval * yval + ypval * ypval) / 2.0;

			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 4) = xAction;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 5) = yAction;
		}
	}

}
