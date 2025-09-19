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
			double x  = part_coord_arr[i][0];
			double xp = part_coord_arr[i][1];
			double y  = part_coord_arr[i][2];
			double yp = part_coord_arr[i][3];
			double z  = part_coord_arr[i][4];
			double Etot = syncPart->getEnergy() + syncPart->getMass();
			double dpp = 1.0 / (beta * beta) * part_coord_arr[i][5] / Etot;

			// Normalize phase space coordinates
            double u1  = matrix[0][0] * x + matrix[0][1] * xp + matrix[0][2] * y + matrix[0][3] * yp + matrix[0][4] * z + matrix[0][5] * dpp;
            double u1p = matrix[1][0] * x + matrix[1][1] * xp + matrix[1][2] * y + matrix[1][3] * yp + matrix[1][4] * z + matrix[1][5] * dpp;
            double u2  = matrix[2][0] * x + matrix[2][1] * xp + matrix[2][2] * y + matrix[2][3] * yp + matrix[2][4] * z + matrix[2][5] * dpp;
            double u2p = matrix[3][0] * x + matrix[3][1] * xp + matrix[3][2] * y + matrix[3][3] * yp + matrix[3][4] * z + matrix[3][5] * dpp;
            double u3  = matrix[4][0] * x + matrix[4][1] * xp + matrix[4][2] * y + matrix[4][3] * yp + matrix[4][4] * z + matrix[4][5] * dpp;
            double u3p = matrix[5][0] * x + matrix[5][1] * xp + matrix[5][2] * y + matrix[5][3] * yp + matrix[5][4] * z + matrix[5][5] * dpp;
			
			// Compute phase advance (mode 1)
			double angle = atan2(u1p, u1);
			if (angle < 0.0) {
				angle += (2.0 * OrbitConst::PI);
			}
			double phase1 = angle;
			double phase1Old = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 0);
			double tune1 = (phase1Old - phase1) / (2.0 * OrbitConst::PI);
			if (tune1 < 0.0) { 
				tune1 += 1.0;
			}

			// Compute phase advance (mode 2)
			angle = atan2(u2p, u2);
			if (angle < 0.0) {
				angle += (2.0 * OrbitConst::PI);
			}
			double phase2 = angle;
			double phase2Old = bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 1);
			double tune2 = (phase2Old - phase2) / (2.0 * OrbitConst::PI);
			if (tune2 < 0.0) {
				tune2 += 1.0;
			}

			// Compute actions
			double action1 = (u1 * u1 + u1p * u1p) / 2.0;
			double action2 = (u2 * u2 + u2p * u2p) / 2.0;

			// Update ParticlePhaseAttributes
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 0) = phase1;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 1) = phase2;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 2) = tune1;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 3) = tune2;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 4) = action1;
			bunch->getParticleAttributes("ParticlePhaseAttributes")->attValue(i, 5) = action2;
		}
	}

}