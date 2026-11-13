#include "BunchTwissAnalysis.hh"

#include <cfloat>
#include <cmath>
#include <iostream>

#include "ParticleMacroSize.hh"
#include "SyncPart.hh"

/** Constructor */
BunchTwissAnalysis::BunchTwissAnalysis() : CppPyWrapper(nullptr)
{

    count = 0;
    _order = 0;

    total_macrosize = 0.;
}

/** Destructor */
BunchTwissAnalysis::~BunchTwissAnalysis() = default;

/** Performs the Twiss analysis of the bunch */
void BunchTwissAnalysis::analyzeBunch(Bunch* bunch)
{

    //initialization
    for (int i = 0; i < 6; i++)
    {
        avg_arr[i] = 0.;
    }

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            corr_arr[i + 6 * j] = 0.;
        }
    }
    count = 0;
    total_macrosize = 0.;

    bunch->compress();
    double m_size = 0.;
    int nParts = bunch->getSize();
    count += nParts;
    double** part_coord_arr = bunch->coordArr();
    int has_msize = bunch->hasParticleAttributes("macrosize");
    if (has_msize > 0)
    {
        ParticleMacroSize* macroSizeAttr =
            (ParticleMacroSize*)bunch->getParticleAttributes("macrosize");
#pragma omp parallel for reduction(+ : total_macrosize) reduction(+ : avg_arr[ : 6])               \
    reduction(+ : corr_arr[ : 36])
        for (int ip = 0; ip < nParts; ip++)
        {
            double const m_size = macroSizeAttr->macrosize(ip);
            total_macrosize += m_size;
            for (int i = 0; i < 6; i++)
            {
                avg_arr[i] += m_size * part_coord_arr[ip][i];
            }

            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < i + 1; j++)
                {
                    corr_arr[i + 6 * j] += m_size * part_coord_arr[ip][i] * part_coord_arr[ip][j];
                }
            }
        }
    }
    else
    {
        m_size = 1.0;
#pragma omp parallel for reduction(+ : avg_arr[ : 6]) reduction(+ : corr_arr[ : 36])
        for (int ip = 0; ip < nParts; ip++)
        {
            for (int i = 0; i < 6; i++)
            {
                avg_arr[i] += part_coord_arr[ip][i];
            }

            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < i + 1; j++)
                {
                    corr_arr[i + 6 * j] += part_coord_arr[ip][i] * part_coord_arr[ip][j];
                }
            }
        }
        total_macrosize += nParts * m_size;
        for (int i = 0; i < 6; i++)
        {
            avg_arr[i] *= m_size;
        }
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < i + 1; j++)
            {
                corr_arr[i + 6 * j] *= m_size;
            }
        }
    }

    ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM,
                        bunch->getMPI_Comm_Local()->comm);
    ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, &total_macrosize, 1, MPI_DOUBLE, MPI_SUM,
                        bunch->getMPI_Comm_Local()->comm);
    ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, avg_arr, 6, MPI_DOUBLE, MPI_SUM,
                        bunch->getMPI_Comm_Local()->comm);
    ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, corr_arr, 36, MPI_DOUBLE, MPI_SUM,
                        bunch->getMPI_Comm_Local()->comm);

    if (std::fabs(total_macrosize) > 0.)
    {
        for (int i = 0; i < 6; i++)
        {
            avg_arr[i] /= total_macrosize;
        }

        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < i + 1; j++)
            {
                corr_arr[i + 6 * j] /= total_macrosize;
                corr_arr[j + 6 * i] = corr_arr[i + 6 * j];
            }
        }
    }

    SyncPart* syncPart = bunch->getSyncPart();

    bunch_momentum = syncPart->getMomentum();
    bunch_beta = syncPart->getBeta();
    bunch_gamma = syncPart->getGamma();
    bunch_kinenergy = syncPart->getEnergy();
    bunch_mass = syncPart->getMass();
}

/** Performs the bunch moments computations */
void BunchTwissAnalysis::computeBunchMoments(Bunch* bunch, int order, int dispersionflag,
                                             int emitnormflag)
{
    _order = order;

    analyzeBunch(bunch);

    const double dispersionCoeff =
        dispersionflag > 0
            ? getDispersion(0) / ((bunch_kinenergy + bunch_mass) * bunch_beta * bunch_beta)
            : 0.0;
    // The corrected centroid follows directly from the weighted coordinate means.
    const double xAvg = getAverage(0) - dispersionCoeff * getAverage(5);
    const double yAvg = getAverage(2);
    const double xbetaterm = std::sqrt(getBeta(0) * (emitnormflag > 0 ? getEmittance(0) : 1.0));
    const double ybetaterm = std::sqrt(getBeta(1) * (emitnormflag > 0 ? getEmittance(1) : 1.0));

    const int stride = _order + 1;
    momentXY.assign(static_cast<std::size_t>(stride) * stride, 0.0);
    if (_order <= 2)
    {
        // analyzeBunch() already reduced all statistics needed for these moments.
        momentXY[0] = 1.;
        if (_order >= 1)
        {
            momentXY[1] = xAvg;
            momentXY[stride] = yAvg;
        }
        if (_order == 2 && std::fabs(total_macrosize) > 0.)
        {
            momentXY[2] =
                (getCorrelation(0, 0) - 2. * dispersionCoeff * getCorrelation(0, 5) +
                 dispersionCoeff * dispersionCoeff * getCorrelation(5, 5)) /
                (xbetaterm * xbetaterm);
            momentXY[stride + 1] =
                (getCorrelation(0, 2) - dispersionCoeff * getCorrelation(2, 5)) /
                (xbetaterm * ybetaterm);
            momentXY[2 * stride] = getCorrelation(2, 2) / (ybetaterm * ybetaterm);
        }
        return;
    }

    const int nParts = bunch->getSize();
    double** part_coord_arr = bunch->coordArr();
    ParticleMacroSize* weights =
        bunch->hasParticleAttributes("macrosize") > 0
            ? static_cast<ParticleMacroSize*>(bunch->getParticleAttributes("macrosize"))
            : nullptr;
    total_macrosize = 0.;
#pragma omp parallel reduction(+ : total_macrosize)
    {
        std::vector<double> momX(stride);
        std::vector<double> momY(stride);
        std::vector<double> localMoments(static_cast<std::size_t>(stride) * stride, 0.0);
        momX[0] = momY[0] = 1.0;

#pragma omp for schedule(static)
        for (int ip = 0; ip < nParts; ip++)
        {
            const double w = weights ? weights->macrosize(ip) : 1.0;
            total_macrosize += w;
            const double dx =
                (part_coord_arr[ip][0] - dispersionCoeff * part_coord_arr[ip][5] - xAvg) /
                xbetaterm;
            const double dy = (part_coord_arr[ip][2] - yAvg) / ybetaterm;

            for (int k = 0; k < _order; k++)
            {
                momX[k + 1] = momX[k] * dx;
                momY[k + 1] = momY[k] * dy;
            }
            for (int j = 0; j <= _order; j++)
            {
                for (int i = 0; i <= _order - j; i++)
                {
                    localMoments[static_cast<std::size_t>(j) * stride + i] += w * momX[i] * momY[j];
                }
            }
        }

// Merge once per thread rather than synchronizing per particle.
#pragma omp critical(BunchTwissAnalysis_moments_merge)
        {
            for (int j = 0; j <= _order; j++)
            {
                for (int i = 0; i <= _order - j; i++)
                {
                    momentXY[static_cast<std::size_t>(j) * stride + i] +=
                        localMoments[static_cast<std::size_t>(j) * stride + i];
                }
            }
        }
    }

    ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, &total_macrosize, 1, MPI_DOUBLE, MPI_SUM,
                        bunch->getMPI_Comm_Local()->comm);

    ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, momentXY.data(), static_cast<int>(momentXY.size()),
                        MPI_DOUBLE, MPI_SUM, bunch->getMPI_Comm_Local()->comm);

    if (std::fabs(total_macrosize) > 0.)
    {
        for (double& moment : momentXY)
        {
            moment /= total_macrosize;
        }
    }

    momentXY[0] = 1.; // 0th moment
    if (_order >= 1)
    {
        momentXY[1] = xAvg;
        momentXY[stride] = yAvg;
    }
}

double BunchTwissAnalysis::getCorrelation(int ic, int jc)
{
    if (ic < 0 || ic > 5 || jc < 0 || jc > 5)
        return 0.;
    return (corr_arr[ic + 6 * jc] - avg_arr[ic] * avg_arr[jc]);
}

/** Returns the XY moment of the beam */
double BunchTwissAnalysis::getBunchMoment(int i, int j)
{
    if (i < 0 || i > _order || j < 0 || j > _order || momentXY.empty())
        return 0.;
    return momentXY[static_cast<std::size_t>(j) * (_order + 1) + i];
}

/** Returns the average value for coordinate with index ic */
double BunchTwissAnalysis::getAverage(int ic)
{
    if (ic < 0 || ic > 5)
        return 0.;
    return avg_arr[ic];
}

/** Returns the total number of analysed macroparticles */
int BunchTwissAnalysis::getGlobalCount()
{
    return count;
}

/** Returns the total macrosize */
double BunchTwissAnalysis::getGlobalMacrosize()
{
    return total_macrosize;
}

/** Returns the emittance for index 0,1,2 - x,y,z planes. */
double BunchTwissAnalysis::getEmittance(int ic)
{
    // for x and y the pure betatron emittance is computed (subtracting the dispersive contribution)
    if (ic < 0 || ic > 2)
        return 0.;
    double x2_avg = std::fabs(this->getCorrelation(2 * ic, 2 * ic));
    double xp2_avg = std::fabs(this->getCorrelation(2 * ic + 1, 2 * ic + 1));
    double x_xp_avg = this->getCorrelation(2 * ic, 2 * ic + 1);
    double x_dE_avg = this->getCorrelation(2 * ic, 5);
    double xp_dE_avg = this->getCorrelation(2 * ic + 1, 5);
    double dE2_avg = std::fabs(this->getCorrelation(5, 5));
    double emitt_rms;
    if (ic == 2 or dE2_avg == 0)
    {
        emitt_rms = std::sqrt(std::fabs(x2_avg * xp2_avg - x_xp_avg * x_xp_avg));
    }
    else
    {
        emitt_rms = std::sqrt(std::fabs((x2_avg - x_dE_avg * x_dE_avg / dE2_avg) *
                                  (xp2_avg - xp_dE_avg * xp_dE_avg / dE2_avg) -
                              (x_xp_avg - x_dE_avg * xp_dE_avg / dE2_avg) *
                                  (x_xp_avg - x_dE_avg * xp_dE_avg / dE2_avg)));
    }
    return emitt_rms;
}

/** Returns the normalized betatron emittance for index 0,1 - x,y planes. */
double BunchTwissAnalysis::getEmittanceNormalized(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    if (ic == 2)
        return this->getEmittance(ic);
    return this->getEmittance(ic) * bunch_gamma * bunch_beta;
}

/** Returns Twiss alpha (without dispersive part for x,y) for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getAlpha(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    double x_xp_avg = this->getCorrelation(2 * ic, 2 * ic + 1);
    double alpha;
    if (ic == 2)
    {
        alpha = -x_xp_avg / this->getEmittance(ic);
    }
    else
    {
        double x_dE_avg = this->getCorrelation(2 * ic, 5);
        double xp_dE_avg = this->getCorrelation(2 * ic + 1, 5);
        double dE2_avg = std::fabs(this->getCorrelation(5, 5));
        alpha = -(x_xp_avg - x_dE_avg * xp_dE_avg / dE2_avg) / this->getEmittance(ic);
    }
    return alpha;
}

/** Returns Twiss beta (without dispersive part for x,y) for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getBeta(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    double x2_avg = std::fabs(this->getCorrelation(2 * ic, 2 * ic));

    double beta;
    if (ic == 2)
    {
        beta = x2_avg / this->getEmittance(ic);
    }
    else
    {
        double x_dE_avg = this->getCorrelation(2 * ic, 5);
        double dE2_avg = std::fabs(this->getCorrelation(5, 5));
        beta = (x2_avg - x_dE_avg * x_dE_avg / dE2_avg) / this->getEmittance(ic);
    }
    return beta;
}

/** Returns Twiss gamma (without dispersive part for x,y) for index 0,1,2 - x,y,z planes.*/
double BunchTwissAnalysis::getGamma(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    double xp2_avg = std::fabs(this->getCorrelation(2 * ic + 1, 2 * ic + 1));
    double gamma;
    if (ic == 2)
    {
        gamma = xp2_avg / this->getEmittance(ic);
    }
    else
    {
        double xp_dE_avg = this->getCorrelation(2 * ic + 1, 5);
        double dE2_avg = std::fabs(this->getCorrelation(5, 5));
        gamma = (xp2_avg - xp_dE_avg * xp_dE_avg / dE2_avg) / this->getEmittance(ic);
    }
    return gamma;
}

/** Returns Twiss dispersion function for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getDispersion(int ic)
{
    if (ic < 0 || ic > 1)
        return 0.;
    double x_dE_avg = this->getCorrelation(2 * ic, 5);
    double dE2_avg = std::fabs(this->getCorrelation(5, 5));
    double dispersion = x_dE_avg / dE2_avg * bunch_momentum * bunch_beta;
    return dispersion;
}

/** Returns Twiss dispersion_prime function for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getDispersionDerivative(int ic)
{
    if (ic < 0 || ic > 1)
        return 0.;
    double xp_dE_avg = this->getCorrelation(2 * ic + 1, 5);
    double dE2_avg = std::fabs(this->getCorrelation(5, 5));
    double dispersion_prime = xp_dE_avg / dE2_avg * bunch_momentum * bunch_beta;
    return dispersion_prime;
}

/** Returns the effective emittance for index 0,1 - x,y planes. */
double BunchTwissAnalysis::getEffectiveEmittance(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    double x_avg = this->getAverage(2 * ic);
    double xp_avg = this->getAverage(2 * ic + 1);
    double x2_avg = std::fabs(this->getCorrelation(2 * ic, 2 * ic));
    double xp2_avg = std::fabs(this->getCorrelation(2 * ic + 1, 2 * ic + 1));
    double x_xp_avg = this->getCorrelation(2 * ic, 2 * ic + 1);
    double emitt_rms = std::sqrt(std::fabs(x2_avg * xp2_avg - x_xp_avg * x_xp_avg));
    return emitt_rms;
}

/** Returns effective Twiss alpha for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getEffectiveAlpha(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    double x_avg = this->getAverage(2 * ic);
    double xp_avg = this->getAverage(2 * ic + 1);
    double x2_avg = std::fabs(this->getCorrelation(2 * ic, 2 * ic));
    double xp2_avg = std::fabs(this->getCorrelation(2 * ic + 1, 2 * ic + 1));
    double x_xp_avg = this->getCorrelation(2 * ic, 2 * ic + 1);
    double emitt2_rms = x2_avg * xp2_avg - x_xp_avg * x_xp_avg;
    if (emitt2_rms <= 0.)
        return 0.;
    double emitt_rms = std::sqrt(emitt2_rms);
    double alpha = -x_xp_avg / emitt_rms;
    return alpha;
}

/** Returns effective Twiss beta for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getEffectiveBeta(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    double x_avg = this->getAverage(2 * ic);
    double xp_avg = this->getAverage(2 * ic + 1);
    double x2_avg = std::fabs(this->getCorrelation(2 * ic, 2 * ic));
    double xp2_avg = std::fabs(this->getCorrelation(2 * ic + 1, 2 * ic + 1));
    double x_xp_avg = this->getCorrelation(2 * ic, 2 * ic + 1);
    double emitt2_rms = x2_avg * xp2_avg - x_xp_avg * x_xp_avg;
    if (emitt2_rms <= 0.)
        return 0.;
    double emitt_rms = std::sqrt(emitt2_rms);
    double beta = x2_avg / emitt_rms;
    return beta;
}

/** Returns effective Twiss gamma for index 0,1 - x,y planes.*/
double BunchTwissAnalysis::getEffectiveGamma(int ic)
{
    if (ic < 0 || ic > 2)
        return 0.;
    double x_avg = this->getAverage(2 * ic);
    double xp_avg = this->getAverage(2 * ic + 1);
    double x2_avg = std::fabs(this->getCorrelation(2 * ic, 2 * ic));
    double xp2_avg = std::fabs(this->getCorrelation(2 * ic + 1, 2 * ic + 1));
    double x_xp_avg = this->getCorrelation(2 * ic, 2 * ic + 1);
    double emitt2_rms = x2_avg * xp2_avg - x_xp_avg * x_xp_avg;
    if (emitt2_rms <= 0.)
        return DBL_MAX;
    double emitt_rms = std::sqrt(emitt2_rms);
    double gamma = xp2_avg / emitt_rms;
    return gamma;
}
