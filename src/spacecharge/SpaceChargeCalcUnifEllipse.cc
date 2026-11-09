/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   SpaceChargeCalcUnifEllipse.cc
//
// AUTHOR
//    A. Shishlo
//
// Created:
//   11/09/10
//
// DESCRIPTION
//  This class calculates the space charge kicks for bunch. It represent the bunch as the set
//  of uniformly charged ellipses in the center of mass of the bunch system.
//  The space charge kick is transformed later into the lab system.
//
/////////////////////////////////////////////////////////////////////////////
#include "SpaceChargeCalcUnifEllipse.hh"
#include "BufferStore.hh"

#include "ParticleMacroSize.hh"

#include <cfloat>
#include <cmath>
#include <iostream>

using namespace OrbitUtils;

SpaceChargeCalcUnifEllipse::SpaceChargeCalcUnifEllipse(int nEllipses_in) : CppPyWrapper(NULL) {
    nEllipses = nEllipses_in;
    ellipsoidCalc_arr = new UniformEllipsoidFieldCalculator *[nEllipses];
    for (int ie = 0; ie < nEllipses; ie++) {
        ellipsoidCalc_arr[ie] = new UniformEllipsoidFieldCalculator();
    }
    macroSizesEll_arr = (double *)malloc(sizeof(double) * nEllipses);
    macroSizesEll_MPI_arr = (double *)malloc(sizeof(double) * nEllipses);
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_arr[ie] = 0.;
        macroSizesEll_MPI_arr[ie] = 0.;
    }
}

SpaceChargeCalcUnifEllipse::~SpaceChargeCalcUnifEllipse() {
    for (int ie = 0; ie < nEllipses; ie++) {
        if (ellipsoidCalc_arr[ie]->getPyWrapper() != NULL) {
            Py_DECREF(ellipsoidCalc_arr[ie]->getPyWrapper());
        } else {
            delete ellipsoidCalc_arr[ie];
        }
    }
    delete[] ellipsoidCalc_arr;

    free(macroSizesEll_arr);
    free(macroSizesEll_MPI_arr);
}

void SpaceChargeCalcUnifEllipse::trackBunch(Bunch *bunch, double length) {

    int nPartsGlobal = bunch->getSizeGlobal();
    if (nPartsGlobal < 3)
        return;

    SyncPart *syncPart = bunch->getSyncPart();
    double beta = syncPart->getBeta();
    double gamma = syncPart->getGamma();

    for (int ie = 0; ie < nEllipses; ie++) {
        ellipsoidCalc_arr[ie]->setQ(0.);
    }

    // analyse the bunch and make the ellipsoid filed sources
    this->bunchAnalysis(bunch);

    // if there is nothing we give up
    if (total_macrosize == 0.)
        return;

    double trans_factor = length * bunch->getClassicalRadius() / (pow(beta, 2) * pow(gamma, 2));
    double long_factor = length * bunch->getClassicalRadius() * bunch->getMass();

    double x, y, z, ex, ey, ez;
    for (int i = 0, n = bunch->getSize(); i < n; i++) {
        x = bunch->x(i) - x_center;
        y = bunch->y(i) - y_center;
        z = (bunch->z(i) - z_center) * gamma;
        this->calculateField(x, y, z, ex, ey, ez);
        // calculate momentum kicks
        bunch->xp(i) += ex * trans_factor;
        bunch->yp(i) += ey * trans_factor;
        bunch->dE(i) += ez * long_factor;
    }
}

/** Analyses the bunch and sets up the ellipsoid filed sources */
void SpaceChargeCalcUnifEllipse::bunchAnalysis(Bunch *bunch) {

    // average values for x,y,z,x2,y2,z2 and total macrosize
    int buff_index0 = 0;
    int buff_index1 = 0;
    double *coord_avg = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0, 7);
    double *coord_avg_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1, 7);
    for (int i = 0; i < 7; i++) {
        coord_avg[i] = 0.;
    }

    // caluclate limits and averages
    double **partArr = bunch->coordArr();
    double *coordArr = NULL;
    bunch->compress();
    double **part_coord_arr = bunch->coordArr();
    int has_msize = bunch->hasParticleAttributes("macrosize");
    if (has_msize > 0) {
        ParticleMacroSize *macroSizeAttr = (ParticleMacroSize *)bunch->getParticleAttributes("macrosize");
        double m_size = 0.;
        for (int ip = 0, n = bunch->getSize(); ip < n; ip++) {
            m_size = macroSizeAttr->macrosize(ip);
            coordArr = partArr[ip];
            coord_avg[0] += m_size * coordArr[0];
            coord_avg[1] += m_size * coordArr[2];
            coord_avg[2] += m_size * coordArr[4];
            coord_avg[3] += m_size * coordArr[0] * coordArr[0];
            coord_avg[4] += m_size * coordArr[2] * coordArr[2];
            coord_avg[5] += m_size * coordArr[4] * coordArr[4];
            coord_avg[6] += m_size;
        }
    } else {
        double m_size = bunch->getMacroSize();
        int nParts = bunch->getSize();
        coord_avg[6] = m_size * nParts;
        for (int ip = 0; ip < nParts; ip++) {
            coordArr = partArr[ip];
            coord_avg[0] += coordArr[0];
            coord_avg[1] += coordArr[2];
            coord_avg[2] += coordArr[4];
            coord_avg[3] += coordArr[0] * coordArr[0];
            coord_avg[4] += coordArr[2] * coordArr[2];
            coord_avg[5] += coordArr[4] * coordArr[4];
        }
        for (int i = 0; i < 6; i++) {
            coord_avg[i] *= m_size;
        }
    }

    // calculates sum over all  CPUs
    ORBIT_MPI_Allreduce(coord_avg, coord_avg_out, 7, MPI_DOUBLE, MPI_SUM, bunch->getMPI_Comm_Local()->comm);

    total_macrosize = coord_avg_out[6];
    if (total_macrosize == 0.) {
        // free resources
        OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
        OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);
        return;
    }

    // calculate the parameters of the biggest ellipse
    x_center = coord_avg_out[0] / total_macrosize;
    y_center = coord_avg_out[1] / total_macrosize;
    z_center = coord_avg_out[2] / total_macrosize;
    x2_avg = fabs(coord_avg_out[3] / total_macrosize - x_center * x_center);
    y2_avg = fabs(coord_avg_out[4] / total_macrosize - y_center * y_center);
    z2_avg = fabs(coord_avg_out[5] / total_macrosize - z_center * z_center);
    a2_ellips = 5.0 * x2_avg;
    b2_ellips = 5.0 * y2_avg;
    c2_ellips = 5.0 * z2_avg;
    a_ellips = sqrt(a2_ellips);
    b_ellips = sqrt(b2_ellips);
    c_ellips = sqrt(c2_ellips);

    // std::cout<<"debug a_ellips="<< a_ellips <<" b_ellips="<< b_ellips <<" c_ellips="<< c_ellips <<std::endl;
    // free resources
    OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
    OrbitUtils::BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);

    // check if the beam size is not zero
    if (x2_avg == 0. || y2_avg == 0. || z2_avg == 0.) {
        int rank = 0;
        ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cerr << "SpaceChargeCalcUnifEllipse::bunchAnalysis(bunch,...)" << std::endl
                      << "The bunch coords min and max sizes are wrong! Cannot calculate space charge!" << std::endl
                      << " x2_rms=" << x2_avg << std::endl
                      << " y2_rms=" << y2_avg << std::endl
                      << " z2_rms=" << z2_avg << std::endl
                      << "Stop." << std::endl;
        }
        ORBIT_MPI_Finalize();
    }

    // relativistic factor gamma
    double gamma = bunch->getSyncPart()->getGamma();

    // if we have only one ellipse we should not distribute anything
    if (nEllipses == 1) {
        macroSizesEll_arr[0] = total_macrosize;
        double r_max = a_ellips;
        if (r_max < b_ellips)
            r_max = b_ellips;
        if (r_max < c_ellips * gamma)
            r_max = c_ellips * gamma;
        ellipsoidCalc_arr[0]->setEllipsoid(a_ellips, b_ellips, c_ellips * gamma, 10. * r_max);
        ellipsoidCalc_arr[0]->setQ(macroSizesEll_arr[0]);
        return;
    }

    // find the distribution of the macrosizes between nEllipses
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_arr[ie] = 0.;
    }

    double pos = 0.;
    int pos_index = 0;
    if (has_msize > 0) {
        ParticleMacroSize *macroSizeAttr = (ParticleMacroSize *)bunch->getParticleAttributes("macrosize");
        double m_size = 0.;
        for (int ip = 0, n = bunch->getSize(); ip < n; ip++) {
            m_size = macroSizeAttr->macrosize(ip);
            coordArr = partArr[ip];
            pos = sqrt(coordArr[0] * coordArr[0] / a2_ellips + coordArr[2] * coordArr[2] / b2_ellips + coordArr[4] * coordArr[4] / c2_ellips);
            pos_index = int(pos * nEllipses);
            if (pos_index < 0)
                pos_index = 0;
            if (pos_index >= nEllipses)
                pos_index = nEllipses - 1;
            macroSizesEll_arr[pos_index] += m_size;
        }
    } else {
        double m_size = bunch->getMacroSize();
        int nParts = bunch->getSize();
        for (int ip = 0, n = bunch->getSize(); ip < n; ip++) {
            coordArr = partArr[ip];
            pos = sqrt(coordArr[0] * coordArr[0] / a2_ellips + coordArr[2] * coordArr[2] / b2_ellips + coordArr[4] * coordArr[4] / c2_ellips);
            pos_index = int(pos * nEllipses) - 1;
            if (pos_index < 0)
                pos_index = 0;
            if (pos_index >= nEllipses)
                pos_index = nEllipses - 1;
            macroSizesEll_arr[pos_index] += m_size;
        }
    }
    // calculates sum over all  CPUs
    ORBIT_MPI_Allreduce(macroSizesEll_arr, macroSizesEll_MPI_arr, nEllipses, MPI_DOUBLE, MPI_SUM, bunch->getMPI_Comm_Local()->comm);
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_arr[ie] = macroSizesEll_MPI_arr[ie];
        // std::cout<<"debug 0 ie ="<< ie <<" macrosize="<< macroSizesEll_MPI_arr[ie] << std::endl;
    }
    // calculate the relative volume density in each region. This density is a sum of all elipsoids
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_MPI_arr[ie] /= ((ie + 2) * (ie + 2) * (ie + 2) - (ie + 1) * (ie + 1) * (ie + 1));
        // std::cout<<"debug 1 ie ="<< ie <<" macrosize="<< macroSizesEll_MPI_arr[ie] << std::endl;
    }
    // calculate the density for each elipsoid
    double rho_sum = 0.;
    for (int ie = (nEllipses - 1); ie >= 0; ie--) {
        macroSizesEll_MPI_arr[ie] -= rho_sum;
        rho_sum += macroSizesEll_MPI_arr[ie];
        // std::cout<<"debug 2 ie ="<< ie <<" macrosize="<< macroSizesEll_MPI_arr[ie] << " rho_sum="<< rho_sum <<std::endl;
    }

    // now set up the relative total charges in ellipsoids
    double q_sum = 0.;
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_MPI_arr[ie] = macroSizesEll_MPI_arr[ie] * (ie + 1) * (ie + 1) * (ie + 1);
        // std::cout<<"debug 3 ie ="<< ie <<" macrosize="<< macroSizesEll_MPI_arr[ie] << std::endl;
        q_sum += macroSizesEll_MPI_arr[ie];
    }
    double q_coeff = total_macrosize / q_sum;
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_arr[ie] = macroSizesEll_MPI_arr[ie] * q_coeff;
    }

    // now we initialize the ellipses filed calculators
    double r_max = a_ellips;
    if (r_max < b_ellips)
        r_max = b_ellips;
    if (r_max < c_ellips * gamma)
        r_max = c_ellips * gamma;
    for (int ie = 0; ie < nEllipses; ie++) {
        double coeff = (ie + 2.) / nEllipses;
        ellipsoidCalc_arr[ie]->setEllipsoid(a_ellips * coeff, b_ellips * coeff, c_ellips * coeff * gamma, 10. * r_max * coeff);
        ellipsoidCalc_arr[ie]->setQ(macroSizesEll_arr[ie]);
        // std::cout<<"debug ie ="<< ie <<" macrosize="<< macroSizesEll_arr[ie]<<" a="
        //          << a_ellips*coeff<< " b="<< b_ellips*coeff<< "  c="<< c_ellips*coeff*gamma
        //				 << " r_max="<< 10.*r_max*coeff << std::endl;
    }
}

/** Calculates the electric filed in the center of the bunch sytem. */
void SpaceChargeCalcUnifEllipse::calculateField(double x, double y, double z,
                                                        double &ex, double &ey, double &ez) {
    ex = 0.;
    ey = 0.;
    ez = 0.;
    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;
    double ex_l, ey_l, ez_l;
    for (int ie = 0; ie < nEllipses; ie++) {
        ellipsoidCalc_arr[ie]->calcField(x, y, z, x2, y2, z2, ex_l, ey_l, ez_l);
        // std::cout<<"debug ie="<<ie<<" ex_l="<<ex_l<<" ey_l="<<ey_l<<" ez_l="<<ez_l<<std::endl;
        ex += ex_l;
        ey += ey_l;
        ez += ez_l;
    }
}

/** Returns the UniformEllipsoidFieldCalculator class instance with a particular index */
UniformEllipsoidFieldCalculator *SpaceChargeCalcUnifEllipse::getEllipsFieldCalculator(int ellipse_index) {
    if (ellipse_index >= 0 && ellipse_index < nEllipses) {
        return ellipsoidCalc_arr[ellipse_index];
    } else {
        return NULL;
    }
}

/** Returns the number of UniformEllipsoidFieldCalculator class instances */
int SpaceChargeCalcUnifEllipse::getNEllipses() {
    return nEllipses;
}

SpaceChargeCalcUnifEllipse2D::SpaceChargeCalcUnifEllipse2D(int nEllipses_in) : SpaceChargeCalcUnifEllipse(nEllipses_in) {
    cos_phi = 1.;
    sin_phi = 0.;
    bunch_length = 0.;
}

void SpaceChargeCalcUnifEllipse2D::trackBunch(Bunch *bunch, double length) {
    int nPartsGlobal = bunch->getSizeGlobal();
    if (nPartsGlobal < 3)
        return;

    for (int ie = 0; ie < nEllipses; ie++) {
        ellipsoidCalc_arr[ie]->setQ(0.);
    }

    this->bunchAnalysis(bunch);
    if (total_macrosize == 0.)
        return;

    SyncPart *syncPart = bunch->getSyncPart();
    double beta = syncPart->getBeta();
    double gamma = syncPart->getGamma();
    double factor = 2. * length * bunch->getClassicalRadius() / (pow(beta, 2) * pow(gamma, 3));

    double ex, ey, ez;
    for (int ip = 0, n = bunch->getSize(); ip < n; ip++) {
        double x = bunch->x(ip) - x_center;
        double y = bunch->y(ip) - y_center;
        this->calculateField(x, y, 0., ex, ey, ez);
        bunch->xp(ip) += ex * factor;
        bunch->yp(ip) += ey * factor;
    }
}

/** Analyses the bunch and sets up the ellipse field sources. */
void SpaceChargeCalcUnifEllipse2D::bunchAnalysis(Bunch *bunch) {
    int buff_index0 = 0;
    int buff_index1 = 0;
    double *coord_avg = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index0, 8);
    double *coord_avg_out = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index1, 8);
    for (int i = 0; i < 8; i++) {
        coord_avg[i] = 0.;
    }

    bunch->compress();
    double **partArr = bunch->coordArr();
    double *coordArr = NULL;
    int has_msize = bunch->hasParticleAttributes("macrosize");
    if (has_msize > 0) {
        ParticleMacroSize *macroSizeAttr = (ParticleMacroSize *)bunch->getParticleAttributes("macrosize");
        for (int ip = 0, n = bunch->getSize(); ip < n; ip++) {
            double m_size = macroSizeAttr->macrosize(ip);
            coordArr = partArr[ip];
            coord_avg[0] += m_size * coordArr[0];
            coord_avg[1] += m_size * coordArr[2];
            coord_avg[2] += m_size * coordArr[4];
            coord_avg[3] += m_size * coordArr[0] * coordArr[0];
            coord_avg[4] += m_size * coordArr[2] * coordArr[2];
            coord_avg[5] += m_size * coordArr[4] * coordArr[4];
            coord_avg[6] += m_size * coordArr[0] * coordArr[2];
            coord_avg[7] += m_size;
        }
    } else {
        double m_size = bunch->getMacroSize();
        int nParts = bunch->getSize();
        coord_avg[7] = m_size * nParts;
        for (int ip = 0; ip < nParts; ip++) {
            coordArr = partArr[ip];
            coord_avg[0] += coordArr[0];
            coord_avg[1] += coordArr[2];
            coord_avg[2] += coordArr[4];
            coord_avg[3] += coordArr[0] * coordArr[0];
            coord_avg[4] += coordArr[2] * coordArr[2];
            coord_avg[5] += coordArr[4] * coordArr[4];
            coord_avg[6] += coordArr[0] * coordArr[2];
        }
        for (int i = 0; i < 7; i++) {
            coord_avg[i] *= m_size;
        }
    }

    ORBIT_MPI_Allreduce(coord_avg, coord_avg_out, 8, MPI_DOUBLE, MPI_SUM, bunch->getMPI_Comm_Local()->comm);

    total_macrosize = coord_avg_out[7];
    if (total_macrosize == 0.) {
        BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
        BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);
        return;
    }

    x_center = coord_avg_out[0] / total_macrosize;
    y_center = coord_avg_out[1] / total_macrosize;
    z_center = coord_avg_out[2] / total_macrosize;
    double cov_xx = coord_avg_out[3] / total_macrosize - x_center * x_center;
    double cov_yy = coord_avg_out[4] / total_macrosize - y_center * y_center;
    double cov_xy = coord_avg_out[6] / total_macrosize - x_center * y_center;
    z2_avg = fabs(coord_avg_out[5] / total_macrosize - z_center * z_center);

    double phi = -0.5 * atan2(2. * cov_xy, cov_xx - cov_yy);
    cos_phi = cos(phi);
    sin_phi = sin(phi);
    double sin_cos_phi = sin_phi * cos_phi;
    x2_avg = fabs(cov_xx * cos_phi * cos_phi + cov_yy * sin_phi * sin_phi - 2. * cov_xy * sin_cos_phi);
    y2_avg = fabs(cov_xx * sin_phi * sin_phi + cov_yy * cos_phi * cos_phi + 2. * cov_xy * sin_cos_phi);
    a2_ellips = 4. * x2_avg;
    b2_ellips = 4. * y2_avg;
    a_ellips = sqrt(a2_ellips);
    b_ellips = sqrt(b2_ellips);
    bunch_length = sqrt(12. * z2_avg);

    BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index0);
    BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index1);

    if (x2_avg == 0. || y2_avg == 0. || z2_avg == 0.) {
        int rank = 0;
        ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cerr << "SpaceChargeCalcUnifEllipse2D::bunchAnalysis(bunch,...)" << std::endl
                      << "The bunch rms sizes are wrong! Cannot calculate space charge!" << std::endl
                      << " x2_rms=" << x2_avg << std::endl
                      << " y2_rms=" << y2_avg << std::endl
                      << " z2_rms=" << z2_avg << std::endl
                      << "Stop." << std::endl;
        }
        ORBIT_MPI_Finalize();
    }

    if (nEllipses == 1) {
        macroSizesEll_arr[0] = total_macrosize;
        ellipsoidCalc_arr[0]->setEllipse(a_ellips, b_ellips);
        ellipsoidCalc_arr[0]->setQ(total_macrosize / bunch_length);
        return;
    }

    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_arr[ie] = 0.;
    }

    for (int ip = 0, n = bunch->getSize(); ip < n; ip++) {
        coordArr = partArr[ip];
        double x = coordArr[0] - x_center;
        double y = coordArr[2] - y_center;
        double x_rot = cos_phi * x - sin_phi * y;
        double y_rot = sin_phi * x + cos_phi * y;
        double pos = sqrt(x_rot * x_rot / a2_ellips + y_rot * y_rot / b2_ellips);
        int pos_index = int(pos * nEllipses) - 1;
        if (pos_index < 0)
            pos_index = 0;
        if (pos_index >= nEllipses)
            pos_index = nEllipses - 1;
        double m_size = bunch->getMacroSize();
        if (has_msize > 0) {
            ParticleMacroSize *macroSizeAttr = (ParticleMacroSize *)bunch->getParticleAttributes("macrosize");
            m_size = macroSizeAttr->macrosize(ip);
        }
        macroSizesEll_arr[pos_index] += m_size;
    }

    ORBIT_MPI_Allreduce(macroSizesEll_arr, macroSizesEll_MPI_arr, nEllipses, MPI_DOUBLE, MPI_SUM, bunch->getMPI_Comm_Local()->comm);
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_arr[ie] = macroSizesEll_MPI_arr[ie];
        macroSizesEll_MPI_arr[ie] /= ((ie + 2) * (ie + 2) - (ie + 1) * (ie + 1));
    }

    double rho_sum = 0.;
    for (int ie = nEllipses - 1; ie >= 0; ie--) {
        macroSizesEll_MPI_arr[ie] -= rho_sum;
        rho_sum += macroSizesEll_MPI_arr[ie];
    }

    double q_sum = 0.;
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_MPI_arr[ie] *= (ie + 1) * (ie + 1);
        q_sum += macroSizesEll_MPI_arr[ie];
    }
    double q_coeff = total_macrosize / q_sum;
    for (int ie = 0; ie < nEllipses; ie++) {
        macroSizesEll_arr[ie] = macroSizesEll_MPI_arr[ie] * q_coeff;
        double coeff = (ie + 2.) / nEllipses;
        ellipsoidCalc_arr[ie]->setEllipse(a_ellips * coeff, b_ellips * coeff);
        ellipsoidCalc_arr[ie]->setQ(macroSizesEll_arr[ie] / bunch_length);
    }
}

/** Calculates the electric field in the bunch-centered system. */
void SpaceChargeCalcUnifEllipse2D::calculateField(double x, double y, double z, double &ex, double &ey, double &ez) {
    double x_rot = cos_phi * x - sin_phi * y;
    double y_rot = sin_phi * x + cos_phi * y;
    double x2 = x_rot * x_rot;
    double y2 = y_rot * y_rot;
    double ex_rot = 0.;
    double ey_rot = 0.;
    for (int ie = 0; ie < nEllipses; ie++) {
        double ex_l, ey_l;
        ellipsoidCalc_arr[ie]->calcField2D(x_rot, y_rot, x2, y2, ex_l, ey_l);
        ex_rot += ex_l;
        ey_rot += ey_l;
    }
    ex = cos_phi * ex_rot + sin_phi * ey_rot;
    ey = -sin_phi * ex_rot + cos_phi * ey_rot;
    ez = 0.;
}
