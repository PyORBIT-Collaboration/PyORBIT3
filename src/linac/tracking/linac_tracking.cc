/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//   linac_tracking.cc
//
// AUTHOR
//   Andrei Shishlo, ORNL, shishlo@ornl.gov
//   2017.01.24
//
// DESCRIPTION
//   Define elementary functions for tracking through
//   specific linac accelerator elements without using
//   TEAPOT algorithms. It is necessary for the cases
//   with a huge energy spread in the bunch.
//
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Include files
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Local Functions:
//
///////////////////////////////////////////////////////////////////////////

#include "OrbitConst.hh"
#include "Bunch.hh"
#include "SyncPart.hh"

#include <complex>

namespace linac_tracking
{

    ///////////////////////////////////////////////////////////////////////////
    // NAME
    //   linac_drift
    //
    // DESCRIPTION
    //   Drifts a particle bunch. Length < 0 is not allowed.
    //
    // PARAMETERS
    //   bunch = reference to the macro-particle bunch
    //   length = length of the drift
    //
    // RETURNS
    //   Nothing
    //
    ///////////////////////////////////////////////////////////////////////////

    void linac_drift(Bunch* bunch, double length)
    {
        if(length <= 0.) return;

        SyncPart* syncPart = bunch->getSyncPart();

        const double beta_s = syncPart->getBeta();
        const double delta_t = length / (OrbitConst::c * beta_s);

        syncPart->setTime(syncPart->getTime() + delta_t);

        const double mass = syncPart->getMass();
        const double Ekin_s = syncPart->getEnergy();
        const double p_s = syncPart->getMomentum();

        //coordinate array [part. index][x,xp,y,yp,z,dE]
        //sqrt(fabs(p_z2)) - fabs function is a protection for case when xp and yp are too big for dE
        //It means the particles are nonphysical and simulations do not make sense
        double** arr = bunch->coordArr();
        const int bunch_size = bunch->getSize();

        #pragma omp parallel for schedule(static)
        for(int i = 0; i < bunch_size; ++i)
        {
            const double dE = arr[i][5];
            const double px = arr[i][1]*p_s;
            const double py = arr[i][3]*p_s;
            const double Ekin = Ekin_s + dE;
            const double Etotal = Ekin_s + dE + mass;
            const double p2 = Ekin*(Etotal + mass);

            const double length_over_pz = length / std::sqrt(std::fabs(p2 - px*px - py*py));

            arr[i][0] += length_over_pz * px;
            arr[i][2] += length_over_pz * py;
            arr[i][4] += length - length_over_pz * beta_s * Etotal;
        }
    }

    ////////////////////////////
    // NAME
    //   linac_quad1
    //
    // DESCRIPTION
    //   Quadrupole element one: linear transport matrix for linac
    //
    // PARAMETERS
    //   bunch  = reference to the macro-particle bunch
    //   length = length of transport
    //   kq = quadrupole field strength [m^(-2)]
    //   kq = G/(B*rho) B*rho = 3.33564*P0/charge [T*m] where P0 in GeV/c
    //   G = dB/dr [T/m] quad gradient
    //
    // RETURNS
    //   Nothing
    //
    ///////////////////////////////////////////////////////////////////////////

    void linac_quad1(Bunch* bunch, double length, double kq, int /*useCharge*/)
    {
        if(kq == 0. || bunch->getCharge() == 0.)
        {
            linac_drift(bunch, length);
            return;
        }

        SyncPart* syncPart = bunch->getSyncPart();
        const double beta_s = syncPart->getBeta();
        if(length > 0.)
        {
            const double delta_t = length / (OrbitConst::c * beta_s);
            syncPart->setTime(syncPart->getTime() + delta_t);
        }

        const double mass = syncPart->getMass();
        const double Ekin_s = syncPart->getEnergy();
        const double p_s = syncPart->getMomentum();
        const bool focus_x = kq > 0.;

        // Coordinate array [particle index][x,xp,y,yp,z,dE].
        double** arr = bunch->coordArr();
        const int nParts = bunch->getSize();

        #pragma omp parallel for schedule(static)
        for(int i = 0; i < nParts; ++i)
        {
            double* part = arr[i];
            const double xp = part[1], yp = part[3];
            const double Ekin = Ekin_s + part[5];
            const double Etotal = Ekin + mass;
            const double p = std::sqrt(Ekin * (Ekin + 2.0 * mass));
            // Rigidity constants and charge cancel: kqc = kq * p_s / p.
            const double coeff = p_s / p;

            // Same time-of-flight expression and entrance momenta as drift.
            // fabs protects nonphysical slopes, matching linac_drift.
            const double px = xp * p_s, py = yp * p_s;
            const double length_over_pz = length / std::sqrt(std::fabs(
                Ekin * (Etotal + mass) - px*px - py*py));
            const double sqrt_kq = std::sqrt(std::fabs(kq * coeff));
            const double phase = sqrt_kq * length;

            double cx, sx, cy, sy, m21, m43;
            if(focus_x)
            {
                cx = std::cos(phase);
                sx = std::sin(phase);
                cy = std::cosh(phase);
                sy = std::sinh(phase);
                m21 = -sx * sqrt_kq;
                m43 = sy * sqrt_kq;
            }
            else
            {
                cx = std::cosh(phase);
                sx = std::sinh(phase);
                cy = std::cos(phase);
                sy = std::sin(phase);
                m21 = sx * sqrt_kq;
                m43 = -sy * sqrt_kq;
            }
            const double m12 = sx / sqrt_kq;
            const double m34 = sy / sqrt_kq;
            const double xp_scaled = xp * coeff;
            const double yp_scaled = yp * coeff;

            const double x = part[0], y = part[2];
            part[0] = x * cx + xp_scaled * m12;
            part[1] = (x * m21 + xp_scaled * cx) / coeff;
            part[2] = y * cy + yp_scaled * m34;
            part[3] = (y * m43 + yp_scaled * cy) / coeff;
            part[4] += length - length_over_pz * beta_s * Etotal;
        }
    }

  ///////////////////////////////////////////////////////////////////////////
  // NAME
  //  linac_ quad2
  //
  // DESCRIPTION
  //   Quadrupole element two: nonlinear piece for linac
  //
  // PARAMETERS
  //   bunch  = reference to the macro-particle bunch
  //   length = length of the element
  //
  // RETURNS
  //   Nothing
  //
  ///////////////////////////////////////////////////////////////////////////

  void linac_quad2(Bunch* /*bunch*/, double /*length*/)
  {
      //there are no non-linear components
  }


    ////////////////////////////
    // NAME
    //   linac_quad3
    //
    // DESCRIPTION
    //   Quadrupole element 3: longitudinal field component of the quad field
    //   Bz(x,y,z) = x*y*dG(z)/dz
    //   This function performs the transverse kicks correction, so the length
    //   is just a parameter.
    //
    // PARAMETERS
    //   bunch  = reference to the macro-particle bunch
    //   length = length of transport
    //   G = dB/dr [T/m] quad gradient
    //   0.299792458 = (1/v_light) * 10^9
    //
    // RETURNS
    //   Nothing
    //
    ///////////////////////////////////////////////////////////////////////////

    void linac_quad3(Bunch* bunch, double length, double dB_dz)
    {
    double charge = bunch->getCharge();

    if(dB_dz == 0. || charge == 0.)
    {
        return;
    }

    SyncPart* syncPart = bunch->getSyncPart();
    //momentum in GeV/c
    double momentum = syncPart->getMomentum();
    double kick_coeff = 0.299792458*charge*dB_dz*length/momentum;

    //coordinate array [part. index][x,xp,y,yp,z,dE]
    double** arr = bunch->coordArr();
    int nParts = bunch->getSize();

    double coef_xy = 0.;

    for(int i = 0; i < nParts; i++){
        coef_xy = kick_coeff*arr[i][0]*arr[i][2];
        arr[i][1]  +=  arr[i][3]*coef_xy;
        arr[i][3]  += -arr[i][1]*coef_xy;
    }
  }

///////////////////////////////////////////////////////////////////////////
// NAME
//   kick
//
// DESCRIPTION
//   Kicks a particle bunch
//
// PARAMETERS
//   bunch = reference to the macro-particle bunch
//   kx = strength of the horizontal kick in rad
//   ky = strength of the vertical kick in rad
//   kE = strength of the energy kick in GeV
//
// RETURNS
//   Nothing
//
///////////////////////////////////////////////////////////////////////////

    void kick(Bunch* bunch, double kx, double ky, double kE, int /*useCharge*/)
    {
        SyncPart* syncPart = bunch->getSyncPart();

        double mass = syncPart->getMass();
        double Ekin_s = syncPart->getEnergy();
        double p_s = syncPart->getMomentum();
        double p2_s = p_s*p_s;
        double dE = 0.;
        double p2 = 0.;
        double p_z = 0.;
        double Ekin = 0.;
        double coeff = 0.;

        //coordinate array [part. index][x,xp,y,yp,z,dE]
        double** arr = bunch->coordArr();
        if(kx != 0.)
        {
            for(int i = 0; i < bunch->getSize(); i++)
            {
                dE = arr[i][5];
                Ekin = Ekin_s+dE;
                p2 = Ekin*(Ekin+2.0*mass);
                p_z = sqrt(p2 - (arr[i][1]*arr[i][1] + arr[i][3]*arr[i][3])*p2_s);
                coeff = p_s/p_z;
                arr[i][1] += kx*coeff;
            }
        }
        if(ky != 0.)
        {
            for(int i = 0; i < bunch->getSize(); i++)
            {
                dE = arr[i][5];
                Ekin = Ekin_s+dE;
                p2 = Ekin*(Ekin+2.0*mass);
                p_z = sqrt(p2 - (arr[i][1]*arr[i][1] + arr[i][3]*arr[i][3])*p2_s);
                coeff = p_s/p_z;
                arr[i][3] += ky*coeff;
            }
        }
        if(kE != 0.)
        {
            for(int i = 0; i < bunch->getSize(); i++)
            {
                arr[i][5] += kE;
            }
        }
    }

}  //end of namespace linac_tracking
