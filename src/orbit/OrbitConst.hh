//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    OrbitConst.hh
//
// AUTHOR
//    A. Shishlo
//
// CREATED
//    06/28/2005
//
// DESCRIPTION
//    Keeps all constants for ORBIT
//
///////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ORBIT_CONSTANTS_H
#define ORBIT_CONSTANTS_H

namespace OrbitConst
{

//-------------------------------------------------------------
//constants
//-------------------------------------------------------------

// PI
constexpr double PI = 3.141592653589793238462643383279502;

//speed of light in [m/sec]
constexpr double c = 2.99792458e+8;

//elementary charge (CGS system)
constexpr double elementary_charge_CGS = 4.8032068e-10;

//elementary charge (MKS system)
constexpr double elementary_charge_MKS = 1.60217565e-19;

// GeV/c to Tm
constexpr double GeV_per_c_to_T_m = 1.0e9 / c;

//coeffitient for magnetic field strength
//to shift from Tesla to H(CGS)/e(SGS) in 1/m^2
constexpr double coeff_Tesla_to_inner = (1.0e+4 / 4.8032068e-10) * 0.01 * (1.0e+6);

//to shift from Volts/m to E in 1/m^2
constexpr double coeff_VoltsPerM_to_inner = 1.0 / ((1.6021773e-19) * (8.98755179e+9));

//coefficient to get volts from our potentials
constexpr double coeff_Phi_to_Volts = (1.6021773e-19) * (8.98755179e+9);

//magnetic permeability
constexpr double permeability = 4 * 3.14159265358979323846264 * 1.0e-007;

//------------------------------------------------------------
//electron parameters
//------------------------------------------------------------

//electron mass in GeV
constexpr double mass_electron = 0.510998910e-3;

//electron classical radius in m
constexpr double classicalRadius_electron = 2.8179409e-15;

//charge of the electron
constexpr double charge_electron = -1.0;

//------------------------------------------------------------
//proton parameters
//------------------------------------------------------------

//proton mass in GeV
constexpr double mass_proton = 0.938272013;

//proton classical radius in m
constexpr double classicalRadius_proton = 1.534698e-18;

//charge of the proton
constexpr double charge_proton = +1.0;

constexpr double tiny = 2.22044604925031e-16;

} // namespace OrbitConst

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
