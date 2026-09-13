#ifndef Dual_Harmonic_Cav_H
#define Dual_Harmonic_Cav_H

//MPI Function Wrappers
#include "mpi/orbit_mpi.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "orbit/Bunch.hh"

//pyORBIT utils


class Dual_Harmonic_Cav
{
  public:
    Dual_Harmonic_Cav(double ZtoPhi, double RFHNum,
                 double RatioRFHNum, double RFVoltage,
                 double RatioVoltage, double RFPhase,
                 double RFPhase2);
    virtual ~Dual_Harmonic_Cav();
    void   setZtoPhi(double ZtoPhi);
    double getZtoPhi();
    void   setRFHNum(double RFHNum);
    double getRatioRFHNum();
    void   setRatioRFHNum(double RatioRFHNum);
    double getRFHNum();
    void   setRFVoltage(double RFVoltage);
    double getRFVoltage();
    void   setRatioVoltage(double RFVoltage);
    double getRatioVoltage();
    void   setRFPhase(double RFPhase);
    double getRFPhase();
    void   setRFPhase2(double RFPhase2);
    double getRFPhase2();
    void   trackBunch(Bunch* bunch);

  private:
  double _ZtoPhi;
  double _RFHNum;
  double _RatioRFHNum;
  double _RFVoltage;
  double _RatioVoltage;
  double _RFPhase;
  double _RFPhase2;

  protected:

};

#endif
