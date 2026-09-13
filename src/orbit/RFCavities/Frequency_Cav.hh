#ifndef FREQUENCY_CAV_H
#define FREQUENCY_CAV_H

//MPI Function Wrappers
#include "mpi/orbit_mpi.hh"

#include <cstdlib>
#include <cmath>

//ORBIT bunch
#include "orbit/Bunch.hh"

//pyORBIT utils


class Frequency_Cav
{
  public:
    Frequency_Cav(double RFFreq, double RFE0TL, double RFPhase);
    virtual ~Frequency_Cav();
    void   setRFFreq(double RFFreq);
    double getRFFreq();
    void   setRFE0TL(double RFE0TL);
    double getRFE0TL();
    void   setRFPhase(double RFPhase);
    double getRFPhase();
    void   trackBunch(Bunch* bunch);

  private:
    double _RFFreq;
    double _RFE0TL;
    double _RFPhase;

  protected:

};

#endif
