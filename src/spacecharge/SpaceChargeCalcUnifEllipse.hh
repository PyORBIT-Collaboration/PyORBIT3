/**
 This class calculates the space charge kicks for bunch. It represent the bunch as the set
 of uniformly charged ellipses in the center of mass of the bunch system.
 The space charge kick is transformed later into the lab system.
*/

#ifndef SPACE_CHARGE_CALC_UNIFORM_ELLIPSE_HH
#define SPACE_CHARGE_CALC_UNIFORM_ELLIPSE_HH

// MPI Function Wrappers
#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <cmath>
#include <cstdlib>

// ORBIT bunch
#include "Bunch.hh"

// pyORBIT utils
#include "CppPyWrapper.hh"

#include "UniformEllipsoidFieldCalculator.hh"

using namespace std;

class SpaceChargeCalcUnifEllipse : public OrbitUtils::CppPyWrapper {
  public:
    /** Constructor with the "x to y ratio" parameter. */
    SpaceChargeCalcUnifEllipse(int nEllipses_in);

    /** Destructor. */
    virtual ~SpaceChargeCalcUnifEllipse();

    /** Calculates space charge and applies the transverse and longitudinal SC kicks to the macro-particles in the bunch. */
    virtual void trackBunch(Bunch *bunch, double length);

    /** Analyzes the bunch and sets up the ellipsoid filed sources. */
    virtual void bunchAnalysis(Bunch *bunch);

    /** Calculates the electric filed in the center of the bunch system. */
    virtual void calculateField(double x, double y, double z, double &ex, double &ey, double &ez);

    /** Returns the UniformEllipsoidFieldCalculator class instance with a particular index. */
    UniformEllipsoidFieldCalculator *getEllipsFieldCalculator(int ellipse_index);

    /** Returns the number of UniformEllipsoidFieldCalculator class instances. */
    int getNEllipses();

  private:
  protected:
    // Number of ellipsoids
    int nEllipses;

    // Total macrosize
    double total_macrosize;

    // Distribution parameters
    double x_center;
    double y_center;
    double z_center;
    double x2_avg;
    double y2_avg;
    double z2_avg;
    double xMin;
    double xMax;
    double yMin;
    double yMax;
    double zMin;
    double zMax;

    // Sizes of the biggest ellipsoid
    double a_ellips;
    double b_ellips;
    double c_ellips;
    double a2_ellips;
    double b2_ellips;
    double c2_ellips;

    // Field calculators
    UniformEllipsoidFieldCalculator **ellipsoidCalc_arr;

    // Total macrosize in each ellipsoid
    double *macroSizesEll_arr;
    double *macroSizesEll_MPI_arr;
};

class SpaceChargeCalcUnifEllipse2D : public SpaceChargeCalcUnifEllipse {
  public:
    /** Constructor. */
    SpaceChargeCalcUnifEllipse2D(int nEllipses_in);

    /** Calculates and applies transverse space-charge kicks. */
    void trackBunch(Bunch *bunch, double length) override;

    /** Analyses the bunch and sets up the ellipse field sources. */
    void bunchAnalysis(Bunch *bunch) override;

    /** Calculates the electric field in the bunch-centered system. */
    void calculateField(double x, double y, double z, double &ex, double &ey, double &ez) override;

  protected:
    double cos_phi;
    double sin_phi;
    double bunch_length;
};
// end of SPACE_CHARGE_CALC_UNIFORM_ELLIPSE_HH

#endif
