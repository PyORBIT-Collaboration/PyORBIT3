#ifndef BUNCH_TWISS_ANALYSIS_H
#define BUNCH_TWISS_ANALYSIS_H

#include "utils/CppPyWrapper.hh"

#include <array>
#include <vector>

class Bunch;
class ParticleMacroSize;

/**
 * Calculates the average of 6D coordinates and their correlations.
 * Returns the Twiss parameters for each plane.
 */
class BunchTwissAnalysis : public OrbitUtils::CppPyWrapper
{
public:
  BunchTwissAnalysis() = default;
  virtual ~BunchTwissAnalysis() = default;

  BunchTwissAnalysis(const BunchTwissAnalysis&) = delete;
  BunchTwissAnalysis& operator=(const BunchTwissAnalysis&) = delete;
  BunchTwissAnalysis(BunchTwissAnalysis&&) = delete;
  BunchTwissAnalysis& operator=(BunchTwissAnalysis&&) = delete;

  /**
   * @brief Performs a default Twiss analysis of the bunch.
   *
   * Equivalent to calling computeBunchMoments() with default parameters (order=2, no
   * normalization, no dispersion correction).
   *
   * After this call the Twiss parameters (emittance, alpha, beta, gamma, dispersion) are
   * available from the corresponding getters.
   *
   * @param bunch The bunch to analyze (compressed before use).
   */
  void analyzeBunch(Bunch* bunch);

  /**
   * @brief Computes central moments of the bunch up to the given order.
   *
   * This is the full-featured entry point. For a simple default analysis, prefer
   * analyzeBunch() instead.
   *
   * The 6D coordinate system is \f$(x, x', y, y', z, \delta E)\f$. Higher-order XY
   * moments \f$\langle \tilde{x}^i \tilde{y}^j \rangle\f$ with
   * \f$i+j \le \mathrm{order}\f$ are available via getBunchMoment(), where
   * \f$\tilde{x}\f$ and \f$\tilde{y}\f$ are the (optionally normalized) dispersive
   * and transverse displacements.
   *
   * @param bunch          The bunch to analyze (compressed before use).
   * @param order          Maximum moment order (default 2).
   * @param normalize      Whether to normalize the moments to \f$\sqrt{\beta_{x,y}}\f$.
   * @param emitnormflag   If true, moments are normalized to \f$\sqrt{\beta_{x,y}
   * \varepsilon_{x,y}}\f$. Implies normalize=true (default false).
   * @param dispersionflag If true, x-coordinates are corrected for dispersion
   *                       (default false).
   */
  void computeBunchMoments(
    Bunch* bunch,
    int order = 2,
    bool normalize = false,
    bool emitnormflag = false,
    bool dispersionflag = false
  );

  /**
   * @brief Returns the centered covariance for a pair of coordinates.
   *
   * Computes \f$\langle (u - \langle u \rangle)(v - \langle v \rangle) \rangle\f$ where
   * \f$u\f$ is coordinate \p i and \f$v\f$ is coordinate \p j.
   *
   * @param i Coordinate index (0-5): 0=x, 1=x', 2=y, 3=y', 4=z, 5=dE.
   * @param j Coordinate index (0-5): 0=x, 1=x', 2=y, 3=y', 4=z, 5=dE.
   * @return The covariance, or 0.0 if \p i or \p j is out of range.
   **/
  double getCovariance(std::size_t i, std::size_t j) const;

  /**
   * @brief Returns the Pearson correlation coefficient for a pair of coordinates.
   *
   * Computes \f$\frac{\langle (u - \langle u \rangle)(v - \langle v \rangle) \rangle}
   * {\sigma_u \sigma_v}\f$ where \f$u\f$ is coordinate \p i and \f$v\f$ is coordinate \p j.
   *
   * @param i Coordinate index (0-5): 0=x, 1=x', 2=y, 3=y', 4=z, 5=dE.
   * @param j Coordinate index (0-5): 0=x, 1=x', 2=y, 3=y', 4=z, 5=dE.
   * @return The correlation in [-1, 1], or 0.0 if \p i or \p j is out of range, or if
   *         either variance is non-positive.
   **/
  double getCorrelation(std::size_t i, std::size_t j) const;

  /** Returns the average value for coordinate with index ic. */
  double getAverage(int ic) const;

  /** Returns the total number of analysed macroparticles. */
  std::size_t getGlobalCount() const;

  /** Returns the total macrosize. */
  double getGlobalMacrosize() const;

  /** Returns the XY moment of the beam. */
  double getBunchMoment(int i, int j) const;

  /** Returns the pure betatron emittance for index 0,1,2 - x,y,z planes. */
  double getEmittance(int ic) const;

  /** Returns the normalized pure betatron emittance for index 0,1 - x,y planes. */
  double getEmittanceNormalized(int ic) const;

  /** Returns Twiss alpha (without dispersive part for x,y) for index 0,1,2 - x,y,z planes. */
  double getAlpha(int ic) const;

  /** Returns Twiss beta (without dispersive part for x,y) for index 0,1,2 - x,y,z planes. */
  double getBeta(int ic) const;

  /** Returns Twiss gamma (without dispersive part for x,y) for index 0,1,2 - x,y,z planes. */
  double getGamma(int ic) const;

  /** Returns Twiss dispersion function for index 0,1 - x,y planes. */
  double getDispersion(int ic) const;

  /** Returns Twiss dispersion_prime function for index 0,1 - x,y planes. */
  double getDispersionDerivative(int ic) const;

  /** Returns the effective emittance for index 0,1 - x,y planes. */
  double getEffectiveEmittance(int ic) const;

  /** Returns effective Twiss alpha for index 0,1 - x,y planes. */
  double getEffectiveAlpha(int ic) const;

  /** Returns effective Twiss beta for index 0,1 - x,y planes. */
  double getEffectiveBeta(int ic) const;

  /** Returns effective Twiss gamma for index 0,1 - x,y planes. */
  double getEffectiveGamma(int ic) const;

private:
  static constexpr int N = 6;      // {x, x', y, y', z, dE}
  static constexpr int NN = N * N; // 36

  using Array2D = std::array<std::array<double, N>, N>;

  template <bool HasMacrosizeAttr, bool IncludeXDispersion>
  void computeBunchMomentsImpl(Bunch* bunch, bool normalize, bool emitnormflag);

  std::size_t count_{};
  double total_macrosize_{};

  std::array<double, N> averages_{};
  Array2D moments_{};

  std::vector<double> momentsXY_;

  int order_{};

  double bunch_momentum_{};
  double bunch_gamma_{};
  double bunch_beta_{};
  double bunch_kinenergy_{};
  double bunch_mass_{};
};

#endif
