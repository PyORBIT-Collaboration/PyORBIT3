#include "orbit/BunchDiagnostics/BunchTwissAnalysis.hh"

#include "orbit/Bunch.hh"
#include "orbit/ParticlesAttributes/ParticleMacroSize.hh"
#include "orbit/SyncPart.hh"

#include <atomic>
#include <cmath>
#include <limits>

/** Performs the Twiss analysis of the bunch */
void BunchTwissAnalysis::analyzeBunch(Bunch* bunch)
{
  static std::atomic_flag warned = ATOMIC_FLAG_INIT;
  if (!warned.test_and_set(std::memory_order_relaxed)) {
    std::cerr << "[WARNING] BunchTwissAnalysis::analyzeBunch is deprecated. Use "
                 "computeBunchMoments() instead.\n";
  }
  computeBunchMoments(bunch);
}

/// @brief Templated implementation of computeBunchMoments specialized on macrosize and dispersion
/// flags.
///
/// Computes low order (order <= 2) moments in the first loop over the bunch followed by an early
/// return. If higher order moments (order >= 3) are requested, then a second loop over the bunch
/// will be triggered to compute the remaining moments.
template <bool HasMacrosizeAttr, bool Dispersion>
void BunchTwissAnalysis::computeBunchMomentsImpl(Bunch* bunch, bool normalize, bool emitnormflag)
{
  const auto size = bunch->getSize();

  if (size == 0) {
    return;
  }

  auto* macrosize_attr =
    HasMacrosizeAttr ? static_cast<ParticleMacroSize*>(bunch->getParticleAttributes("macrosize"))
                     : nullptr;
  auto** coords = bunch->coordArr();

  for (int ip = 0; ip < size; ++ip) {
    const double w = HasMacrosizeAttr ? macrosize_attr->macrosize(ip) : 1.0;

    if constexpr (HasMacrosizeAttr) {
      total_macrosize_ += w;
    }

    for (int i = 0; i < N; ++i) {
      avg_arr[i] += w * coords[ip][i];
      for (int j = 0; j <= i; ++j) {
        cov_arr[covIdx(i, j)] += w * coords[ip][i] * coords[ip][j];
      }
    }
  }

  auto comm = bunch->getMPI_Comm_Local()->comm;

  if constexpr (HasMacrosizeAttr) {
    ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, &total_macrosize_, 1, MPI_DOUBLE, MPI_SUM, comm);
  }
  else {
    total_macrosize_ = bunch->getSizeGlobal();
  }

  ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, avg_arr.data(), N, MPI_DOUBLE, MPI_SUM, comm);
  ORBIT_MPI_Allreduce(ORBIT_MPI_IN_PLACE, cov_arr.data(), NN, MPI_DOUBLE, MPI_SUM, comm);

  // <u - uhat><v - vhat> = <u><v> - uhat*vhat
  for (int i = 0; i < N; ++i) {
    avg_arr[i] /= total_macrosize_;
    for (int j = 0; j < i + 1; ++j) {
      const int idx = covIdx(i, j);
      cov_arr[idx] = cov_arr[idx] / total_macrosize_ - avg_arr[i] * avg_arr[j];
      cov_arr[covIdx(j, i)] = cov_arr[covIdx(i, j)]; // covariance matrix is symmetric.
    }
  }

  double inv_xbt = 1.0;
  double inv_ybt = 1.0;

  if (normalize || emitnormflag) {
    inv_xbt /= std::sqrt(getBeta(0) * (emitnormflag ? getEmittance(0) : 1.0));
    inv_ybt /= std::sqrt(getBeta(1) * (emitnormflag ? getEmittance(1) : 1.0));
  }

  const double disp_scale =
    Dispersion ? getDispersion(0) / (bunch_kinenergy_ + bunch_mass_) / (bunch_beta_ * bunch_beta_)
               : 0.0;
  const double xAvg = avg_arr[0] - disp_scale * avg_arr[5];
  const double yAvg = avg_arr[2];

  const int nMoments = order_ + 1;
  momentXY_[momentIdx(0, 0, nMoments)] = 1.0;
  momentXY_[momentIdx(1, 0, nMoments)] = xAvg;
  momentXY_[momentIdx(0, 1, nMoments)] = yAvg;

  if (order_ < 2) {
    return;
  }

  const double x2 = cov_arr[covIdx(0, 0)] - 2.0 * disp_scale * cov_arr[covIdx(0, 5)] +
                    disp_scale * disp_scale * cov_arr[covIdx(5, 5)];
  const double xy = cov_arr[covIdx(0, 2)] - disp_scale * cov_arr[covIdx(2, 5)];
  const double y2 = cov_arr[covIdx(2, 2)];
  momentXY_[momentIdx(2, 0, nMoments)] = x2 * inv_xbt * inv_xbt;
  momentXY_[momentIdx(1, 1, nMoments)] = xy * inv_xbt * inv_ybt;
  momentXY_[momentIdx(0, 2, nMoments)] = y2 * inv_ybt * inv_ybt;

  if (order_ == 2) {
    return;
  }

  for (int ip = 0; ip < size; ++ip) {
    const double w = HasMacrosizeAttr ? macrosize_attr->macrosize(ip) : 1.0;

    const double dx =
      Dispersion ? coords[ip][0] - xAvg - disp_scale * coords[ip][5] : coords[ip][0] - xAvg;
    const double dy = coords[ip][2] - yAvg;

    const double normX = dx * inv_xbt;
    const double normY = dy * inv_ybt;

    double ny = 1.0;
    for (int j = 0; j < order_; ++j) {
      int i_start = (j < 3) ? 3 - j : 0;
      double nx = 1.0;
      for (int ii = 0; ii < i_start; ++ii) {
        nx *= normX;
      }
      for (int i = i_start; i < nMoments - j; ++i) {
        momentXY_[momentIdx(i, j, nMoments)] += w * nx * ny;
        nx *= normX;
      }
      ny *= normY;
    }
  }

  ORBIT_MPI_Allreduce(
    ORBIT_MPI_IN_PLACE,
    momentXY_.data(),
    nMoments * nMoments,
    MPI_DOUBLE,
    MPI_SUM,
    comm
  );

  for (int i = 0; i < nMoments; ++i) {
    for (int j = 0; j < nMoments - i; ++j) {
      if (i + j <= 2) {
        continue;
      }
      momentXY_[momentIdx(i, j, nMoments)] /= total_macrosize_;
    }
  }
}

void BunchTwissAnalysis::computeBunchMoments(
  Bunch* bunch,
  int order,
  bool normalize,
  bool emitnormflag,
  bool dispersionflag
)
{
  avg_arr.fill(0.0);
  cov_arr.fill(0.0);
  total_macrosize_ = 0.0;

  // the number of unique pairs i, j for moments up to a max order, n:
  // (n+1)(n+2)/2
  momentXY_.assign((order + 1)*(order + 2)/2, 0.0);

  bunch->compress();
  count_ = bunch->getSizeGlobal();
  order_ = order;

  SyncPart* syncPart = bunch->getSyncPart();
  bunch_momentum_ = syncPart->getMomentum();
  bunch_beta_ = syncPart->getBeta();
  bunch_gamma_ = syncPart->getGamma();
  bunch_kinenergy_ = syncPart->getEnergy();
  bunch_mass_ = syncPart->getMass();

  if (bunch->hasParticleAttributes("macrosize")) {
    if (dispersionflag) {
      computeBunchMomentsImpl<true, true>(bunch, normalize, emitnormflag);
    }
    else {
      computeBunchMomentsImpl<true, false>(bunch, normalize, emitnormflag);
    }
  }
  else {
    if (dispersionflag) {
      computeBunchMomentsImpl<false, true>(bunch, normalize, emitnormflag);
    }
    else {
      computeBunchMomentsImpl<false, false>(bunch, normalize, emitnormflag);
    }
  }
}

double BunchTwissAnalysis::getCorrelation(int ic, int jc) const
{
  if (ic < 0 || ic > 5 || jc < 0 || jc > 5) {
    return 0.;
  }
  return cov_arr[covIdx(ic, jc)];
}

double BunchTwissAnalysis::getBunchMoment(int i, int j) const
{
  if (i < 0 || i > order_ || j < 0 || j > order_) {
    return 0.;
  }
  return momentXY_[momentIdx(i, j, order_ + 1)];
}

double BunchTwissAnalysis::getAverage(int ic) const
{
  if (ic < 0 || ic > 5) {
    return 0.;
  }
  return avg_arr[ic];
}

int BunchTwissAnalysis::getGlobalCount() const
{
  return count_;
}

double BunchTwissAnalysis::getGlobalMacrosize() const
{
  return total_macrosize_;
}

double BunchTwissAnalysis::getEmittance(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double x2_avg = std::abs(getCorrelation(2 * ic, 2 * ic));
  double xp2_avg = std::abs(getCorrelation(2 * ic + 1, 2 * ic + 1));
  double x_xp_avg = getCorrelation(2 * ic, 2 * ic + 1);
  double x_dE_avg = getCorrelation(2 * ic, 5);
  double xp_dE_avg = getCorrelation(2 * ic + 1, 5);
  double dE2_avg = std::abs(getCorrelation(5, 5));
  if (ic == 2 || dE2_avg == 0) {
    return std::sqrt(std::abs(x2_avg * xp2_avg - x_xp_avg * x_xp_avg));
  }
  double x1 = x2_avg - x_dE_avg * x_dE_avg / dE2_avg;
  double x2 = xp2_avg - xp_dE_avg * xp_dE_avg / dE2_avg;
  double x3 = x_xp_avg - x_dE_avg * xp_dE_avg / dE2_avg;
  return std::sqrt(std::abs(x1 * x2 - x3 * x3));
}

double BunchTwissAnalysis::getEmittanceNormalized(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  if (ic == 2) {
    return getEmittance(ic);
  }
  return getEmittance(ic) * bunch_gamma_ * bunch_beta_;
}

double BunchTwissAnalysis::getAlpha(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double x_xp_avg = getCorrelation(2 * ic, 2 * ic + 1);
  if (ic == 2) {
    return -x_xp_avg / getEmittance(ic);
  }
  double x_dE_avg = getCorrelation(2 * ic, 5);
  double xp_dE_avg = getCorrelation(2 * ic + 1, 5);
  double dE2_avg = std::abs(getCorrelation(5, 5));
  return -(x_xp_avg - x_dE_avg * xp_dE_avg / dE2_avg) / getEmittance(ic);
}

double BunchTwissAnalysis::getBeta(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double x2_avg = std::abs(getCorrelation(2 * ic, 2 * ic));
  if (ic == 2) {
    return x2_avg / getEmittance(ic);
  }
  double x_dE_avg = getCorrelation(2 * ic, 5);
  double dE2_avg = std::abs(getCorrelation(5, 5));
  return (x2_avg - x_dE_avg * x_dE_avg / dE2_avg) / getEmittance(ic);
}

double BunchTwissAnalysis::getGamma(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double xp2_avg = std::abs(getCorrelation(2 * ic + 1, 2 * ic + 1));
  if (ic == 2) {
    return xp2_avg / getEmittance(ic);
  }
  double xp_dE_avg = getCorrelation(2 * ic + 1, 5);
  double dE2_avg = std::abs(getCorrelation(5, 5));
  return (xp2_avg - xp_dE_avg * xp_dE_avg / dE2_avg) / getEmittance(ic);
}

double BunchTwissAnalysis::getDispersion(int ic) const
{
  if (ic < 0 || ic > 1) {
    return 0.;
  }
  double x_dE_avg = getCorrelation(2 * ic, 5);
  double dE2_avg = std::abs(getCorrelation(5, 5));
  return x_dE_avg / dE2_avg * bunch_momentum_ * bunch_beta_;
}

double BunchTwissAnalysis::getDispersionDerivative(int ic) const
{
  if (ic < 0 || ic > 1) {
    return 0.;
  }
  double xp_dE_avg = getCorrelation(2 * ic + 1, 5);
  double dE2_avg = std::abs(getCorrelation(5, 5));
  return xp_dE_avg / dE2_avg * bunch_momentum_ * bunch_beta_;
}

double BunchTwissAnalysis::getEffectiveEmittance(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double x2_avg = std::abs(getCorrelation(2 * ic, 2 * ic));
  double xp2_avg = std::abs(getCorrelation(2 * ic + 1, 2 * ic + 1));
  double x_xp_avg = getCorrelation(2 * ic, 2 * ic + 1);
  return std::sqrt(std::abs(x2_avg * xp2_avg - x_xp_avg * x_xp_avg));
}

double BunchTwissAnalysis::getEffectiveAlpha(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double x2_avg = std::abs(getCorrelation(2 * ic, 2 * ic));
  double xp2_avg = std::abs(getCorrelation(2 * ic + 1, 2 * ic + 1));
  double x_xp_avg = getCorrelation(2 * ic, 2 * ic + 1);
  double emitt2_rms = x2_avg * xp2_avg - x_xp_avg * x_xp_avg;
  if (emitt2_rms <= 0.) {
    return 0.;
  }
  return -x_xp_avg / std::sqrt(emitt2_rms);
}

double BunchTwissAnalysis::getEffectiveBeta(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double x2_avg = std::abs(getCorrelation(2 * ic, 2 * ic));
  double xp2_avg = std::abs(getCorrelation(2 * ic + 1, 2 * ic + 1));
  double x_xp_avg = getCorrelation(2 * ic, 2 * ic + 1);
  double emitt2_rms = x2_avg * xp2_avg - x_xp_avg * x_xp_avg;
  if (emitt2_rms <= 0.) {
    return 0.;
  }
  return x2_avg / std::sqrt(emitt2_rms);
}

double BunchTwissAnalysis::getEffectiveGamma(int ic) const
{
  if (ic < 0 || ic > 2) {
    return 0.;
  }
  double xp2_avg = std::abs(getCorrelation(2 * ic + 1, 2 * ic + 1));
  double x2_avg = std::abs(getCorrelation(2 * ic, 2 * ic));
  double x_xp_avg = getCorrelation(2 * ic, 2 * ic + 1);
  double emitt2_rms = x2_avg * xp2_avg - x_xp_avg * x_xp_avg;
  if (emitt2_rms <= 0.) {
    return std::numeric_limits<double>::max();
  }
  return xp2_avg / std::sqrt(emitt2_rms);
}
