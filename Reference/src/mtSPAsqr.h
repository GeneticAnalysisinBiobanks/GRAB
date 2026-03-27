#ifndef SPAsqr_H
#define SPAsqr_H

// SPAsqr.h -- SPA-squared: multi-tau wrapper delegating to SPAGRMClass per tau

#include <RcppArmadillo.h>
#include <sstream>
#include <cstdio>
#include "mtSPAGRM.h"

// Per-tau family data passed from R at null-model construction time.

class mtSPAsqrClass {
private:

  const arma::vec m_taus;
  std::vector<mtSPAGRMClass> m_SPAGRMobj_vec;
  const std::vector<double> m_MAF_interval;

  const double m_SPA_Cutoff;
  const double m_zeta;
  const double m_tol;

public:

  mtSPAsqrClass(
    arma::vec taus,
    std::vector<arma::vec> residCols,
    std::vector<nsSPAGRM::FamilyData> tauFamilyData,
    std::vector<double> sum_R_nonOutlier_vec,
    std::vector<double> R_GRM_R_nonOutlier_vec,
    std::vector<double> R_GRM_R_TwoSubjOutlier_vec,
    std::vector<double> R_GRM_R_vec,
    std::vector<double> MAF_interval,
    double SPA_Cutoff,
    double zeta,
    double tol
  )
    : m_taus(std::move(taus)),
      m_MAF_interval(std::move(MAF_interval)),
      m_SPA_Cutoff(SPA_Cutoff),
      m_zeta(zeta),
      m_tol(tol)
  {
    int ntaus = m_taus.n_elem;
    m_SPAGRMobj_vec.reserve(ntaus);

    for (int i = 0; i < ntaus; ++i) {
      m_SPAGRMobj_vec.emplace_back(
        std::move(residCols[i]),
        sum_R_nonOutlier_vec[i],
        R_GRM_R_nonOutlier_vec[i],
        R_GRM_R_TwoSubjOutlier_vec[i],
        R_GRM_R_vec[i],
        m_MAF_interval,
        std::move(tauFamilyData[i]),
        m_SPA_Cutoff,
        m_zeta,
        m_tol
      );
    }
  }

  int get_ntaus() const { return m_taus.n_elem; }
  std::vector<double> getTaus() const { return arma::conv_to<std::vector<double>>::from(m_taus); }

  std::string getHeaderColumns() const {
    std::ostringstream oss;
    auto taus = getTaus();
    for (double tau : taus) {
      char buf[32]; std::snprintf(buf, sizeof(buf), "%.6g", tau);
      oss << "\tZ_tau" << buf;
    }
    for (double tau : taus) {
      char buf[32]; std::snprintf(buf, sizeof(buf), "%.6g", tau);
      oss << "\tP_tau" << buf;
    }
    oss << "\tP_CCT";
    return oss.str();
  }

  arma::vec getMarkerPval(
    const arma::vec& GVec,
    double altFreq,
    arma::vec& zScoreVec
  ) {
    int ntaus = m_taus.n_elem;
    arma::vec pvalVec(ntaus);
    zScoreVec.set_size(ntaus);

    for (int i = 0; i < ntaus; i++) {

      double zScore_i;
      double pval_i = m_SPAGRMobj_vec[i].getMarkerPval(GVec, altFreq, zScore_i);
      pvalVec(i) = pval_i;
      zScoreVec(i) = zScore_i;
    }

    return pvalVec;
  }

  // Fills rv with [z_1..z_ntaus, p_1..p_ntaus, pCCT]
  void getResultVec(const arma::vec& GVec, double altFreq, std::vector<double>& rv) {
    arma::vec zT;
    arma::vec pT = getMarkerPval(GVec, altFreq, zT);
    int nt = m_taus.n_elem;
    rv.clear();
    rv.reserve(2 * nt + 1);
    for (int j = 0; j < nt; ++j) rv.push_back(zT[j]);
    for (int j = 0; j < nt; ++j) rv.push_back(pT[j]);
    // CCT p-value
    std::vector<double> pvals(nt);
    for (int j = 0; j < nt; ++j) pvals[j] = pT[j];
    std::vector<double> pp;
    pp.reserve(nt);
    for (double x : pvals) if (!std::isnan(x)) pp.push_back(x);
    double pCCT = std::numeric_limits<double>::quiet_NaN();
    if (!pp.empty()) {
      double tStat = 0.0;
      bool zero = false;
      for (double x : pp) {
        if (x <= 0.0) { zero = true; break; }
        double xc = (x >= 1.0) ? 0.999 : x;
        tStat += std::tan((0.5 - xc) * M_PI);
      }
      if (zero) pCCT = 0.0;
      else {
        tStat /= static_cast<double>(pp.size());
        pCCT = (tStat > 1e15) ? (1.0 / tStat) / M_PI : 0.5 - std::atan(tStat) / M_PI;
      }
    }
    rv.push_back(pCCT);
  }

  int resultSize() const { return 2 * static_cast<int>(m_taus.n_elem) + 1; }
};

#endif
