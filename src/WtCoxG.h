#ifndef WTCOXG_H
#define WTCOXG_H

// WtCoxG.h -- Weighted Cox-type G-test for batch-effect-aware association

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

namespace WtCoxG {
class WtCoxGClass {
private:
  // Per-marker external-reference metadata.
  struct MarkerInfo {
    double AF_ref;
    double AN_ref;
    double TPR;
    double sigma2;
    double pvalue_bat;
    double w_ext;
    double var_ratio_w0;
    double var_ratio_int;
    double var_ratio_ext;


    MarkerInfo() = default;


    MarkerInfo(double af_ref, double an_ref, double tpr, double sig2,
        double pval_bat, double w_ext_val, double var_ratio_w0_val,
        double var_ratio_int_val, double var_ratio_ext_val)
      : AF_ref(af_ref), AN_ref(an_ref), TPR(tpr), sigma2(sig2),
      pvalue_bat(pval_bat), w_ext(w_ext_val), var_ratio_w0(var_ratio_w0_val),
      var_ratio_int(var_ratio_int_val), var_ratio_ext(var_ratio_ext_val) {}
  };

  std::vector<MarkerInfo> m_markerInfoVec;
  arma::vec m_R;
  arma::vec m_w;
  double m_cutoff;
  double m_SPA_Cutoff;
  arma::vec m_scoreVec;
  arma::vec m_zScoreVec;

public:

  WtCoxGClass(
    const arma::vec& R,
    const arma::vec& w,
    const double cutoff,
    const double SPA_Cutoff
  )
    : m_R(R),
    m_w(w),
    m_cutoff(cutoff),
    m_SPA_Cutoff(SPA_Cutoff) {

  }


  void updateMarkerInfo(const std::vector<double>& AF_ref,
                          const std::vector<double>& AN_ref,
                          const std::vector<double>& TPR,
                          const std::vector<double>& sigma2,
                          const std::vector<double>& pvalue_bat,
                          const std::vector<double>& w_ext,
                          const std::vector<double>& var_ratio_w0,
                          const std::vector<double>& var_ratio_int,
                          const std::vector<double>& var_ratio_ext) {
    m_markerInfoVec.clear();

    size_t n = AF_ref.size();
    if (AN_ref.size() != n || TPR.size() != n || sigma2.size() != n ||
      pvalue_bat.size() != n || w_ext.size() != n || var_ratio_w0.size() != n ||
      var_ratio_int.size() != n || var_ratio_ext.size() != n) {
      throw std::runtime_error("WtCoxG::updateMarkerInfo received vectors with inconsistent lengths.");
    }

    m_markerInfoVec.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      m_markerInfoVec.emplace_back(
        AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
        w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
      );
    }
  }


  arma::vec getpvalVec(const arma::vec& GVec, const int i) {
    arma::vec result(2);
    const MarkerInfo& info = m_markerInfoVec[i];


    arma::vec res_ext = wtCoxGTest(
      GVec, m_R, m_w,
      info.pvalue_bat,
      info.TPR,
      info.sigma2,
      info.w_ext,
      info.var_ratio_int,
      info.var_ratio_w0,
      info.var_ratio_w0,
      info.var_ratio_ext,
      info.var_ratio_ext,
      info.AF_ref,
      info.AN_ref / 2.0,
      m_cutoff);
    result(0) = res_ext(0);


    arma::vec res_noext = wtCoxGTest(
      GVec, m_R, m_w,
      info.pvalue_bat,
      arma::datum::nan,
      arma::datum::nan,
      0.0,
      info.var_ratio_int,
      1.0,
      1.0,
      1.0,
      1.0,
      arma::datum::nan,
      arma::datum::nan,
      m_cutoff);
    result(1) = res_noext(0);


    if (m_scoreVec.n_elem < 2) {
      m_scoreVec.resize(2);
      m_zScoreVec.resize(2);
    }
    m_scoreVec(0) = res_ext(1);
    m_scoreVec(1) = res_noext(1);
    m_zScoreVec(0) = res_ext(2);
    m_zScoreVec(1) = res_noext(2);

    return result;
  }


  arma::vec getScoreVec() const {
    return m_scoreVec;
  }


  arma::vec getZScoreVec() const {
    return m_zScoreVec;
  }

private:


  // ---- Inline helpers ----
  inline arma::vec imputeMissing(const arma::vec& g) {
    arma::vec g_imputed = g;
    arma::uvec missing_idx = arma::find_nonfinite(g);

    if (missing_idx.n_elem > 0) {
      arma::uvec non_missing_idx = arma::find_finite(g);
      if (non_missing_idx.n_elem > 0) {
        double mean_val = arma::mean(g.elem(non_missing_idx));
        g_imputed.elem(missing_idx).fill(mean_val);
      }
    }

    return g_imputed;
  }


  inline double pnormBoost(double x, double mean = 0.0, double sd = 1.0, bool lower_tail = true, bool log_p = false) {
    boost::math::normal dist(mean, sd);
    double result = boost::math::cdf(dist, x);
    if (!lower_tail) result = 1.0 - result;
    if (log_p) result = std::log(result);
    return result;
  }

  inline double qnormBoost(double p, double mean = 0.0, double sd = 1.0, bool lower_tail = true, bool log_p = false) {
    if (log_p) p = std::exp(p);
    if (!lower_tail) p = 1.0 - p;

    p = std::max(1e-300, std::min(1.0 - 1e-15, p));
    boost::math::normal dist(mean, sd);
    return boost::math::quantile(dist, p);
  }

  inline double qchisqBoost(double p, double df, bool lower_tail = true, bool log_p = false) {
    if (log_p) p = std::exp(p);
    if (!lower_tail) p = 1.0 - p;


    p = std::max(1e-300, std::min(1.0 - 1e-15, p));
    boost::math::chi_squared dist(df);
    return boost::math::quantile(dist, p);
  }


  // Moment-generating function components M(t) and cumulant K(t) for genotype.
  inline double mG0(double t, double MAF) {
    return std::pow(1.0 - MAF + MAF * std::exp(t), 2.0);
  }

  inline double mG1(double t, double MAF) {
    return 2.0 * (MAF * std::exp(t)) * (1.0 - MAF + MAF * std::exp(t));
  }

  inline double mG2(double t, double MAF) {
    double maf_exp_t = MAF * std::exp(t);
    return 2.0 * maf_exp_t * maf_exp_t + 2.0 * maf_exp_t * (1.0 - MAF + maf_exp_t);
  }

  inline double kG0(double t, double MAF) {
    return std::log(mG0(t, MAF));
  }

  inline double kG1(double t, double MAF) {
    return mG1(t, MAF) / mG0(t, MAF);
  }

  inline double kG2(double t, double MAF) {
    double m0 = mG0(t, MAF);
    double m1 = mG1(t, MAF);
    double m2 = mG2(t, MAF);
    return (m0 * m2) / (m0 * m0) - std::pow(m1 / m0, 2.0);
  }


  // ---- Core SPA routines ----
  double findRootBrent(std::function<double(double)> f, double a, double b, double tol = 1e-6);

  double hOrg(double t, const arma::vec& R, double MAF, double n_ext,
          double N_all, double sumR, double var_mu_ext,
          double g_var_est, double meanR, double b);

  double h1Adj(double t, const arma::vec& R, double s, double MAF,
          double n_ext, double N_all, double sumR, double var_mu_ext,
          double g_var_est, double meanR, double b);

  double h2(double t, const arma::vec& R, double MAF, double n_ext,
        double N_all, double sumR, double var_mu_ext,
        double g_var_est, double meanR, double b);

  double getProbSpaG(double MAF, const arma::vec& R, double s, double n_ext,
              double N_all, double sumR, double var_mu_ext,
              double g_var_est, double meanR, double b, bool lower_tail);

  arma::vec spaGOneSnpHomo(const arma::vec& g_input, const arma::vec& R,
                  double mu_ext, double n_ext, double b,
                  double sigma2, double var_ratio, double Cutoff,
                  double missing_cutoff, double min_mac);

  arma::vec wtCoxGTest(const arma::vec& g_input, const arma::vec& R, const arma::vec& w,
            double p_bat, double TPR, double sigma2, double b,
            double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
            double var_ratio0, double var_ratio1, double mu_ext,
            double n_ext, double p_cut);
};

}


#endif
