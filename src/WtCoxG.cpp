// WtCoxG.cpp -- WtCoxGClass method implementations

#include "WtCoxG.h"
#include <cmath>
#include <algorithm>
#include <functional>
#include <cstdio>
#include <limits>


static double bvnCdf (double dh, double dk, double r) {
  if (std::abs(r) < 1e-15) {
    return 0.5 * std::erfc(-dh / std::sqrt(2.0)) * 0.5 * std::erfc(-dk / std::sqrt(2.0));
  }

  static const double w6[] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910};
  static const double x6[] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969};

  double h = dh, k = dk, hk = h * k;
  double bvn = 0.0;
  double hs, asr, sn, as_, a, b, c, d, bs, xs;

  if (std::abs(r) < 0.925) {
    hs = (h * h + k * k) / 2.0;
    asr = std::asin(r);
    for (int i = 0; i < 3; ++i) {
      for (int is = -1; is <= 1; is += 2) {
        sn = std::sin(asr * (is * x6[i] + 1.0) / 2.0);
        bvn += w6[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
      }
    }
    bvn *= asr / (4.0 * M_PI);
    bvn += 0.5 * std::erfc(-h / std::sqrt(2.0)) * 0.5 * std::erfc(-k / std::sqrt(2.0));
  } else {
    if (r < 0.0) {
      k = -k; hk = -hk;
    }
    if (std::abs(r) < 1.0) {
      as_ = (1.0 - r) * (1.0 + r);
      a = std::sqrt(as_);
      bs = (h - k) * (h - k);
      c = (4.0 - hk) / 8.0;
      d = (12.0 - hk) / 16.0;
      asr = -(bs / as_ + hk) / 2.0;
      if (asr > -100.0)
        bvn = a * std::exp(asr) * (1.0 - c * (bs - as_) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as_ * as_ / 5.0);
      if (-hk < 100.0) {
        b = std::sqrt(bs);
        bvn -= std::exp(-hk / 2.0) * std::sqrt(2.0 * M_PI) * 0.5 * std::erfc(-b / (std::sqrt(2.0) * a)) *
                  b * (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
      }
      a /= 2.0;
      for (int i = 0; i < 3; ++i) {
        for (int is = -1; is <= 1; is += 2) {
          xs = (a * (is * x6[i] + 1.0));
          xs *= xs;
          asr = -(bs / xs + hk) / 2.0;
          if (asr > -100.0)
          bvn += a * w6[i] * std::exp(asr) * (std::exp(-hk * (1.0 - r) / (2.0 * (1.0 + (is * x6[i] + 1.0) * a * a / as_))) - 1.0 - c * xs * (1.0 + d * xs));
        }
      }
      bvn = -bvn / (2.0 * M_PI);
    }
    if (r > 0.0) {
        double phih = 0.5 * std::erfc(-h / std::sqrt(2.0));
        double phik = 0.5 * std::erfc(-k / std::sqrt(2.0));
        bvn += phih + phik - 1.0 + (phih < phik ? phih : phik);
        if (bvn < 0.0) bvn = 0.0;
    } else {
      bvn = -bvn;
      double phih = 0.5 * std::erfc(-h / std::sqrt(2.0));
      double phik = 0.5 * std::erfc(-k / std::sqrt(2.0));
      if (phih - phik >= 0.0) {
        bvn = phih - phik - bvn;
        if (bvn < 0.0) bvn = 0.0;
      } else {
        bvn = 0.0;
      }
    }
  }
  return bvn;
}

static double pmvnorm2d (double lo1, double hi1,
                          double lo2, double hi2,
                          double var1, double cov12, double var2) {
  double sd1 = std::sqrt(var1);
  double sd2 = std::sqrt(var2);
  double rho = cov12 / (sd1 * sd2);

  auto standardise = [](double v, double sd) -> double {
  if (std::isinf(v)) return v > 0 ? 1e15 : -1e15;
  return v / sd;
  };

  double a1 = standardise(lo1, sd1), b1 = standardise(hi1, sd1);
  double a2 = standardise(lo2, sd2), b2 = standardise(hi2, sd2);

  double p = bvnCdf(b1, b2, rho)
           - bvnCdf(a1, b2, rho)
           - bvnCdf(b1, a2, rho)
           + bvnCdf(a1, a2, rho);
  return std::max(0.0, std::min(1.0, p));
}

namespace WtCoxG {

double WtCoxGClass::findRootBrent(std::function<double(double)> f, double a, double b, double tol) {

  double fa = f(a);
  double fb = f(b);

  if (fa * fb > 0) {
    throw std::runtime_error("Root not bracketed");
  }

  if (std::abs(fa) < std::abs(fb)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }

  double c = a;
  double fc = fa;

  const int max_iter = std::min(50, std::max(15, static_cast<int>(20.0 / tol)));
  tol = std::max(tol, 1e-9);

  double last_improvement = std::abs(b - a);
  int stagnant_iterations = 0;

  for (int iter = 0; iter < max_iter; ++iter) {
    double current_gap = std::abs(b - a);
    if (current_gap < tol || std::abs(fb) < tol * 10) {
      return b;
    }

    if (current_gap >= 0.95 * last_improvement) {
      stagnant_iterations++;
      if (stagnant_iterations > 3)
        break;
    } else {
      stagnant_iterations = 0;
    }
    last_improvement = current_gap;

    double s;
    if (fa != fc && fb != fc && std::abs(fa - fc) > 1e-15 &&
      std::abs(fb - fc) > 1e-15) {
      double denom1 = (fa - fb) * (fa - fc);
      double denom2 = (fb - fa) * (fb - fc);
      double denom3 = (fc - fa) * (fc - fb);

      if (std::abs(denom1) > 1e-12 && std::abs(denom2) > 1e-12 &&
        std::abs(denom3) > 1e-12) {
        s = a * fb * fc / denom1 + b * fa * fc / denom2 +
          c * fa * fb / denom3;
      } else {
        s = b - fb * (b - a) / (fb - fa);
      }
    } else {
      s = b - fb * (b - a) / (fb - fa);
    }

    if (s <= std::min(a, b) || s >= std::max(a, b)) {
      s = (a + b) / 2.0;
    }

    double fs = f(s);
    c = b;
    fc = fb;

    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    if (std::abs(fa) < std::abs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
  }

  return b;
}


double WtCoxGClass::hOrg(
  double t, const arma::vec& R, double MAF, double n_ext,
  double N_all, double sumR, double var_mu_ext,
  double g_var_est, double meanR, double b
) {

  double mu_adj = -2.0 * b * sumR * MAF;
  double var_adj = 4.0 * b * b * sumR * sumR * var_mu_ext;

  double result = 0.0;
  for (size_t i = 0; i < R.n_elem; ++i) {
    result += kG0(t * (R(i) - (1.0 - b) * meanR), MAF);
  }

  return result + mu_adj * t + var_adj * t * t / 2.0;
}

double WtCoxGClass::h1Adj(
  double t, const arma::vec& R, double s, double MAF,
  double n_ext, double N_all, double sumR, double var_mu_ext,
  double g_var_est, double meanR, double b
) {

  double mu_adj = -2.0 * b * sumR * MAF;
  double var_adj = 4.0 * b * b * sumR * sumR * var_mu_ext;

  double result = 0.0;
  for (size_t i = 0; i < R.n_elem; ++i) {
    double R_adj = R(i) - (1.0 - b) * meanR;
    result += R_adj * kG1(t * R_adj, MAF);
  }

  return result + mu_adj + var_adj * t - s;
}

double WtCoxGClass::h2(
  double t, const arma::vec& R, double MAF, double n_ext,
  double N_all, double sumR, double var_mu_ext,
  double g_var_est, double meanR, double b
) {

  double var_adj = n_ext * std::pow(sumR / N_all, 2.0) * 2.0 * MAF * (1.0 - MAF);

  double result = 0.0;
  for (size_t i = 0; i < R.n_elem; ++i) {
    double R_adj = R(i) - (1.0 - b) * meanR;
    result += R_adj * R_adj * kG2(t * R_adj, MAF);
  }

  return result + var_adj;
}


double WtCoxGClass::getProbSpaG(
  double MAF, const arma::vec& R, double s, double n_ext,
  double N_all, double sumR, double var_mu_ext,
  double g_var_est, double meanR, double b, bool lower_tail
) {

  auto h1_func = [&](double t) {
    return h1Adj(t, R, s, MAF, n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
  };

  double zeta;
  try {

    double a = -1.0, b_bound = 1.0;
    double fa = h1_func(a);
    double fb = h1_func(b_bound);

    if (fa * fb > 0) {
      double factor = 2.0;
      int max_extend = 10;
      for (int i = 0; i < max_extend; ++i) {
        if (std::abs(fa) < std::abs(fb)) {

          a = a * factor;
          fa = h1_func(a);
        } else {

          b_bound = b_bound * factor;
          fb = h1_func(b_bound);
        }
        if (fa * fb <= 0) break;
      }
    }
    zeta = findRootBrent(h1_func, a, b_bound, 1e-8);
  } catch (...) {

    std::fprintf(stderr, "[WARN] WtCoxG root finding failed, returning NaN\n");
    std::fflush(stderr);
    return arma::datum::nan;
  }

  double k1 = hOrg(zeta, R, MAF, n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
  double k2 = h2(zeta, R, MAF, n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
  double temp1 = zeta * s - k1;

  double w = (zeta >= 0 ? 1.0 : -1.0) * std::sqrt(2.0 * temp1);
  double v = zeta * std::sqrt(k2);
  double pval = pnormBoost(w + (1.0 / w) * std::log(v / w), 0.0, 1.0, lower_tail, false);

  return pval;
}


arma::vec WtCoxGClass::spaGOneSnpHomo(
  const arma::vec& g_input, const arma::vec& R,
  double mu_ext, double n_ext, double b,
  double sigma2, double var_ratio, double SPA_Cutoff,
  double missing_cutoff, double min_mac
) {

  arma::vec g = imputeMissing(g_input);
  arma::uvec missing_idx = arma::find_nonfinite(g_input);
  double missing_rate = static_cast<double>(missing_idx.n_elem) / g_input.n_elem;

  if (std::isnan(mu_ext)) {
    mu_ext = 0.0;
    n_ext = 0.0;
  }

  double sum_g = arma::sum(g);
  double sum_2_minus_g = arma::sum(2.0 - g);

  if (sum_g < min_mac || sum_2_minus_g < min_mac || missing_rate > missing_cutoff) {
    return arma::vec({arma::datum::nan, arma::datum::nan, arma::datum::nan, arma::datum::nan});
  }

  double N = g.n_elem;
  double mu_int = arma::mean(g) / 2.0;
  double MAF = (1.0 - b) * mu_int + b * mu_ext;
  double sumR = arma::sum(R);
  double N_all = N + n_ext;
  double S = arma::sum(R % (g - 2.0 * MAF));
  double S_raw = S;
  S = S / var_ratio;

  double g_var_est = 2.0 * MAF * (1.0 - MAF);
  double var_mu_ext = (n_ext == 0.0) ? 0.0 : (MAF * (1.0 - MAF) / (2.0 * n_ext) + sigma2);

  double meanR = arma::mean(R);
  arma::vec R_adj = R - (1.0 - b) * meanR;
  double S_var = arma::sum(R_adj % R_adj) * g_var_est + 4.0 * b * b * sumR * sumR * var_mu_ext;
  double z = S / std::sqrt(S_var);

  if (std::abs(z) < SPA_Cutoff) {
    double pval_norm = 2.0 * pnormBoost(-std::abs(z), 0.0, 1.0, true, false);
    pval_norm = std::min(1.0, pval_norm);
    return arma::vec({pval_norm, pval_norm, S_raw, z});
  }
  else {
    double pval1 = getProbSpaG(
      MAF, R, std::abs(S), n_ext, N_all, sumR,
      var_mu_ext, g_var_est, meanR, b, false
    );
    double pval2 = getProbSpaG(
      MAF, R, -std::abs(S), n_ext, N_all, sumR,
      var_mu_ext, g_var_est, meanR, b, true
    );

    double pval_spa = pval1 + pval2;
    pval_spa = std::min(1.0, pval_spa);
    double pval_norm = 2.0 * pnormBoost(-std::abs(z), 0.0, 1.0, true, false);
    pval_norm = std::min(1.0, pval_norm);

    return arma::vec({pval_spa, pval_norm, S_raw, z});
  }
}


arma::vec WtCoxGClass::wtCoxGTest(
  const arma::vec& g_input, const arma::vec& R, const arma::vec& w,
  double p_bat, double TPR, double sigma2, double b,
  double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
  double var_ratio0, double var_ratio1, double mu_ext,
  double n_ext, double p_cut
) {

  arma::vec g = g_input;
  arma::uvec na_indices = arma::find_nonfinite(g);
  double missing_rate = static_cast<double>(na_indices.n_elem) / g.n_elem;

  if (missing_rate != 0.0) {
    arma::uvec valid_indices = arma::find_finite(g);
    if (valid_indices.n_elem > 0) {
      double mean_val = arma::mean(g.elem(valid_indices));
      g.elem(na_indices).fill(mean_val);
    }
  }

  if (std::isnan(mu_ext)) {
    double var_ratio_to_use;
    if (std::isnan(TPR) && std::isnan(sigma2)) {
      var_ratio_to_use = var_ratio_int;
    } else {
      var_ratio_to_use = 1.0;
    }
    arma::vec spa_result = spaGOneSnpHomo(
      g, R, 0.0, 0.0, 0.0, 0.0, var_ratio_to_use, m_SPA_Cutoff, 0.15, 10.0);
    return arma::vec({spa_result(0), spa_result(2), spa_result(3)});
  }

  double sum_g = arma::sum(g);
  double sum_2_minus_g = arma::sum(2.0 - g);

  if (p_bat < p_cut || std::isnan(p_bat) || sum_g < 10 || sum_2_minus_g < 10) {
    return arma::vec({arma::datum::nan, arma::datum::nan, arma::datum::nan});
  }

  double meanR = arma::mean(R);
  double sumR = arma::sum(R);
  double mu_int = arma::mean(g) / 2.0;

  double mu = (1.0 - b) * mu_int + b * mu_ext;
  double S = arma::sum(R % (g - 2.0 * mu));

  arma::vec w1 = w / (2.0 * arma::sum(w));
  double var_mu_ext = mu * (1.0 - mu) / (2.0 * n_ext);
  double var_Sbat = arma::sum(w1 % w1) * 2.0 * mu * (1.0 - mu) + var_mu_ext;

  double qnorm_val = qnormBoost(1.0 - p_cut / 2.0, 0.0, 1.0, true, false);
  double lb = -qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);
  double ub = qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);

  double c = pnormBoost(ub / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
  double d = pnormBoost(lb / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
  double p_deno = TPR * (std::exp(d) * (std::exp(c - d) - 1.0)) + (1.0 - TPR) * (1.0 - p_cut);

  arma::vec spa_result_s0 = spaGOneSnpHomo(
    g, R, mu_ext, n_ext, b, 0.0, var_ratio0, m_SPA_Cutoff, 0.15, 10.0);
  double p_spa_s0 = spa_result_s0(0);
  double qchisq_val = qchisqBoost(p_spa_s0, 1.0, false, false);
  double var_S = S * S / var_ratio0 / qchisq_val;

  arma::vec R_minus_factor = R - (1.0 - b) * meanR;
  double var_int = arma::sum(R_minus_factor % R_minus_factor) * 2.0 * mu * (1.0 - mu);

  double cov_Sbat_S = arma::sum(w1 % R_minus_factor) * 2.0 * mu * (1.0 - mu) +
            2.0 * b * sumR * var_mu_ext;
  double denominator = var_int + 4.0 * b * b * sumR * sumR * var_mu_ext;
  if (denominator <= 0.0) {
    return arma::vec({arma::datum::nan, arma::datum::nan, arma::datum::nan});
  }

  cov_Sbat_S = cov_Sbat_S * std::sqrt(var_S / denominator);
  double z = S / std::sqrt(var_S);

  arma::mat VAR(2, 2);
  VAR(0, 0) = var_S;
  VAR(0, 1) = cov_Sbat_S;
  VAR(1, 0) = cov_Sbat_S;
  VAR(1, 1) = var_Sbat;

  double negInf = -std::numeric_limits<double>::infinity();
  double p0 = pmvnorm2d(
    negInf, -std::abs(S / std::sqrt(var_ratio0)),
    lb / std::sqrt(var_ratio_w0), ub / std::sqrt(var_ratio_w0),
    VAR(0, 0), VAR(0, 1), VAR(1, 1));
  p0 = std::max(0.0, std::min(1.0, p0));

  arma::vec spa_result_s1 = spaGOneSnpHomo(g, R, mu_ext, n_ext, b, sigma2, var_ratio1, m_SPA_Cutoff, 0.15, 10.0);
  double p_spa_s1 = spa_result_s1(0);
  double var_S1 = S * S / var_ratio1 / qchisqBoost(p_spa_s1, 1.0, false, false);

  double cov_Sbat_S1 = arma::sum(w1 % R_minus_factor) * 2.0 * mu * (1.0 - mu) + 2.0 * b * sumR * (var_mu_ext + sigma2);
  cov_Sbat_S1 = cov_Sbat_S1 * std::sqrt(var_S1 / denominator);
  double var_Sbat1 = var_Sbat + sigma2;

  arma::mat VAR1(2, 2);
  VAR1(0, 0) = var_S1;
  VAR1(0, 1) = cov_Sbat_S1;
  VAR1(1, 0) = cov_Sbat_S1;
  VAR1(1, 1) = var_Sbat1;

  double p1 = pmvnorm2d(
    negInf, -std::abs(S / std::sqrt(var_ratio1)),
    lb / std::sqrt(var_ratio_w1), ub / std::sqrt(var_ratio_w1),
    VAR1(0, 0), VAR1(0, 1), VAR1(1, 1)
  );

  p1 = std::max(0.0, std::min(1.0, p1));
  double p_con = 2.0 * (TPR * p1 + (1.0 - TPR) * p0) / p_deno;

  return arma::vec({p_con, S, z});
}

}
