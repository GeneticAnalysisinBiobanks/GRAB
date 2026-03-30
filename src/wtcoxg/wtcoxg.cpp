// wtcoxg.cpp — WtCoxG full implementation (Phases 1–3)

#include "wtcoxg/wtcoxg.hpp"
#include "io/plink.hpp"
#include "io/resid_file.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using NaN = std::numeric_limits<double>;

// ======================================================================
// Anonymous-namespace helpers (SPA internals, ported from Reference)
// ======================================================================

namespace {

double hOrg(double t, const Eigen::VectorXd& R, double MAF,
            double n_ext, double N_all, double sumR,
            double var_mu_ext, double /*g_var_est*/,
            double meanR, double b) {
  double mu_adj  = -2.0 * b * sumR * MAF;
  double var_adj =  4.0 * b * b * sumR * sumR * var_mu_ext;
  double result = 0.0;
  const double bm = (1.0 - b) * meanR;
  for (Eigen::Index i = 0; i < R.size(); ++i)
    result += math::kG0(t * (R[i] - bm), MAF);
  return result + mu_adj * t + var_adj * t * t / 2.0;
}

double h1Adj(double t, const Eigen::VectorXd& R, double s, double MAF,
             double n_ext, double N_all, double sumR,
             double var_mu_ext, double /*g_var_est*/,
             double meanR, double b) {
  double mu_adj  = -2.0 * b * sumR * MAF;
  double var_adj =  4.0 * b * b * sumR * sumR * var_mu_ext;
  double result = 0.0;
  const double bm = (1.0 - b) * meanR;
  for (Eigen::Index i = 0; i < R.size(); ++i) {
    double R_adj = R[i] - bm;
    result += R_adj * math::kG1(t * R_adj, MAF);
  }
  return result + mu_adj + var_adj * t - s;
}

double h2(double t, const Eigen::VectorXd& R, double MAF,
          double n_ext, double N_all, double sumR,
          double var_mu_ext, double /*g_var_est*/,
          double meanR, double b) {
  double var_adj = n_ext * (sumR / N_all) * (sumR / N_all) * 2.0 * MAF * (1.0 - MAF);
  double result = 0.0;
  const double bm = (1.0 - b) * meanR;
  for (Eigen::Index i = 0; i < R.size(); ++i) {
    double R_adj = R[i] - bm;
    result += R_adj * R_adj * math::kG2(t * R_adj, MAF);
  }
  return result + var_adj;
}

double getProbSpaG(double MAF, const Eigen::VectorXd& R, double s,
                   double n_ext, double N_all, double sumR,
                   double var_mu_ext, double g_var_est,
                   double meanR, double b, bool lower_tail) {
  auto h1_func = [&](double t) {
    return h1Adj(t, R, s, MAF, n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
  };
  double zeta;
  try {
    double a = -1.0, bb = 1.0;
    double fa = h1_func(a), fb = h1_func(bb);
    if (fa * fb > 0) {
      double factor = 2.0;
      for (int i = 0; i < 10; ++i) {
        if (std::abs(fa) < std::abs(fb)) { a *= factor; fa = h1_func(a); }
        else                             { bb *= factor; fb = h1_func(bb); }
        if (fa * fb <= 0) break;
      }
    }
    zeta = math::findRootBrent(h1_func, a, bb, 1e-8);
  } catch (...) {
    return NaN::quiet_NaN();
  }

  double k1 = hOrg(zeta, R, MAF, n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
  double k2 = h2(zeta, R, MAF, n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
  double temp1 = zeta * s - k1;
  if (!std::isfinite(zeta) || temp1 < 0.0 || k2 <= 0.0)
    return NaN::quiet_NaN();

  double w = (zeta >= 0 ? 1.0 : -1.0) * std::sqrt(2.0 * temp1);
  double v = zeta * std::sqrt(k2);
  if (w == 0.0 || v == 0.0 || (v / w) <= 0.0)
    return NaN::quiet_NaN();

  return math::pnorm(w + (1.0 / w) * std::log(v / w), 0.0, 1.0, lower_tail, false);
}

struct SpaResult { double pval, pval2, score, zscore; };

SpaResult spaGOneSnpHomo(
    const Eigen::Ref<const Eigen::VectorXd>& g,
    const Eigen::VectorXd& R,
    double mu_ext, double n_ext, double b,
    double sigma2, double var_ratio, double SPA_Cutoff) {

  if (std::isnan(mu_ext)) { mu_ext = 0.0; n_ext = 0.0; }

  const double N = static_cast<double>(g.size());
  double mu_int = g.mean() / 2.0;
  double MAF = std::clamp((1.0 - b) * mu_int + b * mu_ext, 0.0, 1.0);
  double sumR = R.sum();
  double N_all = N + n_ext;
  double S = (R.array() * (g.array() - 2.0 * MAF)).sum();
  double S_raw = S;
  S /= var_ratio;

  double g_var_est = 2.0 * MAF * (1.0 - MAF);
  double var_mu_ext = (n_ext == 0.0) ? 0.0 : (MAF * (1.0 - MAF) / (2.0 * n_ext) + sigma2);

  double meanR = R.mean();
  Eigen::ArrayXd R_adj = R.array() - (1.0 - b) * meanR;
  double S_var = (R_adj * R_adj).sum() * g_var_est + 4.0 * b * b * sumR * sumR * var_mu_ext;

  if (S_var <= 0.0)
    return {NaN::quiet_NaN(), NaN::quiet_NaN(), S_raw, NaN::quiet_NaN()};

  double z = S / std::sqrt(S_var);
  if (std::abs(z) < SPA_Cutoff) {
    double pval_norm = std::min(1.0, 2.0 * math::pnorm(-std::abs(z)));
    return {pval_norm, pval_norm, S_raw, z};
  }

  double pval1 = getProbSpaG(MAF, R,  std::abs(S), n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b, false);
  double pval2 = getProbSpaG(MAF, R, -std::abs(S), n_ext, N_all, sumR, var_mu_ext, g_var_est, meanR, b, true);
  double pval_spa = std::min(1.0, pval1 + pval2);
  return {pval_spa, std::min(1.0, 2.0 * math::pnorm(-std::abs(z))), S_raw, z};
}

} // anon namespace


// ======================================================================
// Phase 1 — File parsing & marker matching
// ======================================================================

std::vector<RefAfRecord> loadRefAfFile(const std::string& filename) {
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open ref-af file: " + filename);

  std::vector<RefAfRecord> recs;
  recs.reserve(100000);
  uint32_t lineNo = 0;
  std::string line;
  while (std::getline(ifs, line)) {
    ++lineNo;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    // Whitespace-delimited parse: CHROM  POS  A1  A2  A1F  N
    const char*       p    = line.c_str();
    const char* const lEnd = p + line.size();
    auto skipWS  = [&]() { while (p < lEnd && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < lEnd && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };
    auto err = [&](const char* msg) -> std::runtime_error {
      return std::runtime_error(filename + " line " + std::to_string(lineNo) + ": " + msg);
    };

    skipWS();
    if (p >= lEnd) continue;  // blank / all-whitespace line

    // Columns: CHROM  POS  A1  A2  A1F  N
    RefAfRecord r;
    r.chrom = nextTok();
    if (r.chrom.empty()) throw err("missing CHROM field");

    char* endPtr;
    const std::string pos_str = nextTok();
    if (pos_str.empty()) throw err("missing POS field");
    r.pos = static_cast<uint32_t>(std::strtoul(pos_str.c_str(), &endPtr, 10));
    if (endPtr == pos_str.c_str()) throw err("invalid POS field");

    r.a1 = nextTok();
    r.a2 = nextTok();
    if (r.a1.empty() || r.a2.empty()) throw err("missing A1 or A2 field");

    const std::string af_str = nextTok();
    if (af_str.empty()) throw err("missing A1F field");
    r.AF_ref = std::strtod(af_str.c_str(), &endPtr);
    if (endPtr == af_str.c_str()) throw err("invalid A1F field");

    const std::string n_str = nextTok();
    if (n_str.empty()) throw err("missing N field");
    double N = std::strtod(n_str.c_str(), &endPtr);
    if (endPtr == n_str.c_str()) throw err("invalid N field");
    r.N_ref = N;

    recs.push_back(std::move(r));
  }
  return recs;
}

std::vector<MatchedMarkerInfo> matchMarkers(
    const PlinkData& plinkData,
    const std::vector<RefAfRecord>& refAf) {

  // Build ref lookup: key = "chr:pos:a1:a2"
  struct RefKey {
    std::string chrom;
    uint32_t    pos;
    std::string a1;
    std::string a2;
  };
  auto makeKey = [](const std::string& chr, uint32_t pos,
                    const std::string& a1, const std::string& a2) -> std::string {
    std::string k;
    k.reserve(chr.size() + 20 + a1.size() + a2.size());
    k += chr; k += ':';
    k += std::to_string(pos); k += ':';
    k += a1; k += ':'; k += a2;
    return k;
  };

  std::unordered_map<std::string, size_t> refMap;
  refMap.reserve(refAf.size());
  for (size_t i = 0; i < refAf.size(); ++i)
    refMap.emplace(makeKey(refAf[i].chrom, refAf[i].pos, refAf[i].a1, refAf[i].a2), i);

  // Match each bim marker against reference
  std::vector<MatchedMarkerInfo> matched;
  matched.reserve(plinkData.markerInfo().size());
  for (const auto& mi : plinkData.markerInfo()) {
    auto key = makeKey(mi.chrom, mi.pos, mi.ref, mi.alt);
    auto it = refMap.find(key);
    if (it == refMap.end()) continue;  // no match → drop
    const auto& ref = refAf[it->second];
    MatchedMarkerInfo m;
    m.genoIndex = mi.genoIndex;
    m.AF_ref    = ref.AF_ref;
    m.N_ref     = ref.N_ref;
    // mu0, mu1, n0, n1 will be filled later during genotype scanning
    m.mu0 = m.mu1 = m.n0 = m.n1 = m.mu_int = 0.0;
    matched.push_back(m);
  }
  return matched;
}


// ======================================================================
// Genotype scanning — compute per-marker case/control allele freq
// ======================================================================

void computeMarkerStats(
    std::vector<MatchedMarkerInfo>& matched,
    const PlinkData& plinkData,
    const ResidData& resid) {

  const uint32_t n = plinkData.nSubjUsed();
  PlinkCursor cursor(plinkData.bedFile(),
                     static_cast<uint32_t>(plinkData.nMarkers()),
                     plinkData.nSubjInFile(),
                     plinkData.samplePosMap(),
                     plinkData.isAltFirst(),
                     plinkData.isIdentityMap());

  if (!matched.empty())
    cursor.beginSequentialBlock(matched.front().genoIndex);

  Eigen::VectorXd gvec(n);

  for (auto& m : matched) {
    cursor.getGenotypesSimple(m.genoIndex, gvec);
    // Accumulate allele freq by case/control (indicator: 1=case, 0=control)
    // Missing genotypes are NaN (from getGenotypesSimple) and are skipped.
    double sum0 = 0.0, sum1 = 0.0;
    double cnt0 = 0.0, cnt1 = 0.0;
    for (uint32_t i = 0; i < n; ++i) {
      if (std::isnan(gvec[i])) continue;
      if (resid.indicator[i] == 1.0) {
        sum1 += gvec[i];
        cnt1 += 1.0;
      } else {
        sum0 += gvec[i];
        cnt0 += 1.0;
      }
    }
    m.mu0 = (cnt0 > 0.0) ? (sum0 / cnt0 / 2.0) : 0.0;  // control allele freq
    m.mu1 = (cnt1 > 0.0) ? (sum1 / cnt1 / 2.0) : 0.0;  // case allele freq
    m.n0  = cnt0;
    m.n1  = cnt1;
    m.mu_int = (cnt0 + cnt1 > 0.0)
             ? (sum0 + sum1) / (2.0 * (cnt0 + cnt1))
             : 0.0;
  }
}


// ======================================================================
// Phase 2 — Batch-effect testing
// ======================================================================

std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo>>
testBatchEffects(
    const std::vector<MatchedMarkerInfo>& matched,
    const ResidData& resid,
    const SparseGRM* grm,
    double refPrevalence,
    double cutoff) {

  const Eigen::Index nSubj = resid.residuals.size();

  // Precompute w1 and R_tilde
  double sumW = resid.weights.sum();
  Eigen::VectorXd w1 = resid.weights / (2.0 * sumW);
  double meanR = resid.residuals.mean();
  Eigen::VectorXd R_tilde = resid.residuals.array() - meanR;

  // --- Variance ratios from sparse GRM ---
  double grm_sum_cov_w   = 0.0;    // sum(GRM_ij * w1_i * w1_j)
  double grm_sum_cov_R   = 0.0;    // sum(GRM_ij * R_tilde_i * R_tilde_j)
  bool   hasGRM = (grm != nullptr);
  if (hasGRM) {
    grm_sum_cov_w = grm->quadForm(w1.data(), static_cast<uint32_t>(nSubj));
    grm_sum_cov_R = grm->quadForm(R_tilde.data(), static_cast<uint32_t>(nSubj));
  }
  double sum_w1_sq    = w1.array().square().sum();
  double sum_Rtilde_sq = R_tilde.array().square().sum();
  double var_ratio_int = hasGRM ? (grm_sum_cov_R / sum_Rtilde_sq) : 1.0;

  // --- Per-marker batch p-value ---
  const size_t nMarkers = matched.size();
  auto refInfoMap = std::make_shared<std::unordered_map<uint64_t, WtCoxGRefInfo>>();
  refInfoMap->reserve(nMarkers);

  // Temporary per-marker storage
  struct MarkerBatchData {
    uint64_t genoIndex;
    double   AF_ref, N_ref;
    double   mu0, mu1, n0, n1, mu_int;
    double   var_ratio_w0;
    double   pvalue_bat;
  };
  std::vector<MarkerBatchData> batchData(nMarkers);

  for (size_t i = 0; i < nMarkers; ++i) {
    const auto& m = matched[i];
    auto& bd = batchData[i];
    bd.genoIndex = m.genoIndex;
    bd.AF_ref = m.AF_ref;  bd.N_ref = m.N_ref;
    bd.mu0 = m.mu0;  bd.mu1 = m.mu1;
    bd.n0 = m.n0;    bd.n1 = m.n1;
    bd.mu_int = m.mu_int;

    double n_ext = m.N_ref;

    // var_ratio_w0 per marker (depends on N_ref)
    bd.var_ratio_w0 = hasGRM
        ? (grm_sum_cov_w + 1.0 / (2.0 * n_ext)) / (sum_w1_sq + 1.0 / (2.0 * n_ext))
        : 1.0;

    // Batch effect p-value (Batcheffect.TestOneMarker)
    double er = m.n1 / (m.n1 + m.n0);
    double w0_val = (1.0 - refPrevalence) / refPrevalence / ((1.0 - er) / er);
    double w1_val = 1.0;
    double weight_maf = (m.mu0 * w0_val * m.n0 + m.mu1 * w1_val * m.n1)
                      / (w0_val * m.n0 + w1_val * m.n1);
    double est_maf = (m.mu0 * w0_val * m.n0 + m.mu1 * w1_val * m.n1 + m.AF_ref * n_ext * w0_val)
                   / (m.n1 * w1_val + m.n0 * w0_val + n_ext * w0_val);
    double v = ((m.n1 * w1_val * w1_val + m.n0 * w0_val * w0_val)
               / (2.0 * std::pow(m.n1 * w1_val + m.n0 * w0_val, 2.0))
               + 1.0 / (2.0 * n_ext))
             * est_maf * (1.0 - est_maf);
    double z = (v > 0.0) ? (weight_maf - m.AF_ref) / std::sqrt(v) : 0.0;
    double z_adj = (bd.var_ratio_w0 > 0.0) ? z / std::sqrt(bd.var_ratio_w0) : z;
    bd.pvalue_bat = 2.0 * math::pnorm(-std::abs(z_adj));
  }

  // --- Estimate TPR, sigma2, w.ext per MAF group ---
  // MAF groups: [−1e-4, 0.05), [0.05, 0.10), … , [0.35, 0.40), [0.40, max_mu_int]
  double max_mu_int = 0.0;
  for (auto& bd : batchData)
    max_mu_int = std::max(max_mu_int, bd.mu_int);

  std::vector<double> mafBreaks;
  for (double x = -1e-4; x <= 0.40 + 1e-9; x += 0.05)
    mafBreaks.push_back(x);
  mafBreaks.push_back(max_mu_int + 1e-6);

  for (size_t grp = 0; grp + 1 < mafBreaks.size(); ++grp) {
    double lo = mafBreaks[grp];
    double hi = mafBreaks[grp + 1];
    double mu = (lo + hi) / 2.0;

    // Collect markers in this group
    std::vector<size_t> idx1;  // narrow group: mu_int in (lo, hi]
    for (size_t i = 0; i < nMarkers; ++i)
      if (batchData[i].mu_int > lo && batchData[i].mu_int <= hi)
        idx1.push_back(i);
    if (idx1.empty()) continue;

    // Wider group for parameter estimation: mu_int in [lo-0.1, hi+0.1]
    std::vector<double> vec_p_bat;
    for (size_t i = 0; i < nMarkers; ++i)
      if (batchData[i].mu_int >= std::max(lo - 0.1, 0.0)
          && batchData[i].mu_int < std::min(1.0, hi + 0.1)
          && !std::isnan(batchData[i].pvalue_bat))
        vec_p_bat.push_back(batchData[i].pvalue_bat);

    if (vec_p_bat.empty()) continue;

    double n_ext = batchData[idx1[0]].N_ref;
    if (std::isnan(n_ext) || n_ext <= 0.0) continue;
    double var_mu_ext = mu * (1.0 - mu) / (2.0 * n_ext);

    // var_Sbat for this MAF group
    double vr_w0 = hasGRM ? batchData[idx1[0]].var_ratio_w0 : 1.0;
    double var_Sbat = hasGRM
        ? vr_w0 * (sum_w1_sq * 2.0 * mu * (1.0 - mu) + var_mu_ext)
        : sum_w1_sq * 2.0 * mu * (1.0 - mu) + var_mu_ext;

    // Empirical pass rates at several cutoffs
    static constexpr double vec_cutoff[] = {0.01, 0.11, 0.21, 0.31};
    static constexpr int    nCut = 4;
    double vec_p_deno[nCut];
    for (int j = 0; j < nCut; ++j) {
      int cnt = 0;
      for (double p : vec_p_bat) if (p > vec_cutoff[j]) ++cnt;
      vec_p_deno[j] = static_cast<double>(cnt) / static_cast<double>(vec_p_bat.size());
    }

    // Nelder-Mead to estimate [TPR, sigma2]
    auto opti_fun = [&](const std::vector<double>& par) -> double {
      double diff = 0.0;
      for (int j = 0; j < nCut; ++j) {
        double p_cut = vec_cutoff[j];
        double q = math::qnorm(1.0 - p_cut / 2.0);
        double lb = -q * std::sqrt(var_Sbat);
        double ub =  q * std::sqrt(var_Sbat);
        double var_Sbat_par2 = var_Sbat + par[1];
        double c_val, d_val;
        if (var_Sbat_par2 >= 0) {
          c_val = math::pnorm(ub, 0.0, std::sqrt(var_Sbat_par2), true, true);
          d_val = math::pnorm(lb, 0.0, std::sqrt(var_Sbat_par2), true, true);
        } else {
          return 1e30;
        }
        double pro_cut = par[0] * (std::exp(d_val) * (std::exp(c_val - d_val) - 1.0))
                       + (1.0 - par[0]) * (1.0 - p_cut);
        double ratio = (vec_p_deno[j] - pro_cut) / (vec_p_deno[j] + 1e-300);
        diff += ratio * ratio;
      }
      return diff;
    };

    auto optResult = math::nelderMead(opti_fun, {0.01, 0.01});
    double TPR    = std::clamp(optResult.par[0], 0.0, 1.0);
    double sigma2 = std::clamp(optResult.par[1], 0.0, 1.0);

    // Optimal external weight via 1-D Brent minimisation
    // The R code does a complex nested optimisation (optim + uniroot + pmvnorm).
    // We replicate the same logic: for a given b, find mu1 via root-finding,
    // then compute the power metric.
    auto fun_optimalWeight = [&](double b) -> double {
      auto p_fun = [&](double mu1_trial) -> double {
        double mu0_val = mu;
        double mu_pop = mu1_trial * refPrevalence + mu0_val * (1.0 - refPrevalence);
        // Build per-subject mu_i
        const auto& R = resid.residuals;
        const auto& y = resid.indicator;
        double nS = static_cast<double>(R.size());
        double meanR_loc = R.mean();
        double sumR_loc  = R.sum();

        // S = sum((R - (1-b)*meanR) * mu_i) - sumR * 2 * b * mu_pop
        double S = 0.0;
        for (Eigen::Index k = 0; k < R.size(); ++k) {
          double mu_i = (y[k] == 1.0) ? 2.0 * mu1_trial : 2.0 * mu0_val;
          S += (R[k] - (1.0 - b) * meanR_loc) * mu_i;
        }
        S -= sumR_loc * 2.0 * b * mu_pop;

        double w1sum2 = sum_w1_sq;  // already computed above
        double mu_local = 0.0;
        for (Eigen::Index k = 0; k < R.size(); ++k)
          mu_local += ((y[k] == 1.0) ? 2.0 * mu1_trial : 2.0 * mu0_val);
        mu_local /= (2.0 * nS);

        double var_mu_ext_loc = mu_local * (1.0 - mu_local) / (2.0 * n_ext);
        double var_Sbat_loc = w1sum2 * 2.0 * mu_local * (1.0 - mu_local) + var_mu_ext_loc;

        double p_cut = 0.1;
        double q = math::qnorm(1.0 - p_cut / 2.0);
        double lb = -q * std::sqrt(var_Sbat_loc);
        double ub =  q * std::sqrt(var_Sbat_loc);
        double c_val = math::pnorm(ub, 0.0, std::sqrt(var_Sbat_loc + sigma2), true, true);
        double d_val = math::pnorm(lb, 0.0, std::sqrt(var_Sbat_loc + sigma2), true, true);
        double p_deno = TPR * (std::exp(d_val) * (std::exp(c_val - d_val) - 1.0))
                      + (1.0 - TPR) * (1.0 - p_cut);

        double var_int = 0.0;
        for (Eigen::Index k = 0; k < R.size(); ++k) {
          double r_adj = R[k] - (1.0 - b) * meanR_loc;
          var_int += r_adj * r_adj;
        }
        var_int *= 2.0 * mu_local * (1.0 - mu_local);
        double var_S = var_int + 4.0 * b * b * sumR_loc * sumR_loc * var_mu_ext_loc;

        double cov_val = 0.0;
        for (Eigen::Index k = 0; k < R.size(); ++k)
          cov_val += w1[k] * (R[k] - (1.0 - b) * meanR_loc);
        cov_val *= 2.0 * mu_local * (1.0 - mu_local);
        cov_val += 2.0 * b * sumR_loc * var_mu_ext_loc;

        // p0 = P(S ≤ −|S|, lb ≤ S_bat ≤ ub) via bivariate normal
        double negInf = -std::numeric_limits<double>::infinity();
        double p0 = math::pmvnorm2d(negInf, -std::abs(S), lb, ub, var_S, cov_val, var_Sbat_loc);
        p0 = std::clamp(p0, 0.0, 1.0);

        // p1 with sigma2
        double var_S1 = var_int + 4.0 * b * b * sumR_loc * sumR_loc * (var_mu_ext_loc + sigma2);
        double cov_val1 = 0.0;
        for (Eigen::Index k = 0; k < R.size(); ++k)
          cov_val1 += w1[k] * (R[k] - (1.0 - b) * meanR_loc);
        cov_val1 *= 2.0 * mu_local * (1.0 - mu_local);
        cov_val1 += 2.0 * b * sumR_loc * (var_mu_ext_loc + sigma2);
        double var_Sbat1 = var_Sbat_loc + sigma2;

        double p1 = math::pmvnorm2d(negInf, -std::abs(S), lb, ub, var_S1, cov_val1, var_Sbat1);
        p1 = std::clamp(p1, 0.0, 1.0);

        double p_con = 2.0 * (TPR * p1 + (1.0 - TPR) * p0) / (p_deno + 1e-300);
        return -std::log10(p_con / 5e-8 + 1e-300);
      };

      // Find mu1 such that p_fun(mu1) == 0 via Brent root finding
      double mu1;
      try {
        mu1 = math::findRootBrent(p_fun, mu, 1.0 - 1e-6, 1e-6);
      } catch (...) {
        mu1 = mu;
      }
      return mu1;  // optim minimises this → lower mu1 = more power
    };

    double w_ext = math::brentMin(fun_optimalWeight, 0.0, 1.0, 1e-6, 200);

    // Compute var_ratio_ext from GRM (if available)
    double var_ratio_ext = 1.0;
    if (hasGRM) {
      Eigen::VectorXd R_tilde_w = resid.residuals.array() - meanR * w_ext;
      double grm_cov_Rext = grm->quadForm(R_tilde_w.data(), static_cast<uint32_t>(nSubj));
      double sumR_sq_over_n = w_ext * w_ext * resid.residuals.sum() * resid.residuals.sum() / n_ext;
      double num = grm_cov_Rext + sumR_sq_over_n;
      double den = R_tilde_w.array().square().sum() + sumR_sq_over_n;
      var_ratio_ext = (den > 0.0) ? num / den : 1.0;
    }

    // Populate refInfoMap for markers in this group
    for (size_t i : idx1) {
      WtCoxGRefInfo ri;
      ri.AF_ref       = batchData[i].AF_ref;
      ri.N_ref        = batchData[i].N_ref;
      ri.TPR          = TPR;
      ri.sigma2       = sigma2;
      ri.pvalue_bat   = batchData[i].pvalue_bat;
      ri.w_ext        = w_ext;
      ri.var_ratio_w0 = batchData[i].var_ratio_w0;
      ri.var_ratio_int = var_ratio_int;
      ri.var_ratio_ext = var_ratio_ext;
      refInfoMap->emplace(batchData[i].genoIndex, ri);
    }
  }

  return refInfoMap;
}


// ======================================================================
// Phase 3 — WtCoxGMethod (MethodBase implementation)
// ======================================================================

WtCoxGMethod::WtCoxGMethod(
    Eigen::VectorXd R, Eigen::VectorXd w,
    double cutoff, double SPA_Cutoff,
    std::shared_ptr<const std::unordered_map<uint64_t, WtCoxGRefInfo>> refMap)
  : m_R(std::move(R)), m_w(std::move(w)),
    m_w1(m_w / (2.0 * m_w.sum())),
    m_meanR(m_R.mean()), m_sumR(m_R.sum()),
    m_cutoff(cutoff), m_SPA_Cutoff(SPA_Cutoff),
    m_refMap(std::move(refMap)) {}

std::unique_ptr<MethodBase> WtCoxGMethod::clone() const {
  auto p = std::make_unique<WtCoxGMethod>(m_R, m_w, m_cutoff, m_SPA_Cutoff, m_refMap);
  return p;
}

std::string WtCoxGMethod::getHeaderColumns() const {
  return "\tWtCoxG.ext\tWtCoxG.noext\tzscore.ext\tzscore.noext\tAF_ref\tN_ref\tpvalue_bat";
}

void WtCoxGMethod::prepareChunk(const std::vector<uint64_t>& gIndices) {
  size_t n = gIndices.size();
  m_chunkRefInfo.resize(n);
  for (size_t i = 0; i < n; ++i) {
    WtCoxGRefInfo ri;  // defaults to NaN
    if (m_refMap) {
      auto it = m_refMap->find(gIndices[i]);
      if (it != m_refMap->end()) ri = it->second;
    }
    m_chunkRefInfo[i] = ri;
  }
}

void WtCoxGMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double /*altFreq*/, int markerInChunkIdx,
    bool /*flipped*/, std::vector<double>& result) {

  const auto& info = m_chunkRefInfo[markerInChunkIdx];

  // With external reference
  WtResult res_ext = wtCoxGTest(
      GVec,
      info.pvalue_bat, info.TPR, info.sigma2, info.w_ext,
      info.var_ratio_int, info.var_ratio_w0, info.var_ratio_w0,
      info.var_ratio_ext, info.var_ratio_ext,
      info.AF_ref, info.N_ref, m_cutoff);

  // Without external reference
  WtResult res_noext = wtCoxGTest(
      GVec,
      info.pvalue_bat, NaN::quiet_NaN(), NaN::quiet_NaN(), 0.0,
      info.var_ratio_int, 1.0, 1.0, 1.0, 1.0,
      NaN::quiet_NaN(), NaN::quiet_NaN(), m_cutoff);

  result.push_back(res_ext.pval);
  result.push_back(res_noext.pval);
  result.push_back(res_ext.zscore);
  result.push_back(res_noext.zscore);
  result.push_back(info.AF_ref);
  result.push_back(info.N_ref);
  result.push_back(info.pvalue_bat);
}

WtCoxGMethod::DualResult WtCoxGMethod::computeDual(
    Eigen::Ref<Eigen::VectorXd> GVec, int markerInChunkIdx) {

  const auto& info = m_chunkRefInfo[markerInChunkIdx];

  WtResult res_ext = wtCoxGTest(
      GVec,
      info.pvalue_bat, info.TPR, info.sigma2, info.w_ext,
      info.var_ratio_int, info.var_ratio_w0, info.var_ratio_w0,
      info.var_ratio_ext, info.var_ratio_ext,
      info.AF_ref, info.N_ref, m_cutoff);

  WtResult res_noext = wtCoxGTest(
      GVec,
      info.pvalue_bat, NaN::quiet_NaN(), NaN::quiet_NaN(), 0.0,
      info.var_ratio_int, 1.0, 1.0, 1.0, 1.0,
      NaN::quiet_NaN(), NaN::quiet_NaN(), m_cutoff);

  return {res_ext.pval, res_noext.pval, res_ext.score, res_noext.score};
}

WtCoxGMethod::WtResult WtCoxGMethod::wtCoxGTest(
    const Eigen::Ref<const Eigen::VectorXd>& g_input,
    double p_bat, double TPR, double sigma2, double b,
    double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
    double var_ratio0, double var_ratio1,
    double mu_ext, double n_ext, double p_cut) const {

  // No external info → delegate to SPA-only
  if (std::isnan(mu_ext)) {
    double vr = (std::isnan(TPR) && std::isnan(sigma2)) ? var_ratio_int : 1.0;
    SpaResult spa = spaGOneSnpHomo(g_input, m_R, 0.0, 0.0, 0.0, 0.0, vr, m_SPA_Cutoff);
    return {spa.pval, spa.score, spa.zscore};
  }

  double sum_g = g_input.sum();
  double sum_2mg = (2.0 - g_input.array()).sum();
  if (p_bat < p_cut || std::isnan(p_bat) || sum_g < 10 || sum_2mg < 10)
    return {NaN::quiet_NaN(), NaN::quiet_NaN(), NaN::quiet_NaN()};

  double mu_int = g_input.mean() / 2.0;
  double mu = (1.0 - b) * mu_int + b * mu_ext;
  double S = (m_R.array() * (g_input.array() - 2.0 * mu)).sum();

  double var_mu_ext = mu * (1.0 - mu) / (2.0 * n_ext);
  double var_Sbat = (m_w1.array().square()).sum() * 2.0 * mu * (1.0 - mu) + var_mu_ext;

  double qnorm_val = math::qnorm(1.0 - p_cut / 2.0);
  double lb = -qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);
  double ub =  qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);

  double c_val = math::pnorm(ub / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
  double d_val = math::pnorm(lb / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
  double p_deno = TPR * (std::exp(d_val) * (std::exp(c_val - d_val) - 1.0))
                + (1.0 - TPR) * (1.0 - p_cut);

  // Internal SPA (no sigma2)
  SpaResult spa_s0 = spaGOneSnpHomo(g_input, m_R, mu_ext, n_ext, b, 0.0, var_ratio0, m_SPA_Cutoff);
  double qchi = math::qchisq(spa_s0.pval, 1.0, false, false);
  double var_S = (qchi > 0.0) ? S * S / var_ratio0 / qchi : NaN::quiet_NaN();

  // Covariance between S_bat and S
  Eigen::ArrayXd R_adj = m_R.array() - (1.0 - b) * m_meanR;
  double var_int_denom = (R_adj * R_adj).sum() * 2.0 * mu * (1.0 - mu)
                       + 4.0 * b * b * m_sumR * m_sumR * var_mu_ext;
  if (var_int_denom <= 0.0)
    return {NaN::quiet_NaN(), NaN::quiet_NaN(), NaN::quiet_NaN()};

  double cov_val = (m_w1.array() * R_adj).sum() * 2.0 * mu * (1.0 - mu)
                 + 2.0 * b * m_sumR * var_mu_ext;
  cov_val *= std::sqrt(var_S / var_int_denom);
  double z = S / std::sqrt(var_S);

  double negInf = -std::numeric_limits<double>::infinity();
  double p0 = math::pmvnorm2d(
      negInf, -std::abs(S / std::sqrt(var_ratio0)),
      lb / std::sqrt(var_ratio_w0), ub / std::sqrt(var_ratio_w0),
      var_S, cov_val, var_Sbat);
  p0 = std::clamp(p0, 0.0, 1.0);

  // External SPA (with sigma2)
  SpaResult spa_s1 = spaGOneSnpHomo(g_input, m_R, mu_ext, n_ext, b, sigma2, var_ratio1, m_SPA_Cutoff);
  double var_S1 = S * S / var_ratio1 / math::qchisq(spa_s1.pval, 1.0, false, false);
  double cov_val1 = (m_w1.array() * R_adj).sum() * 2.0 * mu * (1.0 - mu)
                  + 2.0 * b * m_sumR * (var_mu_ext + sigma2);
  cov_val1 *= std::sqrt(var_S1 / var_int_denom);
  double var_Sbat1 = var_Sbat + sigma2;

  double p1 = math::pmvnorm2d(
      negInf, -std::abs(S / std::sqrt(var_ratio1)),
      lb / std::sqrt(var_ratio_w1), ub / std::sqrt(var_ratio_w1),
      var_S1, cov_val1, var_Sbat1);
  p1 = std::clamp(p1, 0.0, 1.0);

  double p_con = 2.0 * (TPR * p1 + (1.0 - TPR) * p0) / p_deno;
  return {p_con, S, z};
}


// ======================================================================
// Top-level orchestration
// ======================================================================

void runWtCoxG(
    const std::string& residFile,
    const std::string& bfilePrefix,
    const std::string& refAfFile,
    const std::string& sparseGrmFile,
    const std::string& outputFile,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff) {

  // ---- Load & match ----
  infoMsg("Loading resid file: %s", residFile.c_str());
  ResidData resid = loadResidFile(residFile);
  infoMsg("  %zu subjects loaded", resid.subjects.size());

  infoMsg("Loading ref-af file: %s", refAfFile.c_str());
  auto refAf = loadRefAfFile(refAfFile);
  infoMsg("  %zu reference records loaded", refAf.size());

  infoMsg("Loading PLINK data: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      resid.subjects,
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);
  infoMsg("  %u subjects matched in PLINK, %u markers available",
          plinkData.nSubjUsed(), plinkData.nMarkers());

  infoMsg("Matching markers against reference allele frequencies...");
  auto matched = matchMarkers(plinkData, refAf);
  infoMsg("  %zu markers matched", matched.size());

  infoMsg("Computing per-marker case/control allele frequencies...");
  computeMarkerStats(matched, plinkData, resid);

  // ---- Batch-effect testing ----
  infoMsg("Batch-effect testing and parameter estimation...");
  std::unique_ptr<SparseGRM> grm;
  if (!sparseGrmFile.empty()) {
    infoMsg("  Loading sparse GRM: %s", sparseGrmFile.c_str());
    grm = std::make_unique<SparseGRM>(sparseGrmFile, resid.subjects);
    infoMsg("  Sparse GRM: %u subjects, %zu non-zeros",
            grm->nSubjects(), grm->nnz());
  }
  auto refInfoMap = testBatchEffects(
      matched, resid, grm.get(), refPrevalence, cutoff);
  infoMsg("  %zu markers retained after batch-effect QC", refInfoMap->size());

  // ---- Marker-level SPA tests ----
  infoMsg("Running marker-level WtCoxG tests (%d thread(s))...", nthread);
  WtCoxGMethod method(
      resid.residuals, resid.weights,
      cutoff, spaCutoff, refInfoMap);
  markerEngine(plinkData, method, outputFile,
               nthread,
               missingCutoff,
               minMafCutoff,
               minMacCutoff,
               /*exactHwe=*/false);
}
