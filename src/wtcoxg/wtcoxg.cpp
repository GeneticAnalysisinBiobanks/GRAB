// wtcoxg.cpp — WtCoxG full implementation (Phases 1–3)

#include "wtcoxg/wtcoxg.hpp"
#include "io/plink.hpp"
#include "io/subject_data.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/text_scanner.hpp"

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
            double obs_ct, double N_all, double sumR,
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
             double obs_ct, double N_all, double sumR,
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
          double obs_ct, double N_all, double sumR,
          double var_mu_ext, double /*g_var_est*/,
          double meanR, double b) {
  double var_adj = obs_ct * (sumR / N_all) * (sumR / N_all) * MAF * (1.0 - MAF);
  double result = 0.0;
  const double bm = (1.0 - b) * meanR;
  for (Eigen::Index i = 0; i < R.size(); ++i) {
    double R_adj = R[i] - bm;
    result += R_adj * R_adj * math::kG2(t * R_adj, MAF);
  }
  return result + var_adj;
}

double getProbSpaG(double MAF, const Eigen::VectorXd& R, double s,
                   double obs_ct, double N_all, double sumR,
                   double var_mu_ext, double g_var_est,
                   double meanR, double b, bool lower_tail) {
  auto h1_func = [&](double t) {
    return h1Adj(t, R, s, MAF, obs_ct, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
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

  double k1 = hOrg(zeta, R, MAF, obs_ct, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
  double k2 = h2(zeta, R, MAF, obs_ct, N_all, sumR, var_mu_ext, g_var_est, meanR, b);
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
    double mu_ext, double obs_ct, double b,
    double sigma2, double var_ratio, double SPA_Cutoff) {

  if (std::isnan(mu_ext)) { mu_ext = 0.0; obs_ct = 0.0; }

  const double N = static_cast<double>(g.size());
  double mu_int = g.mean() / 2.0;
  double MAF = std::clamp((1.0 - b) * mu_int + b * mu_ext, 0.0, 1.0);
  double sumR = R.sum();
  double N_all = N + obs_ct / 2.0;
  double S = (R.array() * (g.array() - 2.0 * MAF)).sum();
  double S_raw = S;
  S /= var_ratio;

  double g_var_est = 2.0 * MAF * (1.0 - MAF);
  double var_mu_ext = (obs_ct == 0.0) ? 0.0 : (MAF * (1.0 - MAF) / obs_ct + sigma2);

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

  double pval1 = getProbSpaG(MAF, R,  std::abs(S), obs_ct, N_all, sumR, var_mu_ext, g_var_est, meanR, b, false);
  double pval2 = getProbSpaG(MAF, R, -std::abs(S), obs_ct, N_all, sumR, var_mu_ext, g_var_est, meanR, b, true);
  double pval_spa = std::min(1.0, pval1 + pval2);
  return {pval_spa, std::min(1.0, 2.0 * math::pnorm(-std::abs(z))), S_raw, z};
}

} // anon namespace


// ======================================================================
// Phase 1 — File parsing & marker matching
// ======================================================================

std::vector<RefAfRecord> loadRefAfFile(const std::string& filename,
                                       bool* isNumericFallback) {
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open ref-af file: " + filename);

  std::vector<RefAfRecord> recs;
  recs.reserve(100000);
  std::string line;

  // ---- Detect column positions from header line ----
  int colChrom = -1, colId = -1, colRef = -1, colAlt = -1;
  int colAltFreqs = -1, colObsCt = -1;

  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line[0] != '#') break;  // first non-header line
    // Parse header columns
    std::istringstream hss(line);
    std::string tok;
    int col = 0;
    while (hss >> tok) {
      // Strip leading '#' from the first token
      if (col == 0 && !tok.empty() && tok[0] == '#')
        tok = tok.substr(1);
      if (tok == "CHROM")              colChrom    = col;
      else if (tok == "ID")            colId       = col;
      else if (tok == "REF")           colRef      = col;
      else if (tok == "ALT" || tok == "ALT1")
                                       colAlt      = col;
      else if (tok == "ALT_FREQS" || tok == "ALT1_FREQ")
                                       colAltFreqs = col;
      else if (tok == "OBS_CT")        colObsCt    = col;
      ++col;
    }
    // Keep reading — last '#' line wins (plink2 puts one header line)
  }

  if (colChrom < 0 || colId < 0 || colRef < 0 || colAlt < 0 ||
      colAltFreqs < 0 || colObsCt < 0) {
    // ---- Two-column numeric fallback ----
    // No valid header found.  If the first data line has exactly two
    // numeric tokens, treat the whole file as (ALT_FREQS  OBS_CT) rows
    // assumed to be in .bim order.
    if (!line.empty() && line[0] != '#') {
      char* end1 = nullptr;
      char* end2 = nullptr;
      std::istringstream probe(line);
      std::string t1, t2, t3;
      if ((probe >> t1 >> t2) && !(probe >> t3)) {
        std::strtod(t1.c_str(), &end1);
        std::strtod(t2.c_str(), &end2);
        if (end1 != t1.c_str() && end2 != t2.c_str()) {
          // Confirmed two-column numeric format
          if (isNumericFallback) *isNumericFallback = true;
          uint32_t lineNo = 0;
          auto parseNumLine = [&](const std::string& ln) {
            ++lineNo;
            if (ln.empty()) return;
            std::istringstream iss(ln);
            RefAfRecord r;
            if (!(iss >> r.alt_freq >> r.obs_ct))
              throw std::runtime_error(
                  filename + " line " + std::to_string(lineNo) +
                  ": expected 2 numeric columns (ALT_FREQS OBS_CT)");
            recs.push_back(std::move(r));
          };
          parseNumLine(line);
          while (std::getline(ifs, line)) {
            if (text::skipLine(line)) continue;
            parseNumLine(line);
          }
          return recs;
        }
      }
    }
    throw std::runtime_error(
        filename + ": missing required header columns "
        "(need #CHROM, ID, REF, ALT, ALT_FREQS, OBS_CT)");
  }

  if (isNumericFallback) *isNumericFallback = false;

  const int maxCol = std::max({colChrom, colId, colRef, colAlt,
                               colAltFreqs, colObsCt});

  // ---- Parse data lines ----
  // `line` already holds the first non-header line from the header scan
  uint32_t lineNo = 0;
  auto parseLine = [&](const std::string& ln) {
    ++lineNo;
    if (ln.empty()) return;
    // Tokenise
    std::vector<std::string> tokens;
    tokens.reserve(maxCol + 2);
    std::istringstream iss(ln);
    std::string t;
    while (iss >> t) tokens.push_back(std::move(t));
    if (static_cast<int>(tokens.size()) <= maxCol)
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
          ": expected at least " + std::to_string(maxCol + 1) +
          " columns, got " + std::to_string(tokens.size()));

    RefAfRecord r;
    r.chrom      = std::move(tokens[colChrom]);
    r.id         = std::move(tokens[colId]);
    r.ref_allele = std::move(tokens[colRef]);
    r.alt_allele = std::move(tokens[colAlt]);
    // Uppercase alleles for consistent matching
    for (auto& ch : r.ref_allele) ch = static_cast<char>(std::toupper(ch));
    for (auto& ch : r.alt_allele) ch = static_cast<char>(std::toupper(ch));

    char* endPtr;
    r.alt_freq = std::strtod(tokens[colAltFreqs].c_str(), &endPtr);
    if (endPtr == tokens[colAltFreqs].c_str())
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
          ": invalid ALT_FREQS value");
    r.obs_ct = std::strtod(tokens[colObsCt].c_str(), &endPtr);
    if (endPtr == tokens[colObsCt].c_str())
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
          ": invalid OBS_CT value");
    recs.push_back(std::move(r));
  };

  // Process the first data line that was read during header scan
  if (!line.empty() && line[0] != '#')
    parseLine(line);
  while (std::getline(ifs, line)) {
    if (text::skipLine(line)) continue;
    parseLine(line);
  }
  return recs;
}

std::vector<MatchedMarkerInfo> matchMarkers(
    const PlinkData& plinkData,
    const std::vector<RefAfRecord>& refAf) {

  // Build ref lookup: key = "chrom:id" → index into refAf
  auto makeKey = [](const std::string& chr, const std::string& id) -> std::string {
    std::string k;
    k.reserve(chr.size() + 1 + id.size());
    k += chr; k += ':'; k += id;
    return k;
  };

  std::unordered_map<std::string, size_t> refMap;
  refMap.reserve(refAf.size());
  for (size_t i = 0; i < refAf.size(); ++i)
    refMap.emplace(makeKey(refAf[i].chrom, refAf[i].id), i);

  // Match each bim marker against reference by (CHROM, ID),
  // then check allele orientation.
  // With ref-first allele order: mi.ref = bim col5, mi.alt = bim col6.
  std::vector<MatchedMarkerInfo> matched;
  matched.reserve(plinkData.markerInfo().size());
  for (const auto& mi : plinkData.markerInfo()) {
    auto key = makeKey(mi.chrom, mi.id);
    auto it = refMap.find(key);
    if (it == refMap.end()) continue;  // no match by CHROM+ID

    const auto& ref = refAf[it->second];
    MatchedMarkerInfo m;
    m.genoIndex = mi.genoIndex;

    //  Case 1: afreq ALT == bim col5 AND afreq REF == bim col6
    //          → same orientation, AF_ref = ALT_FREQS
    //  Case 2: afreq REF == bim col5 AND afreq ALT == bim col6
    //          → flipped, AF_ref = 1 - ALT_FREQS
    if (ref.alt_allele == mi.ref && ref.ref_allele == mi.alt) {
      m.AF_ref = ref.alt_freq;
    } else if (ref.ref_allele == mi.ref && ref.alt_allele == mi.alt) {
      m.AF_ref = 1.0 - ref.alt_freq;
    } else {
      continue;  // alleles don't match → drop
    }
    m.obs_ct = ref.obs_ct;

    // mu0, mu1, n0, n1 will be filled later during genotype scanning
    m.mu0 = m.mu1 = m.n0 = m.n1 = m.mu_int = 0.0;
    matched.push_back(m);
  }
  return matched;
}

std::vector<MatchedMarkerInfo> matchMarkersNumeric(
    const PlinkData& plinkData,
    const std::vector<RefAfRecord>& refAf) {

  const auto& markers = plinkData.markerInfo();
  if (refAf.size() != markers.size())
    throw std::runtime_error(
        "ref-af numeric fallback: row count (" +
        std::to_string(refAf.size()) + ") != bim marker count (" +
        std::to_string(markers.size()) + ")");

  std::vector<MatchedMarkerInfo> matched;
  matched.reserve(markers.size());
  for (size_t i = 0; i < markers.size(); ++i) {
    MatchedMarkerInfo m;
    m.genoIndex = markers[i].genoIndex;
    m.AF_ref    = refAf[i].alt_freq;       // col 5 frequency directly
    m.obs_ct    = refAf[i].obs_ct;
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
    const Eigen::VectorXd& indicator) {

  const uint32_t n = plinkData.nSubjUsed();
  PlinkCursor cursor(plinkData.bedFile(),
                     static_cast<uint32_t>(plinkData.nMarkers()),
                     plinkData.nSubjInFile(),
                     plinkData.usedMask(),
                     plinkData.nSubjUsed(),
                     plinkData.isAltFirst(),
                     plinkData.allUsed());

  if (!matched.empty())
    cursor.beginSequentialBlock(matched.front().genoIndex);

  Eigen::VectorXd gvec(n);

  for (auto& m : matched) {
    cursor.getGenotypesSimple(m.genoIndex, gvec);
    double sum0 = 0.0, sum1 = 0.0;
    double cnt0 = 0.0, cnt1 = 0.0;
    for (uint32_t i = 0; i < n; ++i) {
      if (std::isnan(gvec[i])) continue;
      if (indicator[i] == 1.0) {
        sum1 += gvec[i];
        cnt1 += 1.0;
      } else {
        sum0 += gvec[i];
        cnt0 += 1.0;
      }
    }
    m.mu0 = (cnt0 > 0.0) ? (sum0 / cnt0 / 2.0) : 0.0;
    m.mu1 = (cnt1 > 0.0) ? (sum1 / cnt1 / 2.0) : 0.0;
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
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& weights,
    const Eigen::VectorXd& indicator,
    const SparseGRM* grm,
    double refPrevalence,
    double cutoff) {

  const Eigen::Index nSubj = residuals.size();

  // Precompute w1 and R_tilde
  double sumW = weights.sum();
  Eigen::VectorXd w1 = weights / (2.0 * sumW);
  double meanR = residuals.mean();
  Eigen::VectorXd R_tilde = residuals.array() - meanR;

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
    double   AF_ref, obs_ct;
    double   mu0, mu1, n0, n1, mu_int;
    double   var_ratio_w0;
    double   pvalue_bat;
  };
  std::vector<MarkerBatchData> batchData(nMarkers);

  for (size_t i = 0; i < nMarkers; ++i) {
    const auto& m = matched[i];
    auto& bd = batchData[i];
    bd.genoIndex = m.genoIndex;
    bd.AF_ref = m.AF_ref;  bd.obs_ct = m.obs_ct;
    bd.mu0 = m.mu0;  bd.mu1 = m.mu1;
    bd.n0 = m.n0;    bd.n1 = m.n1;
    bd.mu_int = m.mu_int;

    // var_ratio_w0 per marker (depends on obs_ct)
    bd.var_ratio_w0 = hasGRM
        ? (grm_sum_cov_w + 1.0 / m.obs_ct) / (sum_w1_sq + 1.0 / m.obs_ct)
        : 1.0;

    // Batch effect p-value (Batcheffect.TestOneMarker)
    double er = m.n1 / (m.n1 + m.n0);
    double w0_val = (1.0 - refPrevalence) / refPrevalence / ((1.0 - er) / er);
    double w1_val = 1.0;
    double weight_maf = (m.mu0 * w0_val * m.n0 + m.mu1 * w1_val * m.n1)
                      / (w0_val * m.n0 + w1_val * m.n1);
    double est_maf = (m.mu0 * w0_val * m.n0 + m.mu1 * w1_val * m.n1 + m.AF_ref * (m.obs_ct / 2.0) * w0_val)
                   / (m.n1 * w1_val + m.n0 * w0_val + (m.obs_ct / 2.0) * w0_val);
    double v = ((m.n1 * w1_val * w1_val + m.n0 * w0_val * w0_val)
               / (2.0 * std::pow(m.n1 * w1_val + m.n0 * w0_val, 2.0))
               + 1.0 / m.obs_ct)
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

    double obs_ct_ext = batchData[idx1[0]].obs_ct;
    if (std::isnan(obs_ct_ext) || obs_ct_ext <= 0.0) continue;
    double var_mu_ext = mu * (1.0 - mu) / obs_ct_ext;

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
        const auto& R = residuals;
        const auto& y = indicator;
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

        double var_mu_ext_loc = mu_local * (1.0 - mu_local) / obs_ct_ext;
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
      Eigen::VectorXd R_tilde_w = residuals.array() - meanR * w_ext;
      double grm_cov_Rext = grm->quadForm(R_tilde_w.data(), static_cast<uint32_t>(nSubj));
      double sumR_sq_over_n = w_ext * w_ext * residuals.sum() * residuals.sum() * 2.0 / obs_ct_ext;
      double num = grm_cov_Rext + sumR_sq_over_n;
      double den = R_tilde_w.array().square().sum() + sumR_sq_over_n;
      var_ratio_ext = (den > 0.0) ? num / den : 1.0;
    }

    // Populate refInfoMap for markers in this group
    for (size_t i : idx1) {
      WtCoxGRefInfo ri;
      ri.AF_ref       = batchData[i].AF_ref;
      ri.obs_ct       = batchData[i].obs_ct;
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
  return "\tp_ext\tp_noext\tz_ext\tz_noext\tp_batch";
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
      info.AF_ref, info.obs_ct, m_cutoff);

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
      info.AF_ref, info.obs_ct, m_cutoff);

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
    double mu_ext, double obs_ct, double p_cut) const {

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

  double var_mu_ext = mu * (1.0 - mu) / obs_ct;
  double var_Sbat = (m_w1.array().square()).sum() * 2.0 * mu * (1.0 - mu) + var_mu_ext;

  double qnorm_val = math::qnorm(1.0 - p_cut / 2.0);
  double lb = -qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);
  double ub =  qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);

  double c_val = math::pnorm(ub / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
  double d_val = math::pnorm(lb / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
  double p_deno = TPR * (std::exp(d_val) * (std::exp(c_val - d_val) - 1.0))
                + (1.0 - TPR) * (1.0 - p_cut);

  // Internal SPA (no sigma2)
  SpaResult spa_s0 = spaGOneSnpHomo(g_input, m_R, mu_ext, obs_ct, b, 0.0, var_ratio0, m_SPA_Cutoff);
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
  SpaResult spa_s1 = spaGOneSnpHomo(g_input, m_R, mu_ext, obs_ct, b, sigma2, var_ratio1, m_SPA_Cutoff);
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
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
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
  auto famIIDs = parseFamIIDs(bfilePrefix + ".fam");
  SubjectData sd(std::move(famIIDs));
  sd.loadResidWtCoxG(residFile);
  sd.finalize();
  infoMsg("  %u subjects loaded", sd.nUsed());

  infoMsg("Loading ref-af file: %s", refAfFile.c_str());
  bool refAfNumeric = false;
  auto refAf = loadRefAfFile(refAfFile, &refAfNumeric);
  infoMsg("  %zu reference records loaded%s", refAf.size(),
          refAfNumeric ? " (numeric fallback)" : "");

  infoMsg("Loading PLINK data: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      sd.usedMask(),
      sd.nFam(),
      sd.nUsed(),
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);
  infoMsg("  %u subjects matched in PLINK, %u markers available",
          plinkData.nSubjUsed(), plinkData.nMarkers());

  infoMsg("Matching markers against reference allele frequencies...");
  auto matched = refAfNumeric
      ? matchMarkersNumeric(plinkData, refAf)
      : matchMarkers(plinkData, refAf);
  infoMsg("  %zu markers matched", matched.size());

  infoMsg("Computing per-marker case/control allele frequencies...");
  computeMarkerStats(matched, plinkData, sd.indicator());

  // ---- Batch-effect testing ----
  infoMsg("Batch-effect testing and parameter estimation...");
  std::unique_ptr<SparseGRM> grm;
  if (!spgrmGctaFile.empty() || !spgrmGrabFile.empty()) {
    infoMsg("  Loading sparse GRM...");
    grm = std::make_unique<SparseGRM>(
        SparseGRM::load(spgrmGrabFile, spgrmGctaFile, sd.usedIIDs(),
                        sd.famIIDs()));
    infoMsg("  Sparse GRM: %u subjects, %zu non-zeros",
            grm->nSubjects(), grm->nnz());
  }
  auto refInfoMap = testBatchEffects(
      matched, sd.residuals(), sd.weights(), sd.indicator(), grm.get(), refPrevalence, cutoff);
  infoMsg("  %zu markers retained after batch-effect QC", refInfoMap->size());

  // ---- Marker-level SPA tests ----
  infoMsg("Running marker-level WtCoxG tests (%d thread(s))...", nthread);
  WtCoxGMethod method(
      sd.residuals(), sd.weights(),
      cutoff, spaCutoff, refInfoMap);
  markerEngine(plinkData, method, outputFile,
               nthread,
               missingCutoff,
               minMafCutoff,
               minMacCutoff,
               /*exactHwe=*/false);
}


// ══════════════════════════════════════════════════════════════════════
// runWtCoxGPheno — new --pheno path: internal regression
// ══════════════════════════════════════════════════════════════════════

#include "wtcoxg/regression.hpp"

void runWtCoxGPheno(
    const std::string& phenoFile,
    const std::string& covarFile,
    const std::vector<std::string>& covarNames,
    const std::string& binaryPheno,
    const std::string& survPheno,
    const std::string& bfilePrefix,
    const std::string& refAfFile,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const std::string& outputFile,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff) {

  const bool isSurv = !survPheno.empty();

  // ---- Load phenotype / covariate data ----
  infoMsg("Loading phenotype file: %s", phenoFile.c_str());
  auto famIIDs = parseFamIIDs(bfilePrefix + ".fam");
  SubjectData sd(std::move(famIIDs));
  sd.loadPhenoFile(phenoFile);
  if (!covarFile.empty()) {
    infoMsg("Loading covariate file: %s", covarFile.c_str());
    sd.loadCovar(covarFile);
  }
  sd.finalize();
  infoMsg("  %u subjects loaded", sd.nUsed());

  const Eigen::Index N = static_cast<Eigen::Index>(sd.nUsed());

  // ---- Extract response ----
  Eigen::VectorXd indicator;   // case/control or event
  Eigen::VectorXd survTime;    // only for Cox

  if (isSurv) {
    std::string timeCol, eventCol;
    regression::parseSurvSpec(survPheno, timeCol, eventCol);
    survTime  = sd.getColumn(timeCol);
    indicator = sd.getColumn(eventCol);
    infoMsg("  Survival phenotype: time=%s, event=%s", timeCol.c_str(), eventCol.c_str());
  } else {
    indicator = sd.getColumn(binaryPheno);
    infoMsg("  Binary phenotype: %s", binaryPheno.c_str());
  }

  // ---- Build design matrices ----
  // Logistic: [1 | covariates] (intercept needed)
  // Cox PH:   [covariates]     (no intercept — absorbed into baseline hazard)
  const int nCov = covarNames.empty() ? 0 : static_cast<int>(covarNames.size());
  Eigen::MatrixXd covarMat;
  if (nCov > 0) {
    covarMat = sd.getColumns(covarNames);
    infoMsg("  %d covariate(s) from --covar-name", nCov);
  }

  // ---- Compute regression weights ----
  Eigen::VectorXd regrWeight = regression::calRegrWeight(refPrevalence, indicator);
  infoMsg("  Regression weights computed (prevalence=%.6f)", refPrevalence);

  // ---- Fit null model and compute residuals ----
  Eigen::VectorXd resid;
  if (isSurv) {
    infoMsg("Fitting weighted Cox PH model...");
    resid = regression::coxResiduals(survTime, indicator, covarMat, regrWeight);
  } else {
    // Logistic needs intercept
    Eigen::MatrixXd designMat(N, 1 + nCov);
    designMat.col(0).setOnes();
    if (nCov > 0) designMat.rightCols(nCov) = covarMat;
    infoMsg("Fitting weighted logistic regression...");
    resid = regression::logisticResiduals(indicator, designMat, regrWeight);
  }
  infoMsg("  Residuals computed (N=%d)", static_cast<int>(N));

  // ---- Set on SubjectData ----
  sd.setResidWeightIndicator(std::move(resid), std::move(regrWeight), indicator);

  // ---- From here, same pipeline as runWtCoxG ----
  infoMsg("Loading ref-af file: %s", refAfFile.c_str());
  bool refAfNumeric = false;
  auto refAf = loadRefAfFile(refAfFile, &refAfNumeric);
  infoMsg("  %zu reference records loaded%s", refAf.size(),
          refAfNumeric ? " (numeric fallback)" : "");

  infoMsg("Loading PLINK data: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      sd.usedMask(),
      sd.nFam(),
      sd.nUsed(),
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);
  infoMsg("  %u subjects matched in PLINK, %u markers available",
          plinkData.nSubjUsed(), plinkData.nMarkers());

  infoMsg("Matching markers against reference allele frequencies...");
  auto matched = refAfNumeric
      ? matchMarkersNumeric(plinkData, refAf)
      : matchMarkers(plinkData, refAf);
  infoMsg("  %zu markers matched", matched.size());

  infoMsg("Computing per-marker case/control allele frequencies...");
  computeMarkerStats(matched, plinkData, sd.indicator());

  // ---- Batch-effect testing ----
  infoMsg("Batch-effect testing and parameter estimation...");
  std::unique_ptr<SparseGRM> grm;
  if (!spgrmGctaFile.empty() || !spgrmGrabFile.empty()) {
    infoMsg("  Loading sparse GRM...");
    grm = std::make_unique<SparseGRM>(
        SparseGRM::load(spgrmGrabFile, spgrmGctaFile, sd.usedIIDs(),
                        sd.famIIDs()));
    infoMsg("  Sparse GRM: %u subjects, %zu non-zeros",
            grm->nSubjects(), grm->nnz());
  }
  auto refInfoMap = testBatchEffects(
      matched, sd.residuals(), sd.weights(), sd.indicator(), grm.get(),
      refPrevalence, cutoff);
  infoMsg("  %zu markers retained after batch-effect QC", refInfoMap->size());

  // ---- Marker-level SPA tests ----
  infoMsg("Running marker-level WtCoxG tests (%d thread(s))...", nthread);
  WtCoxGMethod method(
      sd.residuals(), sd.weights(),
      cutoff, spaCutoff, refInfoMap);
  markerEngine(plinkData, method, outputFile,
               nthread,
               missingCutoff,
               minMafCutoff,
               minMacCutoff,
               /*exactHwe=*/false);
}
