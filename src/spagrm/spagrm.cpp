// spagrm.cpp — nsSPAGRM free functions, SPAGRMClass, and runSPAGRM

#include "spagrm/spagrm.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spagrm/grm_null.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace nsSPAGRM {

// ──────────────────────────────────────────────────────────────────────
// MGF and its first two derivatives
// ──────────────────────────────────────────────────────────────────────

void mgf(
    double t,
    const Eigen::VectorXd &resid_unrelated_outliers,
    const std::vector<std::array<double, 2> > &TwoSubj_resid,
    const std::vector<std::vector<double> > &TwoSubj_rho,
    const std::vector<UpdatedThreeSubj> &threeSubj,
    double MAF,
    MgfWorkspace &ws
) {
    Eigen::Index w = 0;
    const Eigen::Index nu = resid_unrelated_outliers.size();

    // ── Unrelated outliers ─────────────────────────────────────────────
    if (nu > 0) {
        ws.ul_lambda = (t * resid_unrelated_outliers).array().exp().matrix();
        ws.ul_alpha = (Eigen::VectorXd::Constant(nu, 1.0 - MAF) + MAF * ws.ul_lambda);
        ws.ul_alpha_1 = MAF * resid_unrelated_outliers.cwiseProduct(ws.ul_lambda);
        ws.ul_alpha_2 = resid_unrelated_outliers.cwiseProduct(ws.ul_alpha_1);

        ws.MGF0.segment(0, nu) = ws.ul_alpha.cwiseProduct(ws.ul_alpha);
        ws.MGF1.segment(0, nu) = 2.0 * ws.ul_alpha.cwiseProduct(ws.ul_alpha_1);
        ws.MGF2.segment(0, nu) =
            2.0 * (ws.ul_alpha_1.cwiseProduct(ws.ul_alpha_1) + ws.ul_alpha.cwiseProduct(ws.ul_alpha_2));
        w += nu;
    }

    // ── Two-subject families ───────────────────────────────────────────
    {
        double *g0 = ws.MGF0.data();
        double *g1 = ws.MGF1.data();
        double *g2 = ws.MGF2.data();
        const double maf01 = MAF * (1.0 - MAF);
        const double one_minus_MAF = 1.0 - MAF;

        const int n1 = static_cast<int>(TwoSubj_resid.size());
        for (int i = 0; i < n1; ++i) {
            const double *rho = TwoSubj_rho[i].data();
            const size_t ki = TwoSubj_rho[i].size();
            const double R1 = TwoSubj_resid[i][0];
            const double R2 = TwoSubj_resid[i][1];
            const double etR1 = std::exp(t * R1);
            const double etR2 = std::exp(t * R2);
            const double Rsum = R1 + R2;
            const double R1sq = R1 * R1;
            const double R2sq = R2 * R2;
            const double Rssq = Rsum * Rsum;
            const double etR12 = etR1 * etR2;

            for (size_t j = 0; j < ki; ++j) {
                const double tj = (1.0 - rho[j]) * maf01;
                const double m1j = etR1 * tj;
                const double m2j = etR2 * tj;
                const double m3j = etR12 * (MAF - tj);
                const size_t idx = static_cast<size_t>(w) + j;
                g0[idx] = m1j + m2j + m3j - tj + one_minus_MAF;
                g1[idx] = R1 * m1j + R2 * m2j + Rsum * m3j;
                g2[idx] = R1sq * m1j + R2sq * m2j + Rssq * m3j;
            }
            w += static_cast<Eigen::Index>(ki);
        }
    }

    // ── Three-subject families ─────────────────────────────────────────
    const int n2 = static_cast<int>(threeSubj.size());
    for (int i = 0; i < n2; ++i) {
        const double *ss = threeSubj[i].stand_S.data();
        const double *ap = threeSubj[i].arr_prob.data();
        const size_t ns = threeSubj[i].stand_S.size();
        double s0 = 0.0, s1 = 0.0, s2 = 0.0;
        for (size_t j = 0; j < ns; ++j) {
            const double m0j = std::exp(t * ss[j]) * ap[j];
            const double m1j = ss[j] * m0j;
            s0 += m0j;
            s1 += m1j;
            s2 += ss[j] * m1j;
        }
        ws.MGF0[w] = s0;
        ws.MGF1[w] = s1;
        ws.MGF2[w] = s2;
        ++w;
    }
}

// ──────────────────────────────────────────────────────────────────────
// Newton-Raphson root finder
// ──────────────────────────────────────────────────────────────────────

static inline double signOf(double x) {
    return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
}

double fastGetRoot(
    const Eigen::VectorXd &resid_unrelated_outliers,
    const std::vector<std::array<double, 2> > &TwoSubj_resid,
    const std::vector<std::vector<double> > &TwoSubj_rho,
    const std::vector<UpdatedThreeSubj> &threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score,
    double MAF,
    double init_t,
    double tol,
    MgfWorkspace &ws,
    int maxiter
) {
    double t = init_t;
    double CGF1 = 0.0, CGF2 = 0.0;
    double diff_t = std::numeric_limits<double>::infinity();

    const double mean = 2.0 * MAF * sum_R_nonOutlier;
    const double var = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

    for (int iter = 1; iter < maxiter; ++iter) {
        const double old_t = t;
        const double old_diff_t = diff_t;
        const double old_CGF1 = CGF1;

        mgf(t, resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj, MAF, ws);

        // temp = MGF1 ./ MGF0
        const Eigen::Index sz = ws.MGF0.size();
        for (Eigen::Index k = 0; k < sz; ++k)
            ws.temp[k] = ws.MGF1[k] / ws.MGF0[k];

        CGF1 = ws.temp.sum() + mean + var * t - Score;
        // CGF2 = sum(MGF2/MGF0) - sum(temp^2) + var
        double s1 = 0.0, s2 = 0.0;
        for (Eigen::Index k = 0; k < sz; ++k) {
            s1 += ws.MGF2[k] / ws.MGF0[k];
            s2 += ws.temp[k] * ws.temp[k];
        }
        CGF2 = s1 - s2 + var;

        diff_t = -CGF1 / CGF2;

        if (std::isnan(diff_t) || std::isinf(CGF2)) {
            t = t / 2.0;
            diff_t = std::min(std::abs(t), 1.0) * signOf(Score);
            continue;
        }

        if (std::isnan(old_CGF1) || (signOf(old_CGF1) != 0 && signOf(CGF1) != signOf(old_CGF1))) {
            if (std::abs(diff_t) < tol) {
                t = old_t + diff_t;
                break;
            } else {
                for (int bisect = 0; bisect < 200 && std::abs(old_diff_t) > tol &&
                     std::abs(diff_t) > std::abs(old_diff_t) - tol; ++bisect)
                    diff_t /= 2.0;
                t = old_t + diff_t;
                continue;
            }
        }

        if (signOf(Score) != signOf(old_t + diff_t) && (signOf(old_CGF1) == 0 || signOf(CGF1) == signOf(old_CGF1))) {
            for (int bisect = 0; bisect < 200 && signOf(Score) != signOf(old_t + diff_t); ++bisect)
                diff_t /= 2.0;
            t = old_t + diff_t;
            continue;
        }

        t = old_t + diff_t;
        if (std::abs(diff_t) < tol) break;
    }

    // Final mgf evaluation at convergence.
    mgf(t, resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj, MAF, ws);
    return t;
}

// ──────────────────────────────────────────────────────────────────────
// Saddlepoint approximation tail probability
// ──────────────────────────────────────────────────────────────────────

double getProbSpa(
    const Eigen::VectorXd &resid_unrelated_outliers,
    const std::vector<std::array<double, 2> > &TwoSubj_resid,
    const std::vector<std::vector<double> > &TwoSubj_rho,
    const std::vector<UpdatedThreeSubj> &threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score,
    double MAF,
    bool lower_tail,
    double zeta,
    double tol,
    MgfWorkspace &ws
) {
    zeta = fastGetRoot(resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj, sum_R_nonOutlier,
                       R_GRM_R_nonOutlier, Score, MAF, zeta, tol, ws);

    const double mean = 2.0 * MAF * sum_R_nonOutlier;
    const double var = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

    const Eigen::Index sz = ws.MGF0.size();
    for (Eigen::Index k = 0; k < sz; ++k)
        ws.temp[k] = ws.MGF1[k] / ws.MGF0[k];

    // CGF0 = sum(log(MGF0)) + mean*zeta + 0.5*var*zeta^2
    double logSum = 0.0;
    for (Eigen::Index k = 0; k < sz; ++k)
        logSum += std::log(ws.MGF0[k]);
    const double CGF0 = logSum + mean * zeta + 0.5 * var * zeta * zeta;

    // CGF2 = sum(MGF2/MGF0) - sum(temp^2) + var
    double s1 = 0.0, s2 = 0.0;
    for (Eigen::Index k = 0; k < sz; ++k) {
        s1 += ws.MGF2[k] / ws.MGF0[k];
        s2 += ws.temp[k] * ws.temp[k];
    }
    const double CGF2_val = s1 - s2 + var;

    const double w_val = signOf(zeta) * std::sqrt(2.0 * (zeta * Score - CGF0));
    const double v_val = zeta * std::sqrt(CGF2_val);
    const double u = w_val + (1.0 / w_val) * std::log(v_val / w_val);

    return math::pnorm(u, 0.0, 1.0, lower_tail);
}

} // namespace nsSPAGRM

// ══════════════════════════════════════════════════════════════════════
// SPAGRMClass
// ══════════════════════════════════════════════════════════════════════

SPAGRMClass::SPAGRMClass(
    Eigen::VectorXd resid,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double R_GRM_R_TwoSubjOutlier,
    double R_GRM_R,
    std::vector<double> MAF_interval,
    nsSPAGRM::FamilyData fam,
    double SPA_Cutoff,
    double zeta,
    double tol
)
{
    auto sd = std::make_shared<SharedData>();
    sd->resid = std::move(resid);
    sd->resid_unrelated_outliers = std::move(fam.resid_unrelated_outliers);
    sd->sum_unrelated_outliers2 = sd->resid_unrelated_outliers.squaredNorm();
    sd->sum_R_nonOutlier = sum_R_nonOutlier;
    sd->R_GRM_R_nonOutlier = R_GRM_R_nonOutlier;
    sd->R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier;
    sd->R_GRM_R = R_GRM_R;
    sd->resid_sum = sd->resid.sum();
    sd->MAF_interval = std::move(MAF_interval);
    sd->TwoSubj_resid_list = std::move(fam.twoSubj_resid);
    sd->TwoSubj_rho_list = std::move(fam.twoSubj_rho);
    sd->ThreeSubj_standS_list = std::move(fam.threeSubj_standS);
    sd->ThreeSubj_CLT_list = std::move(fam.threeSubj_CLT);
    sd->SPA_Cutoff = SPA_Cutoff;
    sd->zeta = zeta;
    sd->tol = tol;
    m_shared = std::move(sd);
    rebuildScratch();
}

SPAGRMClass::SPAGRMClass(const SPAGRMClass &o)
    : m_shared(o.m_shared)
{
    rebuildScratch();
}

void SPAGRMClass::rebuildScratch() {
    const auto &s = *m_shared;
    const auto n_unrel = s.resid_unrelated_outliers.size();
    const auto mgfSz = static_cast<Eigen::Index>(
        nsSPAGRM::mgfOutputSize(static_cast<size_t>(n_unrel), s.TwoSubj_rho_list, s.ThreeSubj_standS_list.size()));
    m_workspace = nsSPAGRM::MgfWorkspace(mgfSz, n_unrel);

    const int n3 = static_cast<int>(s.ThreeSubj_standS_list.size());
    m_threeSubj_scratch.resize(n3);
    for (int i = 0; i < n3; ++i) {
        m_threeSubj_scratch[i].stand_S = s.ThreeSubj_standS_list[i];
        m_threeSubj_scratch[i].arr_prob.resize(s.ThreeSubj_standS_list[i].size());
    }
}

double SPAGRMClass::getMarkerPval(
    const Eigen::VectorXd &GVec,
    double altFreq,
    double &zScore,
    double gMean
) {
    const auto &s = *m_shared;
    const double MAF = std::min(altFreq, 1.0 - altFreq);
    if (std::isnan(gMean))
        gMean = GVec.mean();
    const double Score = GVec.dot(s.resid) - gMean * s.resid_sum;
    const double G_var = 2.0 * MAF * (1.0 - MAF);
    const double Score_var = G_var * s.R_GRM_R;

    // Guard: monomorphic or degenerate marker — skip SPA entirely
    if (Score_var <= 0.0 || MAF <= 0.0) {
        zScore = 0.0;
        return std::numeric_limits<double>::quiet_NaN();
    }

    zScore = Score / std::sqrt(Score_var);

    if (!std::isfinite(zScore)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (std::abs(zScore) <= s.SPA_Cutoff) {
        return 2.0 * math::pnorm(std::abs(zScore), 0.0, 1.0, /*lower_tail=*/ false);
    }

    int order2 =
        static_cast<int>(std::lower_bound(s.MAF_interval.begin(), s.MAF_interval.end(), MAF) - s.MAF_interval.begin());
    // Clamp to valid interpolation range
    const int nBins = static_cast<int>(s.MAF_interval.size());
    if (order2 <= 0) order2 = 1;
    if (order2 >= nBins) order2 = nBins - 1;
    const int order1 = order2 - 1;
    const double MAF_ratio = (s.MAF_interval[order2] - MAF) / (s.MAF_interval[order2] - s.MAF_interval[order1]);
    const double one_minus_mr = 1.0 - MAF_ratio;

    // Update arr_prob in-place for three-subject families.
    double Var_ThreeOutlier = 0.0;
    const int n3 = static_cast<int>(m_threeSubj_scratch.size());
    for (int i = 0; i < n3; ++i) {
        const double *c1 = s.ThreeSubj_CLT_list[i].col(order1).data();
        const double *c2 = s.ThreeSubj_CLT_list[i].col(order2).data();
        const double *ss = m_threeSubj_scratch[i].stand_S.data();
        double *ap = m_threeSubj_scratch[i].arr_prob.data();
        const size_t sz = m_threeSubj_scratch[i].stand_S.size();
        double s1 = 0.0, s2 = 0.0;
        for (size_t k = 0; k < sz; ++k) {
            ap[k] = MAF_ratio * c1[k] + one_minus_mr * c2[k];
            const double tmp = ss[k] * ap[k];
            s1 += tmp;
            s2 += ss[k] * tmp;
        }
        Var_ThreeOutlier += s2 - s1 * s1;
    }

    const double EmpVar =
        G_var * (s.R_GRM_R_nonOutlier + s.sum_unrelated_outliers2 + s.R_GRM_R_TwoSubjOutlier) + Var_ThreeOutlier;
    const double Score_adj = Score * std::sqrt(EmpVar / Score_var);

    double zeta1 = std::abs(Score_adj) / Score_var;
    zeta1 = std::min(zeta1, 1.2);
    const double zeta2 = -std::abs(s.zeta);

    const double pval1 = nsSPAGRM::getProbSpa(s.resid_unrelated_outliers, s.TwoSubj_resid_list, s.TwoSubj_rho_list,
                                              m_threeSubj_scratch, s.sum_R_nonOutlier, s.R_GRM_R_nonOutlier,
                                              std::abs(Score_adj), MAF, false, zeta1, 1e-4, m_workspace);
    const double pval2 = nsSPAGRM::getProbSpa(s.resid_unrelated_outliers, s.TwoSubj_resid_list, s.TwoSubj_rho_list,
                                              m_threeSubj_scratch, s.sum_R_nonOutlier, s.R_GRM_R_nonOutlier,
                                              -std::abs(Score_adj), MAF, true, zeta2, s.tol, m_workspace);
    return pval1 + pval2;
}

double SPAGRMClass::getMarkerPvalFromScore(
    double Score,
    double altFreq,
    double &zScore
) {
    const auto &s = *m_shared;
    const double MAF = std::min(altFreq, 1.0 - altFreq);
    const double G_var = 2.0 * MAF * (1.0 - MAF);
    const double Score_var = G_var * s.R_GRM_R;

    if (Score_var <= 0.0 || MAF <= 0.0) {
        zScore = 0.0;
        return std::numeric_limits<double>::quiet_NaN();
    }

    zScore = Score / std::sqrt(Score_var);

    if (!std::isfinite(zScore)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (std::abs(zScore) <= s.SPA_Cutoff) {
        return 2.0 * math::pnorm(std::abs(zScore), 0.0, 1.0, /*lower_tail=*/ false);
    }

    int order2 =
        static_cast<int>(std::lower_bound(s.MAF_interval.begin(), s.MAF_interval.end(), MAF) - s.MAF_interval.begin());
    const int nBins = static_cast<int>(s.MAF_interval.size());
    if (order2 <= 0) order2 = 1;
    if (order2 >= nBins) order2 = nBins - 1;
    const int order1 = order2 - 1;
    const double MAF_ratio = (s.MAF_interval[order2] - MAF) / (s.MAF_interval[order2] - s.MAF_interval[order1]);
    const double one_minus_mr = 1.0 - MAF_ratio;

    double Var_ThreeOutlier = 0.0;
    const int n3 = static_cast<int>(m_threeSubj_scratch.size());
    for (int i = 0; i < n3; ++i) {
        const double *c1 = s.ThreeSubj_CLT_list[i].col(order1).data();
        const double *c2 = s.ThreeSubj_CLT_list[i].col(order2).data();
        const double *ss = m_threeSubj_scratch[i].stand_S.data();
        double *ap = m_threeSubj_scratch[i].arr_prob.data();
        const size_t sz = m_threeSubj_scratch[i].stand_S.size();
        double s1 = 0.0, s2 = 0.0;
        for (size_t k = 0; k < sz; ++k) {
            ap[k] = MAF_ratio * c1[k] + one_minus_mr * c2[k];
            const double tmp = ss[k] * ap[k];
            s1 += tmp;
            s2 += ss[k] * tmp;
        }
        Var_ThreeOutlier += s2 - s1 * s1;
    }

    const double EmpVar =
        G_var * (s.R_GRM_R_nonOutlier + s.sum_unrelated_outliers2 + s.R_GRM_R_TwoSubjOutlier) + Var_ThreeOutlier;
    const double Score_adj = Score * std::sqrt(EmpVar / Score_var);

    double zeta1 = std::abs(Score_adj) / Score_var;
    zeta1 = std::min(zeta1, 1.2);
    const double zeta2 = -std::abs(s.zeta);

    const double pval1 = nsSPAGRM::getProbSpa(s.resid_unrelated_outliers, s.TwoSubj_resid_list, s.TwoSubj_rho_list,
                                              m_threeSubj_scratch, s.sum_R_nonOutlier, s.R_GRM_R_nonOutlier,
                                              std::abs(Score_adj), MAF, false, zeta1, 1e-4, m_workspace);
    const double pval2 = nsSPAGRM::getProbSpa(s.resid_unrelated_outliers, s.TwoSubj_resid_list, s.TwoSubj_rho_list,
                                              m_threeSubj_scratch, s.sum_R_nonOutlier, s.R_GRM_R_nonOutlier,
                                              -std::abs(Score_adj), MAF, true, zeta2, s.tol, m_workspace);
    return pval1 + pval2;
}

// ══════════════════════════════════════════════════════════════════════
// Helper — build edges, components, family entries from GRM entries
// ══════════════════════════════════════════════════════════════════════
namespace {

struct GRMTopology {
    std::unordered_set<uint32_t> singletonSet;
    std::vector<std::vector<uint32_t> > families;
    std::vector<std::vector<SparseGRM::Entry> > familyEntries;
};

GRMTopology buildTopology(
    uint32_t N,
    const std::vector<SparseGRM::Entry> &entries
) {
    std::vector<std::pair<uint32_t, uint32_t> > edges;
    {
        std::unordered_set<uint64_t> seen;
        for (const auto &e : entries) {
            if (e.row == e.col) continue;
            uint32_t lo = std::min(e.row, e.col);
            uint32_t hi = std::max(e.row, e.col);
            uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
            if (seen.insert(key).second) edges.push_back({lo, hi});
        }
    }
    auto components = nsGRMNull::getComponents(N, edges);

    GRMTopology topo;
    std::vector<std::vector<uint32_t> > singletons;
    for (auto &comp : components) {
        if (comp.size() == 1)
            singletons.push_back(std::move(comp));
        else
            topo.families.push_back(std::move(comp));
    }
    for (const auto &s : singletons)
        topo.singletonSet.insert(s[0]);

    std::unordered_map<uint32_t, size_t> subjToFamily;
    for (size_t fi = 0; fi < topo.families.size(); ++fi)
        for (uint32_t idx : topo.families[fi])
            subjToFamily[idx] = fi;

    topo.familyEntries.resize(topo.families.size());
    for (const auto &e : entries) {
        auto it = subjToFamily.find(e.row);
        if (it != subjToFamily.end()) topo.familyEntries[it->second].push_back(e);
    }
    return topo;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// runSPAGRM — entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAGRM(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile
) {
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    sd.loadResidOne(phenoFile, residNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();
    const uint32_t N = sd.nUsed();
    infoMsg("  %u subjects in union mask", N);

    auto subjIDs = sd.usedIIDs();
    auto subjIdMap = text::buildIIDMap(subjIDs);

    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjIDs, sd.famIIDs());
    infoMsg("Sparse GRM: %zu entries (diagonal + off-diag)", grm.nnz());

    infoMsg("Loading pairwise IBD: %s", pairwiseIBDFile.c_str());
    auto ibdEntries = nsGRMNull::loadIndexedIBD(pairwiseIBDFile, subjIdMap);
    infoMsg("Loaded %zu IBD records", ibdEntries.size());

    const auto &allEntries = grm.entries();
    const auto &grmDiag = grm.diagonal();

    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    infoMsg("Genotype data: %u markers, %u subjects", genoData->nMarkers(), genoData->nSubjUsed());

    // ---- Build per-phenotype tasks ----
    auto phenoInfos = sd.buildPerColumnMasks();
    const int K = sd.residOneCols();
    if (K > 1) infoMsg("Multi-column residual file: %d phenotypes", K);

    std::vector<PhenoTask> tasks(K);

    if (K == 1) {
        // Single phenotype: use union structures directly
        auto topo = buildTopology(N, allEntries);
        auto ibdPairMap = nsGRMNull::buildIBDPairMap(ibdEntries);
        infoMsg("Building null model for '%s'...", phenoInfos[0].name.c_str());
        SPAGRMClass sg = nsGRMNull::buildSPAGRMNullModel(sd.residuals(), N, topo.singletonSet, grmDiag, topo.families,
                                                         topo.familyEntries, allEntries, ibdEntries, ibdPairMap,
                                                         spaCutoff, minMafCutoff, minMacCutoff, nthreads);
        tasks[0].phenoName = phenoInfos[0].name;
        tasks[0].method = std::make_unique<SPAGRMMethod>(std::move(sg));
        tasks[0].unionToLocal = phenoInfos[0].unionToLocal;
        tasks[0].nUsed = phenoInfos[0].nUsed;
        infoMsg("  Phenotype '%s': %u subjects", phenoInfos[0].name.c_str(), phenoInfos[0].nUsed);
    } else {
        for (int rc = 0; rc < K; ++rc) {
            const auto &pi = phenoInfos[rc];
            const auto &u2l = pi.unionToLocal;
            const uint32_t Np = pi.nUsed;

            // Re-index GRM entries to per-phenotype dense indices
            std::vector<SparseGRM::Entry> pEntries;
            std::vector<double> pDiag(Np, 1.0);
            for (const auto &e : allEntries) {
                uint32_t li = u2l[e.row], lj = u2l[e.col];
                if (li != UINT32_MAX && lj != UINT32_MAX) {
                    pEntries.push_back({li, lj, e.value});
                    if (li == lj) pDiag[li] = e.value;
                }
            }

            auto topo = buildTopology(Np, pEntries);

            // Re-index IBD entries
            std::vector<nsGRMNull::IndexedIBD> pIbd;
            for (const auto &ibd : ibdEntries) {
                uint32_t li = u2l[ibd.idx1], lj = u2l[ibd.idx2];
                if (li != UINT32_MAX && lj != UINT32_MAX) pIbd.push_back({li, lj, ibd.pa, ibd.pb, ibd.pc});
            }
            auto pIbdMap = nsGRMNull::buildIBDPairMap(pIbd);

            Eigen::VectorXd phenoResid = extractPhenoVec(sd.residMatrix().col(rc), pi);

            infoMsg("Building null model for '%s' (%u subjects)...", pi.name.c_str(), Np);
            SPAGRMClass sg = nsGRMNull::buildSPAGRMNullModel(phenoResid, Np, topo.singletonSet, pDiag, topo.families,
                                                             topo.familyEntries, pEntries, pIbd, pIbdMap, spaCutoff,
                                                             minMafCutoff, minMacCutoff, nthreads);

            tasks[rc].phenoName = pi.name;
            tasks[rc].method = std::make_unique<SPAGRMMethod>(std::move(sg));
            tasks[rc].unionToLocal = pi.unionToLocal;
            tasks[rc].nUsed = pi.nUsed;
            infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
        }
    }

    infoMsg("Running SPAGRM marker tests (%d thread(s), %d phenotype(s))...", nthreads, K);
    multiPhenoEngine(*genoData, tasks, outPrefix, "SPAGRM", compression, compressionLevel, nthreads,
                     missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
}
