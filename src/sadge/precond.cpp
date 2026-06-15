// precond.cpp — SADGE stage-2 per-phenotype preconditioning.
#include "sadge/precond.hpp"

#include "sadge/impute_kernel.hpp"

#include <cmath>

namespace sadge {

void varGenoPar(double mu, double H, double I, double J,
                std::array<double, kNStruct> &out) {
    const double temp = mu * (1.0 - mu);
    out[static_cast<int>(FamStruct::s0par1sib)] = 2.0 * temp;
    out[static_cast<int>(FamStruct::s1par1sib)] = H * temp;
    out[static_cast<int>(FamStruct::s2par1sib)] = 4.0 * temp;
    out[static_cast<int>(FamStruct::s0par2sib)] = J * temp;
    out[static_cast<int>(FamStruct::s1par2sib)] = I * temp;
    out[static_cast<int>(FamStruct::s2par2sib)] = 4.0 * temp;
}

double computeRho(double mu, const std::array<double, kNStruct> &sumSquare) {
    PerMarkerCoeffs c;
    precomputeCoeffs(mu, c);
    std::array<double, kNStruct> vgp;
    varGenoPar(mu, c.H, c.I, c.J, vgp);
    const double covGGpar = 2.0 * mu * (1.0 - mu);
    double num = 0.0, den = 0.0;
    for (int fs = 0; fs < kNStruct; ++fs) {
        const double cov = covGGpar - 0.5 * vgp[fs];
        num += sumSquare[fs] * cov;
        den += sumSquare[fs] * vgp[fs];
    }
    return num / den;
}

PerPhenoData buildPerPhenoData(
    const Eigen::VectorXd &residUnionRaw,
    const std::vector<FamilyGroup> &families,
    uint32_t nUnion,
    const std::vector<int> &unionToNonSing,
    int nNonSing
) {
    PerPhenoData d;
    d.residUnion = Eigen::VectorXd::Zero(nUnion);
    d.residSing = Eigen::VectorXd::Zero(nUnion);
    d.residNonSingComp = Eigen::VectorXd::Zero(nNonSing);

    for (const auto &fam : families) {
        // Present (non-NaN) union siblings of this family for this phenotype.
        int present[2];
        double r[2];
        int np = 0;
        for (int k = 0; k < fam.n; ++k) {
            const int u = fam.u[k];
            const double v = residUnionRaw[u];
            if (!std::isnan(v)) {
                present[np] = u;
                r[np] = v;
                ++np;
            }
        }
        if (np == 0) continue;

        const FamStruct fs2 = famStruct(fam.genoNPar, np);
        const int fsi = static_cast<int>(fs2);
        const bool singleton = sadge::isSingleton(fs2);

        for (int k = 0; k < np; ++k) {
            const int u = present[k];
            const double v = r[k];
            d.residUnion[u] = v;
            d.sumR += v;
            d.sumSquare[fsi] += v * v;
            if (singleton) {
                d.residSing[u] = v;
                d.sumRSing += v;
            } else {
                // fs2-non-singleton ⟹ fs1-non-singleton, so the compact index
                // is always valid here.
                d.residNonSingComp[unionToNonSing[u]] = v;
            }
        }
        if (np == 2) {
            // sum_prod: product of the two present siblings' residuals (2-sib
            // structures only; a 2-present family is never a singleton).
            d.sumProd[fsi] += r[0] * r[1];
        }
    }

    d.Rsq = 0.0;
    d.Rpr = 0.0;
    for (int fs = 0; fs < kNStruct; ++fs) {
        d.Rsq += d.sumSquare[fs];
        d.Rpr += d.sumProd[fs];
    }

    for (int k = 0; k < 1000; ++k) {
        const double mu = 0.5 * static_cast<double>(k) / 999.0;
        d.rhoGrid[k] = computeRho(mu, d.sumSquare);
    }

    return d;
}

} // namespace sadge
