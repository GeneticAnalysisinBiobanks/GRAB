// impute_kernel.cpp — SADGE stage-1 parental-genotype imputation kernels.
#include "sadge/impute_kernel.hpp"

namespace sadge {

void precomputeCoeffs(double mu, PerMarkerCoeffs &c) {
    const double mu2 = mu * mu;
    const double mu3 = mu2 * mu;
    const double mu4 = mu3 * mu;
    c.mu = mu;
    c.A = mu / (2.0 - mu);
    c.B = (3.0 * mu + 1.0) / (1.0 + mu);
    c.C = 1.0 - mu;
    c.D = 2.0 - mu;
    c.E = 1.0 + mu;
    c.F = (1.0 + 4.0 * mu - 2.0 * mu2) / (1.0 + mu - mu2);
    c.G = 1.0 + mu - mu2;
    c.H = mu2 - mu + 3.0;
    c.I = (-2.0 * mu4 + 4.0 * mu3 - 3.0 * mu2 + mu - 14.0) /
          (2.0 * (mu - 2.0) * (mu + 1.0));
    c.J = (4.0 * mu4 - 8.0 * mu3 - mu2 + 5.0 * mu + 6.0) /
          (mu4 - 2.0 * mu3 - 2.0 * mu2 + 3.0 * mu + 2.0);
}

double imputeGpar(
    FamStruct fs,
    double gSelf,
    double gCoSib,
    double gPar1,
    double gPar2,
    const PerMarkerCoeffs &c
) {
    const double mu = c.mu;
    switch (fs) {
    case FamStruct::s2par1sib:
    case FamStruct::s2par2sib:
        // G_par = father + mother (observed dosages).
        return gPar1 + gPar2;

    case FamStruct::s1par1sib: {
        double pp[3], ps[3];
        dsg2prb(gPar1, pp[0], pp[1], pp[2]); // lone parent
        dsg2prb(gSelf, ps[0], ps[1], ps[2]);
        const double E = c.E;
        const double M[3][3] = {
            {mu, E, 4.0 * mu},
            {E, mu + E, 1.0 + E},
            {4.0 * mu, 1.0 + E, 2.0 + E}};
        double s = 0.0;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                s += pp[a] * M[a][b] * ps[b];
        return s;
    }

    case FamStruct::s0par2sib: {
        double q[3], r[3];
        dsg2prb(gSelf, q[0], q[1], q[2]);
        dsg2prb(gCoSib, r[0], r[1], r[2]);
        const double D = c.D, E = c.E, F = c.F;
        const double M[3][3] = {
            {2.0 * mu / D, 2.0 / D, 2.0},
            {2.0 / D, F, 2.0 * (mu + E) / E},
            {2.0, 2.0 * (mu + E) / E, 2.0 * (2.0 * mu + E) / E}};
        double s = 0.0;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                s += q[a] * M[a][b] * r[b];
        return s;
    }

    case FamStruct::s1par2sib: {
        double pp[3], q[3], r[3];
        dsg2prb(gPar1, pp[0], pp[1], pp[2]); // lone parent
        dsg2prb(gSelf, q[0], q[1], q[2]);
        dsg2prb(gCoSib, r[0], r[1], r[2]);
        const double A = c.A, B = c.B, E = c.E;
        const double t = 2.0 * mu;
        // W[parentGeno][s1][s2]; default 2*mu, with the specified overrides
        // (generate_Gpar_1par2sib).
        double W[3][3][3];
        for (int a = 0; a < 3; ++a)
            for (int x = 0; x < 3; ++x)
                for (int y = 0; y < 3; ++y)
                    W[a][x][y] = t;
        // parent = 0
        W[0][0][0] = A;
        W[0][0][1] = 1.0;
        W[0][1][0] = 1.0;
        W[0][1][1] = B;
        // parent = 1
        W[1][0][0] = 1.0 + A;
        W[1][0][1] = E;
        W[1][1][0] = E;
        W[1][0][2] = 2.0;
        W[1][2][0] = 2.0;
        W[1][1][1] = mu + E;
        W[1][1][2] = 1.0 + E;
        W[1][2][1] = 1.0 + E;
        W[1][2][2] = 1.0 + B;
        // parent = 2
        W[2][1][1] = 2.0 + A;
        W[2][1][2] = 3.0;
        W[2][2][1] = 3.0;
        W[2][2][2] = 2.0 + B;
        double s = 0.0;
        for (int a = 0; a < 3; ++a) {
            double inner = 0.0;
            for (int x = 0; x < 3; ++x)
                for (int y = 0; y < 3; ++y)
                    inner += q[x] * W[a][x][y] * r[y];
            s += pp[a] * inner;
        }
        return s;
    }

    case FamStruct::s0par1sib:
    default:
        // Singleton: G_par = G_observed + 2*mu, applied by the caller.
        return 0.0;
    }
}

} // namespace sadge
