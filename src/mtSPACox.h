#ifndef SPACOX_H
#define SPACOX_H

// SPACox.h -- Saddlepoint approximation for Cox proportional hazards model

#include <RcppArmadillo.h>
#include <vector>

class mtSPACoxClass {
private:

  // ---- Interpolation tables (merged from approxfunClass) ----
  // All three cumulant functions share the same x-grid.
  const arma::vec m_xGrid;
  const int m_nGrid;
  const arma::vec m_yK0, m_slopeK0;
  const arma::vec m_yK1, m_slopeK1;
  const arma::vec m_yK2, m_slopeK2;

  // Binary-search piecewise-linear interpolation on shared x-grid.
  double interp(const double* yp, const double* sp, double v) const;
  double interpK0(double v) const { return interp(m_yK0.memptr(), m_slopeK0.memptr(), v); }
  double interpK1(double v) const { return interp(m_yK1.memptr(), m_slopeK1.memptr(), v); }
  double interpK2(double v) const { return interp(m_yK2.memptr(), m_slopeK2.memptr(), v); }

  // ---- Cumulant evaluation (scalar loops, zero heap alloc) ----
  // adjG:  full-length normalized genotype vector (raw pointer)
  // idx:   indices of non-zero subjects; nullptr → scan 0..n-1
  // n:     number of entries to scan
  double evalK0(double t, int N0, double adjG0,
                const double* adjG, const arma::uword* idx, int n) const;
  double evalK1(double t, int N0, double adjG0,
                const double* adjG, const arma::uword* idx, int n, double q2) const;
  double evalK2(double t, int N0, double adjG0,
                const double* adjG, const arma::uword* idx, int n) const;

  // ---- SPA root-finding and probability ----
  struct RootResult { double root; bool converge; double K2; };
  RootResult fastGetRootK1(double initX, int N0, double adjG0,
                           const double* adjG, const arma::uword* idx, int n,
                           double q2) const;
  double getProbSpa(double adjG0, const double* adjG, const arma::uword* idx, int n,
                    int N0, double q2, bool lowerTail) const;

  // ---- Model data ----
  const arma::vec m_mresid;
  const double m_varResid;
  const arma::mat m_XinvXX, m_tX;
  const int m_N;
  const double m_pVal_covaAdj_Cutoff;
  const double m_SPA_Cutoff;

public:
  mtSPACoxClass(
    arma::mat cumul,
    arma::vec mresid,
    arma::mat XinvXX,
    arma::mat tX,
    int N,
    double pVal_covaAdj_Cutoff,
    double SPA_Cutoff
  );

  double getMarkerPval(const arma::vec& GVec, double MAF, double& zScore);
  void getRegionPVec(const arma::vec& GVec, double& zScore, double& pval0, double& pval1,
                     arma::vec& P1Vec, arma::vec& P2Vec);

  // Fills rv with [pval, zScore]
  void getResultVec(const arma::vec& GVec, double altFreq, std::vector<double>& rv) {
    double zScore;
    double pval = getMarkerPval(GVec, altFreq, zScore);
    rv.clear();
    rv.push_back(pval);
    rv.push_back(zScore);
  }

  static int resultSize() { return 2; }

  std::string getHeaderColumns() const {
    return "\tPvalue\tzScore";
  }
};

#endif
