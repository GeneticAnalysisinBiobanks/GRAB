#ifndef SPAMIXPLUS_H
#define SPAMIXPLUS_H

// SPAmixPlus.h -- Extended SPA for admixed populations with sparse GRM and AF model

#include <RcppArmadillo.h>
#include <boost/math/distributions/normal.hpp>
#include <vector>
#include <tuple>

#include "mtSPAmix.h"

class mtSPAmixPlusClass {
public:

  struct OutlierData {
    arma::uvec posValue;
    arma::uvec posOutlier;
    arma::uvec posNonOutlier;
    arma::vec resid;
    arma::vec residOutlier;
    arma::vec residNonOutlier;
    arma::vec resid2NonOutlier;
  };

  // ---- Construction ----
  mtSPAmixPlusClass(
    arma::vec residuals,
    arma::mat pcs,
    int sampleSize,
    double spaCutoff,
    OutlierData outlier,
    std::vector<std::tuple<int, int, double>> sparseTriplets,
    std::string afFilePath = "",
    std::string afFilePrecision = "double"
  );

  ~mtSPAmixPlusClass() = default;

  // ---- Marker-level testing ----
  double getMarkerPval(const arma::vec& GVec, double altFreq);

  // ---- Accessors ----

  void getResultVec(const arma::vec& GVec, double altFreq, std::vector<double>& rv) {
    getMarkerPval(GVec, altFreq);
    rv.clear();
    rv.push_back(m_pval);
    rv.push_back(m_zScore);
  }

  int resultSize() const { return 2; }

  std::string getHeaderColumns() const {
    return "\tPvalue\tzScore";
  }

private:

  // ---- Private types ----
  struct AFModelInfo {
    int status;
    arma::vec betas;
  };

  // ---- Private helpers ----
  arma::vec getMafEst(const arma::vec& GVec, double altFreq);
  AFModelInfo computeAFModel(const arma::vec& GVec, double altFreq);
  arma::vec getAFFromModel(const AFModelInfo& model, double altFreq);
  double getMarkerPvalFromModel(const arma::vec& GVec, const AFModelInfo& model, double altFreq);
  arma::vec fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues);
  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);
  double calculateSparseVariance(const arma::vec& R_new);
  int getNPCs() const { return m_PCs.n_cols; }

  // ---- Members ----
  const int m_N;
  const double m_SPA_Cutoff;

  const arma::vec m_resid;
  const arma::mat m_PCs;
  const arma::mat m_onePlusPCs;
  const arma::vec m_sqrt_XTX_inv_diag;

  const OutlierData m_outlier;
  const std::vector<std::tuple<int, int, double>> m_sparseTriplets;
  // Sparse triplets pre-filtered to only (i,j) pairs where both are in posValue.
  // Stored as (posValue-local-index-of-i, posValue-local-index-of-j, grmValue)
  // so calculateSparseVariance needs no set lookup.
  struct FilteredTriplet { arma::uword li; arma::uword lj; double grm; };
  const std::vector<FilteredTriplet> m_filteredTriplets;
  const std::string m_afFilePath;
  const std::string m_afFilePrecision;

  double m_pval;
  double m_zScore;
  double m_Beta;
  double m_S;
  double m_Smean;
  double m_VarS;
  arma::vec m_MAFVec;

  // Pre-allocated scratch buffers (sized at construction, reused across markers)
  arma::vec m_scratch_posValue;       // size = posValue.n_elem
  arma::vec m_scratch_posOutlier;     // size = posOutlier.n_elem
  arma::vec m_scratch_posNonOutlier;  // size = posNonOutlier.n_elem

};

#endif
