#ifndef SPAMIXPLUS_H
#define SPAMIXPLUS_H

// SPAmixPlus.h -- Extended SPA for admixed populations with sparse GRM and AF model

#include <RcppArmadillo.h>
#include <boost/math/distributions/normal.hpp>
#include <vector>
#include <tuple>
#include <unordered_set>

#include "mtSPAmix.h"

class mtSPAmixPlusClass {
public:

  struct PhenoOutlierData {
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
    arma::mat residuals,
    arma::mat pcs,
    int sampleSize,
    double spaCutoff,
    std::vector<PhenoOutlierData> outlierList,
    std::vector<std::tuple<int, int, double>> sparseTriplets,
    const std::string& afFilePath = "",
    const std::string& afFilePrecision = "double"
  );

  ~mtSPAmixPlusClass() = default;

  // ---- Marker-level testing ----
  double getMarkerPval(arma::vec GVec, double altFreq);

  // ---- Accessors ----
  int getNpheno() const { return m_Npheno; }
  arma::vec getpvalVec() const { return m_pvalVec; }
  arma::vec getzScoreVec() const { return m_zScoreVec; }

private:

  // ---- Private types ----
  struct AFModelInfo {
    int status;
    arma::vec betas;
  };

  // ---- Private helpers ----
  arma::vec getMafEst(arma::vec GVec, double altFreq);
  AFModelInfo computeAFModel(arma::vec GVec, double altFreq);
  arma::vec getAFFromModel(AFModelInfo model, double altFreq);
  double getMarkerPvalFromModel(arma::vec GVec, AFModelInfo model, double altFreq);
  arma::vec fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues);
  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);
  double calculateSparseVariance(const arma::vec& R_new, const arma::uvec& posValue);
  int getNPCs() const { return m_PCs.n_cols; }

  // ---- Members ----
  int m_N;
  int m_Npheno;
  double m_SPA_Cutoff;

  arma::mat m_resid;
  arma::mat m_PCs;
  arma::mat m_onePlusPCs;
  arma::vec m_sqrt_XTX_inv_diag;

  std::vector<PhenoOutlierData> m_outlierList;
  std::vector<std::tuple<int, int, double>> m_sparseTriplets;
  std::string m_afFilePath;
  std::string m_afFilePrecision;

  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;
  arma::vec m_BetaVec;
  arma::vec m_SVec;
  arma::vec m_SmeanVec;
  arma::vec m_VarSVec;
  arma::vec m_MAFVec;

};

#endif
