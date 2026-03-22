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
    const std::string& afFilePath = "",
    const std::string& afFilePrecision = "double"
  );

  ~mtSPAmixPlusClass() = default;

  // ---- Marker-level testing ----
  double getMarkerPval(arma::vec GVec, double altFreq);

  // ---- Accessors ----

  void getResultVec(arma::vec& GVec, double altFreq, std::vector<double>& rv) {
    getMarkerPval(std::move(GVec), altFreq);
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
  double m_SPA_Cutoff;

  arma::vec m_resid;
  arma::mat m_PCs;
  arma::mat m_onePlusPCs;
  arma::vec m_sqrt_XTX_inv_diag;

  OutlierData m_outlier;
  std::vector<std::tuple<int, int, double>> m_sparseTriplets;
  std::string m_afFilePath;
  std::string m_afFilePrecision;

  double m_pval;
  double m_zScore;
  double m_Beta;
  double m_S;
  double m_Smean;
  double m_VarS;
  arma::vec m_MAFVec;

};

#endif
