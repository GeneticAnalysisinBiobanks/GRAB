#ifndef SPAMIXPLUS_H
#define SPAMIXPLUS_H

// SPAmixPlus.h -- Extended SPA for admixed populations with sparse GRM and AF model

#include <RcppArmadillo.h>
#include <boost/math/distributions/normal.hpp>
#include <vector>
#include <tuple>
#include <unordered_set>

namespace SPAmixPlus {


class SPAmixPlusClass {
public:

  // ---- Public types ----
  struct AFModelInfo {
    int status;
    arma::vec betas;
  };


  struct PhenoOutlierData {
    arma::uvec posValue;
    arma::uvec posOutlier;
    arma::uvec posNonOutlier;
    arma::vec resid;
    arma::vec residOutlier;
    arma::vec residNonOutlier;
    arma::vec resid2NonOutlier;
  };


  struct RootResult {
    double root;
    int iter;
    bool converge;
    double K2;
  };


  // ---- Construction ----
  SPAmixPlusClass(
    arma::mat residuals,
    arma::mat pcs,
    int sampleSize,
    double spaCutoff,
    std::vector<PhenoOutlierData> outlierList,
    std::vector<std::tuple<int, int, double>> sparseTriplets,
    const std::string& afFilePath = "",
    const std::string& afFilePrecision = "double"
  );


  ~SPAmixPlusClass() = default;


  // ---- Marker-level testing ----
  double getMarkerPval(arma::vec GVec, double altFreq);


  arma::vec getMafEst(arma::vec GVec, double altFreq);


  AFModelInfo computeAFModel(arma::vec GVec, double altFreq);


  arma::vec getAFFromModel(AFModelInfo model, double altFreq);


  double getMarkerPvalFromModel(arma::vec GVec, AFModelInfo model, double altFreq);


  // ---- Regression helpers ----
  arma::vec fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues);


  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);


  arma::vec fit_logistic(const arma::mat& X, const arma::vec& y);


  double calculateSparseVariance(const arma::vec& R_new,
                                  const arma::uvec& posValue);


  // ---- Accessors ----
  int getNpheno() const { return m_Npheno; }
  int getNPCs() const { return m_PCs.n_cols; }


  arma::vec getpvalVec() const { return m_pvalVec; }
  arma::vec getzScoreVec() const { return m_zScoreVec; }


private:


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
  arma::vec m_diffTime1;
  arma::vec m_diffTime2;


};


}

#endif
