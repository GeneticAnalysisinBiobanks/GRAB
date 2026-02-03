#ifndef SPAmixPlus_HPP
#define SPAmixPlus_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "UTIL.h"

#include <vector>
#include <tuple>
#include <unordered_set>

namespace SPAmixPlus {

class SPAmixPlusClass {
public:
    struct AFModelInfo {
      int status;     // 0: Mean, 1: Linear, 2: Logistic, 3: Mean (fallback), 4: Linear (fallback) - simplified to 0,1,2
      arma::vec betas;
    };

    // Constructor matching Main.cpp
    SPAmixPlusClass(arma::mat t_resid,
                    arma::mat t_PCs,
                    int t_N,
                    double t_SPA_Cutoff,
                    Rcpp::List t_outlierList,
                    Rcpp::DataFrame t_sparseGRM,
                    Rcpp::DataFrame t_ResidMat);

    // Core method
    double getMarkerPval(arma::vec t_GVec, double t_altFreq);
    double getMarkerPvalFromModel(arma::vec t_GVec, AFModelInfo t_model); // Optimized version
    
    // Internal helpers
    arma::vec getMAFest(arma::vec t_GVec, double t_altFreq);
    AFModelInfo computeAFModel(arma::vec t_GVec, double t_altFreq); // Logic extraction
    arma::vec getAFFromModel(AFModelInfo t_model, double t_altFreq); // Reconstruction

    double getMarkerPvalFromModel(arma::vec t_GVec, AFModelInfo t_model, double t_altFreq);

    double calculateSparseVariance(const arma::vec& R_new, const arma::uvec& posValue);
    arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);
    arma::vec fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues); // Helper returning beta


    // Getters
    arma::vec getpvalVec() { return m_pvalVec; };
    arma::vec getzScoreVec() { return m_zScoreVec; };
    arma::vec getBetaVec() { return m_BetaVec; };
    arma::vec getSVec() { return m_SVec; };
    arma::vec getSmeanVec() { return m_SmeanVec; };
    arma::vec getVarSVec() { return m_VarSVec; };
    int getNpheno() { return m_Npheno; };
    int getNPCs() { return m_PCs.n_cols; };

private:
    int m_N;
    int m_Npheno;
    double m_SPA_Cutoff;

    arma::mat m_resid;
    arma::mat m_PCs;
    arma::mat m_onePlusPCs;
    arma::vec m_sqrt_XTX_inv_diag;

    Rcpp::List m_outlierList;
    arma::mat m_ResidMat;
    arma::ivec m_subjIndices;
    std::vector<std::tuple<int, int, double>> m_sparseTriplets;

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
