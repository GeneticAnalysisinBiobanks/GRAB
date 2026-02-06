/**
 * @file SPAmixPlus.h
 * @brief SPAmixPlus Step 2: Marker-level association testing with pre-computed AF models
 *
 * This file implements Step 2 of SPAmixPlus workflow:
 * - Reads pre-computed AF models from Step 0 (exportAFModelInCPP)
 * - Performs marker-level association testing
 * - Uses sparse GRM for kinship correction
 * - Applies SPA correction for accurate p-values
 *
 * @author GRAB Development Team
 * @date 2026
 */

#ifndef SPAMIXPLUS_H
#define SPAMIXPLUS_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "UTIL.h"

#include <vector>
#include <tuple>
#include <unordered_set>

namespace SPAmixPlus {



/**
 * @class SPAmixPlusClass
 * @brief Main class for SPAmixPlus Step 2 marker testing
 *
 * Workflow:
 * 1. Constructor initializes with null model results and optional AF file
 * 2. getMarkerPval() reads AF model from file and computes p-value
 * 3. Results stored in member vectors for output
 */
class SPAmixPlusClass {
public:
  /**
   * @struct AFModelInfo
   * @brief Stores allele frequency estimation model parameters
   */
  struct AFModelInfo {
    int status;       ///< Model type: 0=mean, 1=linear, 2=logistic
    arma::vec betas;  ///< Regression coefficients [intercept, PC1, PC2, ...]
  };

  // ==========================================================================
  // Constructor and Destructor
  // ==========================================================================
  
  /**
   * @brief Construct SPAmixPlus analysis object
   *
   * @param residuals Residual matrix from null model (N x nPheno)
   * @param pcs Principal component matrix for ancestry (N x nPCs)
   * @param sampleSize Total number of samples
   * @param spaCutoff Z-score threshold for applying SPA correction
   * @param outlierList List containing outlier information per phenotype
   * @param sparseGRM Sparse GRM data frame (ID1_Index, ID2_Index, Value)
   * @param afFilePath Path to pre-computed AF model file (empty string if none)
   * @param afFilePrecision Precision of AF file ("double", "single", or "text")
   */
  SPAmixPlusClass(
    const arma::mat& residuals,
    const arma::mat& pcs,
    int sampleSize,
    double spaCutoff,
    const Rcpp::List& outlierList,
    const Rcpp::DataFrame& sparseGRM,
    const std::string& afFilePath = "",
    const std::string& afFilePrecision = "double"
  );

  /**
   * @brief Destructor
   */
  ~SPAmixPlusClass() = default;

  // ==========================================================================
  // Core Analysis Methods
  // ==========================================================================
  
  /**
   * @brief Compute p-value for a marker
   *
   * @param genotypeVec Genotype dosage vector (0/1/2 or continuous)
   * @param altFreq Population-level alternate allele frequency
   * @return P-value (also stored internally in member vectors)
   */
  double getMarkerPval(arma::vec t_GVec, double t_altFreq);

  // ==========================================================================
  // Helper Methods
  // ==========================================================================
  
  /**
   * @brief Estimate MAF vector for a marker
   *
   * @param t_GVec Genotype vector
   * @param t_altFreq Alternate allele frequency
   * @return Vector of individual-specific AF estimates
   */
  arma::vec getMAFest(arma::vec t_GVec, double t_altFreq);
  
  /**
   * @brief Compute AF model for a marker
   *
   * @param t_GVec Genotype vector
   * @param t_altFreq Alternate allele frequency
   * @return AF model structure
   */
  AFModelInfo computeAFModel(arma::vec t_GVec, double t_altFreq);
  
  /**
   * @brief Get AF vector from model
   *
   * @param t_model AF model
   * @param t_altFreq Alternate allele frequency (fallback)
   * @return Vector of individual-specific AF estimates
   */
  arma::vec getAFFromModel(AFModelInfo t_model, double t_altFreq);
  
  /**
   * @brief Compute p-value using pre-computed AF model
   *
   * @param t_GVec Genotype vector
   * @param t_model AF model
   * @param t_altFreq Alternate allele frequency
   * @return P-value
   */
  double getMarkerPvalFromModel(arma::vec t_GVec, AFModelInfo t_model, double t_altFreq);
  
  /**
   * @brief Fit linear model and get beta coefficients
   *
   * @param g Genotype vector
   * @param pvalues Output: p-values for each PC
   * @return Coefficient vector [intercept, PC1, PC2, ...]
   */
  arma::vec fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues);
  
  /**
   * @brief Fit linear model
   *
   * @param g Genotype vector
   * @param pvalues Output: p-values for each PC
   * @return Fitted values
   */
  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);
  
  /**
   * @brief Fit logistic regression using IRLS
   *
   * @param X Predictor matrix
   * @param y Binary response vector
   * @return Coefficient vector [intercept, PC1, PC2, ...]
   */
  arma::vec fit_logistic(const arma::mat& X, const arma::vec& y);

  /**
   * @brief Calculate variance accounting for sparse kinship matrix
   *
   * @param R_new Weighted residuals
   * @param posValue Indices of non-missing subjects
   * @return Variance estimate
   */
  double calculateSparseVariance(const arma::vec& R_new,
                                  const arma::uvec& posValue);

  // ==========================================================================
  // Accessor Methods
  // ==========================================================================
  
  int getNpheno() const { return m_Npheno; }
  int getNPCs() const { return m_PCs.n_cols; }
  
  // Legacy method names for compatibility with existing R interface
  arma::vec getpvalVec() const { return m_pvalVec; }
  arma::vec getzScoreVec() const { return m_zScoreVec; }
  // arma::vec getSVec() const { return m_SVec; }
  // arma::vec getSmeanVec() const { return m_SmeanVec; }
  // arma::vec getVarSVec() const { return m_VarSVec; }

private:
  // ==========================================================================
  // Member Variables
  // ==========================================================================
  
  // Dimensions
  int m_N;                   ///< Total number of samples
  int m_Npheno;              ///< Number of phenotypes analyzed
  double m_SPA_Cutoff;       ///< |Z-score| threshold for applying SPA

  // Null model components
  arma::mat m_resid;         ///< Residual matrix (N x nPheno)
  arma::mat m_PCs;           ///< Principal components (N x nPCs)
  arma::mat m_onePlusPCs;    ///< [1, PCs] for AF reconstruction (N x (nPCs+1))
  arma::vec m_sqrt_XTX_inv_diag; ///< sqrt(diag((X'X)^-1)) for SE calculation

  // Kinship and outliers
  Rcpp::List m_outlierList;  ///< Per-phenotype outlier subject info
  std::vector<std::tuple<int, int, double>> m_sparseTriplets; ///< Sparse GRM triplets

  // AF file info (for pre-computed models)
  std::string m_afFilePath;      ///< Path to AF model file (empty if not using file)
  std::string m_afFilePrecision; ///< Precision format ("double", "single", or "text")

  // Result storage (per phenotype)
  arma::vec m_pvalVec;       ///< P-values
  arma::vec m_zScoreVec;     ///< Z-scores
  arma::vec m_BetaVec;       ///< Effect size estimates
  arma::vec m_SVec;          ///< Score statistics
  arma::vec m_SmeanVec;      ///< Expected scores under null
  arma::vec m_VarSVec;       ///< Score variances
  
  // Temporary storage
  arma::vec m_MAFVec;        ///< Individual-specific AF estimates (temp)
  arma::vec m_diffTime1;     ///< Timing diagnostic
  arma::vec m_diffTime2;     ///< Timing diagnostic


// ============================================================================
// Helper Functions for Saddle Point Approximation
// ============================================================================

arma::vec momentGeneratingFunction(const arma::vec& t, const arma::vec& maf) {
  return arma::pow((1.0 - maf + maf % arma::exp(t)), 2);
}

arma::vec momentGeneratingFunctionDeriv1(const arma::vec& t, const arma::vec& maf) {
  return 2.0 * (maf % arma::exp(t)) % (1.0 - maf + maf % arma::exp(t));
}

arma::vec momentGeneratingFunctionDeriv2(const arma::vec& t, const arma::vec& maf) {
  arma::vec maf_exp_t = maf % arma::exp(t);
  return 2.0 * arma::pow(maf_exp_t, 2) + 2.0 * maf_exp_t % (1.0 - maf + maf_exp_t);
}

arma::vec cumulantGeneratingFunction(const arma::vec& t, const arma::vec& maf) {
  return arma::log(momentGeneratingFunction(t, maf));
}

arma::vec cumulantGeneratingFunctionDeriv1(const arma::vec& t, const arma::vec& maf) {
  return momentGeneratingFunctionDeriv1(t, maf) / momentGeneratingFunction(t, maf);
}

arma::vec cumulantGeneratingFunctionDeriv2(const arma::vec& t, const arma::vec& maf) {
  arma::vec mgf = momentGeneratingFunction(t, maf);
  arma::vec mgfD1 = momentGeneratingFunctionDeriv1(t, maf);
  arma::vec mgfD2 = momentGeneratingFunctionDeriv2(t, maf);
  return (mgf % mgfD2 - arma::pow(mgfD1, 2)) / arma::pow(mgf, 2);
}

arma::vec computeHorgH2(double t, const arma::vec& residuals, const arma::vec& mafVec) {
  arma::vec result(2);
  arma::vec tR = t * residuals;
  arma::vec expTR = arma::exp(tR);
  arma::vec mafExpTR = mafVec % expTR;
  arma::vec mgf = arma::pow((1.0 - mafVec + mafExpTR), 2);
  arma::vec mgfD1 = 2.0 * mafExpTR % (1.0 - mafVec + mafExpTR);
  arma::vec mgfD2 = 2.0 * arma::pow(mafExpTR, 2) + 2.0 * mafExpTR % (1.0 - mafVec + mafExpTR);
  arma::vec cgf = arma::log(mgf);
  arma::vec cgfD2 = (mgf % mgfD2 - arma::pow(mgfD1, 2)) / arma::pow(mgf, 2);
  
  result(0) = arma::sum(cgf);  // Horg
  result(1) = arma::sum(arma::pow(residuals, 2) % cgfD2);  // H2
  return result;
}

arma::vec computeH1AdjH2(double t, const arma::vec& residuals, double s, const arma::vec& mafVec) {
  arma::vec result(2);
  arma::vec tR = t * residuals;
  arma::vec expTR = arma::exp(tR);
  arma::vec mafExpTR = mafVec % expTR;
  arma::vec mgf = arma::pow((1.0 - mafVec + mafExpTR), 2);
  arma::vec mgfD1 = 2.0 * mafExpTR % (1.0 - mafVec + mafExpTR);
  arma::vec mgfD2 = 2.0 * arma::pow(mafExpTR, 2) + 2.0 * mafExpTR % (1.0 - mafVec + mafExpTR);
  arma::vec cgfD1 = mgfD1 / mgf;
  arma::vec cgfD2 = (mgf % mgfD2 - arma::pow(mgfD1, 2)) / arma::pow(mgf, 2);
  
  result(0) = arma::sum(residuals % cgfD1) - s;  // H1_adj
  result(1) = arma::sum(arma::pow(residuals, 2) % cgfD2);  // H2
  return result;
}

Rcpp::List fastRootK1(double initX, double s, const arma::vec& mafOutlier,
                      double meanNonOutlier, double varNonOutlier,
                      const arma::vec& residOutlier) {
  const double kTolerance = 0.001;
  const int kMaxIter = 100;
  
  double x = initX;
  double oldX, oldDiffX;
  double k1 = 0.0, k2 = 0.0, oldK1;
  double diffX = arma::datum::inf;
  bool converge = true;
  int iter = 0;
  
  for (iter = 0; iter < kMaxIter; ++iter) {
    oldX = x;
    oldDiffX = diffX;
    oldK1 = k1;
    
    arma::vec h1AdjH2 = computeH1AdjH2(x, residOutlier, s, mafOutlier);
    k1 = h1AdjH2(0) + meanNonOutlier + varNonOutlier * x;
    k2 = h1AdjH2(1) + varNonOutlier;
    diffX = -k1 / k2;
    
    if (!std::isfinite(k1)) {
      x = arma::datum::inf;
      k2 = 0.0;
      break;
    }
    
    // Damping if sign changed
    if (arma::sign(k1) != arma::sign(oldK1)) {
      while (std::abs(diffX) > std::abs(oldDiffX) - kTolerance) {
        diffX *= 0.5;
      }
    }
    
    if (std::abs(diffX) < kTolerance) break;
    x = oldX + diffX;
  }
  
  if (iter == kMaxIter) converge = false;
  
  return Rcpp::List::create(
    Rcpp::Named("root") = x,
    Rcpp::Named("iter") = iter,
    Rcpp::Named("converge") = converge,
    Rcpp::Named("K2") = k2
  );
}

double getSPAProbability(const arma::vec& mafOutlier, const arma::vec& residOutlier,
                         double s, bool lowerTail, double meanNonOutlier,
                         double varNonOutlier) {
  Rcpp::List rootList = fastRootK1(0.0, s, mafOutlier, meanNonOutlier, 
                                    varNonOutlier, residOutlier);
  double zeta = rootList["root"];
  
  arma::vec k12 = computeHorgH2(zeta, residOutlier, mafOutlier);
  double k1 = k12(0) + meanNonOutlier * zeta + 0.5 * varNonOutlier * std::pow(zeta, 2);
  double k2 = k12(1) + varNonOutlier;
  
  double temp1 = zeta * s - k1;
  double w = arma::sign(zeta) * std::sqrt(2.0 * temp1);
  double v = zeta * std::sqrt(k2);
  
  double pval = R::pnorm(arma::sign(lowerTail - 0.5) * (w + std::log(v / w) / w), 
                         0.0, 1.0, 1, 0);
  return pval;
}


}; // end class SPAmixPlusClass


}  // namespace SPAmixPlus

#endif  // SPAMIXPLUS_H
