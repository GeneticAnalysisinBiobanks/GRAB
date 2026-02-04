#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <RcppArmadillo.h>
#include "Main.h"
#include "UTIL.h"

// ============================================================================
// External C++ Functions for R Interface
// ============================================================================

static arma::mat g_pcs;
static int g_npcs;
static int g_N;

static arma::mat g_X;           // Design matrix [1, g_pcs]
static arma::vec g_sqrtXtXInvDiag;  // sqrt(diag((X'X)^-1))

static const double kTolerance = 1e-6;
static const int kMaxIter = 100;


static arma::vec fitLogisticRegression(const arma::mat& X, const arma::vec& y) {

  // 1. Start with initial beta estimates
  arma::vec beta(g_X.n_cols, arma::fill::zeros);
  arma::vec mu(g_N);
  arma::mat weightedX(g_N, g_X.n_cols);
  
  for (int iter = 0; iter < kMaxIter; ++iter) {

    // 2. Compute predicted probabilities (mu)
    mu = 1.0 / (1.0 + arma::exp(-g_X * beta));

    // 3. Compute weights = mu * (1-mu)  
    arma::vec weights = mu % (1.0 - mu);
    arma::vec workingResp = g_X * beta + (y - mu) / weights;
    
    for (int j = 0; j < g_X.n_cols; ++j) {
      weightedX.col(j) = g_X.col(j) % weights;
    }
    // 4. Solve weighted least squares: beta_new = solve(X'W X, X'W z)
    arma::vec betaNew = arma::solve(
      g_X.t() * weightedX, 
      g_X.t() * (weights % workingResp)
    );
    
    if (arma::norm(betaNew - beta) < kTolerance) {
      // 5. Repeat until convergence
      return betaNew;
    }
    beta = betaNew;
  }
  
  return beta;
}


static void computeAFoneMarker(
    const arma::vec& genotypeVec,
    double& altFreq,
    int& status,
    arma::vec& betas
  ) {
    
  const double kMacCutoff = 20.0;
  const double kPCPvalueCutoff = 0.05;
  const double kNegativeRatioCutoff = 0.1;
  
  // Initialize outputs
  status = 0;  // Default: mean-based
  betas.zeros(g_npcs + 1);
  
  double mac = altFreq * 2.0 * g_N;
  
  // Low MAC: use mean-based estimation
  if (mac <= kMacCutoff) {
    return;
  }
  
  // Try linear model using precomputed design matrix
  arma::vec coefLinear = arma::solve(g_X, genotypeVec);
  arma::vec fittedLinear = g_X * coefLinear / 2.0;  // Scale to [0,1]
  
  // Check boundary violations
  int nErrors = arma::sum(fittedLinear < 0.0) + arma::sum(fittedLinear > 1.0);
  double errorProp = static_cast<double>(nErrors) / g_N;
  
  if (errorProp < kNegativeRatioCutoff) {
    // Linear model is good
    status = 1;
    betas = coefLinear;
    return;
  }
  
  // Calculate p-values for PCs to check significance
  arma::vec fittedValues = g_X * coefLinear;
  double residualSS = arma::sum(arma::pow(genotypeVec - fittedValues, 2));
  double sigma2 = residualSS / (g_N - g_npcs - 1);
  arma::vec se = g_sqrtXtXInvDiag * std::sqrt(sigma2);
  arma::vec tStats = coefLinear / se;
  
  arma::vec pvalues(g_npcs);
  for (int i = 0; i < g_npcs; ++i) {
    pvalues(i) = 2.0 * R::pt(std::abs(tStats(i + 1)), g_N - g_npcs - 1, 0, 0);
  }
  
  // Linear model has issues, check for significant PCs
  arma::uvec significantPCs = arma::find(pvalues < kPCPvalueCutoff);
  
  if (significantPCs.n_elem == 0) {
    // No significant PCs: fall back to mean
    return;
  }
  
  // Try logistic regression with significant PCs
  arma::mat pcsSig = g_pcs.cols(significantPCs);
  arma::vec genotypesBinary(g_N, arma::fill::zeros);
  genotypesBinary.elem(arma::find(genotypeVec > 0.5)).ones();
  
  double macAfter = arma::sum(genotypesBinary);
  if (macAfter <= kMacCutoff) {
    // Still low MAC after dichotomization
    return;
  }
  
  // Fit logistic regression
  status = 2;
  arma::vec betaLogistic = fitLogisticRegression(pcsSig, genotypesBinary);
  
  // Map back to full coefficient vector
  betas(0) = betaLogistic(0);  // Intercept
  for (size_t i = 0; i < significantPCs.n_elem; ++i) {
    betas(significantPCs(i) + 1) = betaLogistic(i + 1);
  }
}


/**
 * @brief Export allele frequency models to file (Step 0)
 *
 * Pre-computes AF estimation models for all markers and saves to disk.
 * These models can be reused in Step 2 for efficient analysis.
 * This is a standalone function that does not require any global state.
 *
 * @param genoType Genotype file format ("PLINK" or "BGEN")
 * @param genoIndex Vector of genotype indices to process
 * @param outputFile Path to output file
 * @param pcs Principal component matrix (N x g_npcs)
 * @param outputFormat Output format:
 *   - "binary64": int32 status + float64 betas (8 bytes/coef, ~15-16 digits precision)
 *   - "binary32": int8 status + float32 betas (4 bytes/coef, ~6-7 digits precision, 50% smaller)
 *   - "text6g": TSV format with 6 significant digits (human-readable)
 */
// [[Rcpp::export]]
void exportAFModelInCPP(std::string genoType,
                        const std::vector<uint64_t>& genoIndex,
                        std::string afFileOutput,
                        const arma::mat& t_pcs,
                        std::string afFilePrecision) {

  g_pcs = t_pcs;
  g_N = g_pcs.n_rows;
  g_npcs = g_pcs.n_cols;

  // Precompute constant matrices
  g_X = arma::join_horiz(arma::ones(g_N), g_pcs);
  arma::mat XtXInv = arma::inv(g_X.t() * g_X);
  g_sqrtXtXInvDiag = arma::sqrt(XtXInv.diag());

  int numMarkers = genoIndex.size();
  int progressStep = std::max(1, numMarkers / 100);

  long long recordSize = 0;
  std::fstream outFileBinary;
  std::ofstream outFileText;

  // Open file handle based on format
  if (afFilePrecision == "double") {
    // Create file if needed
    {
      std::ifstream testFile(afFileOutput);
      if (!testFile.good()) {
        std::ofstream createFile(afFileOutput, std::ios::binary);
        createFile.close();
      }
    }
    recordSize = sizeof(int) + static_cast<long long>(g_npcs + 1) * sizeof(double);
    outFileBinary.open(afFileOutput, std::ios::binary | std::ios::in | std::ios::out);
    if (!outFileBinary) {
      Rcpp::stop("Failed to open output file: " + afFileOutput);
    }
  } else if (afFilePrecision == "single") {
    // Create file if needed
    {
      std::ifstream testFile(afFileOutput);
      if (!testFile.good()) {
        std::ofstream createFile(afFileOutput, std::ios::binary);
        createFile.close();
      }
    }
    recordSize = sizeof(int) + static_cast<long long>(g_npcs + 1) * sizeof(float);
    outFileBinary.open(afFileOutput, std::ios::binary | std::ios::in | std::ios::out);
    if (!outFileBinary) {
      Rcpp::stop("Failed to open output file: " + afFileOutput);
    }
  } else if (afFilePrecision == "text") {  // text
    outFileText.open(afFileOutput);
    if (!outFileText) {
      Rcpp::stop("Failed to open output file: " + afFileOutput);
    }
    // Write header
    outFileText << "Marker\tStatus";
    for (int j = 0; j <= g_npcs; ++j) {
      outFileText << "\tBeta" << j;
    }
    outFileText << "\n";

  } else {
    Rcpp::stop("Invalid afFilePrecision: " + afFilePrecision + 
               ". Must be 'double', 'single', or 'text'.");
  }
  
  // =========== Process markers ===========
  for (int i = 0; i < numMarkers; ++i) {
    if (i % progressStep == 0) {
      Rcpp::checkUserInterrupt();
      Rcpp::Rcout << "Processed " << i << " / " << numMarkers << " markers (" 
          << std::fixed << std::setprecision(1) << 100.0 * i / numMarkers << "%)" << "\r";

    }
    
    uint64_t markerIndex = genoIndex[i];
    
    // Get marker genotypes
    std::string marker, chr, ref, alt;
    uint32_t position;
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> missingIndices, nonZeroIndices;
    
    arma::vec genotypeVec = Unified_getOneMarker(
      genoType, markerIndex, ref, alt, marker, position, chr,
      altFreq, altCounts, missingRate, imputeInfo,
      true, missingIndices, false, nonZeroIndices
    );
    
    // Impute genotypes
    imputeGenoAndFlip(genotypeVec, altFreq, missingIndices, missingRate, 
                      "mean", "SPAmixPlus");
    
    // Compute AF model
    int status;
    arma::vec betas(g_npcs + 1);
    computeAFoneMarker(genotypeVec, altFreq, status, betas);
    
    // Write to file based on format
    if (afFilePrecision == "double") {
      // int32 status + double betas
      long long filePos = static_cast<long long>(markerIndex) * recordSize;
      outFileBinary.seekp(filePos, std::ios::beg);
      outFileBinary.write(reinterpret_cast<const char*>(&status), sizeof(int));
      outFileBinary.write(reinterpret_cast<const char*>(betas.memptr()), 
                           (g_npcs + 1) * sizeof(double));
    } else if (afFilePrecision == "single") {
      // int32 status + float32 betas
      long long filePos = static_cast<long long>(markerIndex) * recordSize;
      outFileBinary.seekp(filePos, std::ios::beg);
      outFileBinary.write(reinterpret_cast<const char*>(&status), sizeof(int));
      arma::fvec betas_float32 = arma::conv_to<arma::fvec>::from(betas);
      outFileBinary.write(reinterpret_cast<const char*>(betas_float32.memptr()), 
                           (g_npcs + 1) * sizeof(float));
    } else {  // text
      // TSV format: markerIndex, status, beta0, beta1, ...
      char buffer[32];
      snprintf(buffer, sizeof(buffer), "%s\t%d", marker.c_str(), status);
      outFileText << buffer;
      for (int j = 0; j <= g_npcs; ++j) {
        snprintf(buffer, sizeof(buffer), "\t%.6g", betas(j));
        outFileText << buffer;
      }
      outFileText << "\n";
    }
  }
  
  // Close file handle based on format
  if (afFilePrecision == "double" || afFilePrecision == "single") {
    outFileBinary.close();
  } else if (afFilePrecision == "text") {
    outFileText.close();
  }
  
  // For final message:
  Rcpp::Rcout << "\nCompleted processing " << numMarkers << " markers" << std::endl;
}
