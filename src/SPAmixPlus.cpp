// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map> 
#include <vector>        
#include <unordered_set>

#include "SPAmixPlus.h"

namespace SPAmixPlus {

// ==========================================
// Helper functions (copied from GRAB SPAmixPlusV4.hpp/cpp)
// ==========================================

  // The MGF of G (genotype)
  arma::vec M_G0(arma::vec t, arma::vec MAF){
    arma::vec re = pow((1 - MAF + MAF % arma::exp(t)), 2);
    return re;
  }
  
  // The first derivative of the MGF of G (genotype)
  arma::vec M_G1(arma::vec t, arma::vec MAF){
    arma::vec re = 2 * (MAF % arma::exp(t)) % (1 - MAF + MAF % arma::exp(t));
    return re;                           
  }
  
  // The second derivative of the MGF of G (genotype)
  arma::vec M_G2(arma::vec t, arma::vec MAF){
    arma::vec re = 2 * pow(MAF % arma::exp(t), 2) + 2 * (MAF % arma::exp(t)) % (1 - MAF + MAF % arma::exp(t));
    return re;
  }
  
  // The CGF of G (genotype)
  arma::vec K_G0(arma::vec t, arma::vec MAF){
    arma::vec re = arma::log(M_G0(t, MAF));
    return re;
  }
  
  // The first derivative of the CGF of G (genotype)
  arma::vec K_G1(arma::vec t, arma::vec MAF){
    arma::vec re = M_G1(t, MAF) / M_G0(t, MAF);
    return re;
  }
  
  // The second derivative of the CGF of G (genotype)
  arma::vec K_G2(arma::vec t, arma::vec MAF){
    arma::vec re = (M_G0(t, MAF) % M_G2(t, MAF) - pow(M_G1(t, MAF), 2)) / pow(M_G0(t, MAF), 2);
    return re;
  }

  // partial normal distribution approximation
  
  arma::vec Horg_H2(double t, arma::vec R, const arma::vec MAFVec)
  {
    arma::vec Horg_H2_vec(2);
    arma::vec t_R = t * R;
    arma::vec exp_tR = arma::exp(t_R);
    arma::vec MAF_exp_tR = MAFVec % exp_tR;
    arma::vec M_G0_vec = pow((1 - MAFVec + MAF_exp_tR), 2);;
    arma::vec M_G1_vec = 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec M_G2_vec = 2 * pow(MAF_exp_tR, 2) + 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec K_G0_vec = arma::log(M_G0_vec);
    arma::vec K_G2_vec = (M_G0_vec % M_G2_vec - pow(M_G1_vec, 2)) / pow(M_G0_vec, 2);
    double Horg = sum(K_G0_vec);
    double H2 = sum(pow(R, 2) % K_G2_vec);
    Horg_H2_vec.at(0) = Horg;
    Horg_H2_vec.at(1) = H2;
    return Horg_H2_vec;
  }
  
  arma::vec H1_adj_H2(double t, arma::vec R, double s, const arma::vec MAFVec)
  {
    arma::vec H1_adj_H2_vec(2);
    arma::vec t_R = t * R;
    arma::vec exp_tR = arma::exp(t_R);
    arma::vec MAF_exp_tR = MAFVec % exp_tR;
    arma::vec M_G0_vec = pow((1 - MAFVec + MAF_exp_tR), 2);;
    arma::vec M_G1_vec = 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec M_G2_vec = 2 * pow(MAF_exp_tR, 2) + 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec K_G1_vec = M_G1_vec / M_G0_vec;
    arma::vec K_G2_vec = (M_G0_vec % M_G2_vec - pow(M_G1_vec, 2)) / pow(M_G0_vec, 2);
    double H1_adj = sum(R % K_G1_vec) - s;
    double H2 = sum(pow(R, 2) % K_G2_vec);
    H1_adj_H2_vec.at(0) = H1_adj;
    H1_adj_H2_vec.at(1) = H2;
    return H1_adj_H2_vec;
  }
  
  Rcpp::List fastgetroot_K1(double t_initX,
                            const double& s,
                            const arma::vec MAF_outlier,
                            double mean_nonOutlier,
                            double var_nonOutlier,
                            const arma::vec residOutlier)
  {
    double x = t_initX, oldX;
    double K1 = 0, K2 = 0, oldK1;
    double diffX = arma::datum::inf, oldDiffX;
    bool converge = true;
    double tol = 0.001;
    int maxiter = 100;
    int iter = 0;
    
    for(iter = 0; iter < maxiter; iter ++){
      
      oldX = x;
      oldDiffX = diffX;
      oldK1 = K1;
      
      arma::vec H1_adj_H2_vec = H1_adj_H2(x, residOutlier, s, MAF_outlier);
      
      K1 = H1_adj_H2_vec.at(0) + mean_nonOutlier + var_nonOutlier * x;
      K2 = H1_adj_H2_vec.at(1) + var_nonOutlier;
      
      diffX = -1 * K1 / K2;
      
      if(!std::isfinite(K1)){
        x = arma::datum::inf;
        K2 = 0;
        break;
      }
      
      if(arma::sign(K1) != arma::sign(oldK1)){
        while(std::abs(diffX) > std::abs(oldDiffX) - tol){
          diffX = diffX / 2;
        }
      }
      
      if(std::abs(diffX) < tol) break;
      
      x = oldX + diffX;
    }
    
    if(iter == maxiter) 
      converge = false;
    
    Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                          Rcpp::Named("iter") = iter,
                                          Rcpp::Named("converge") = converge,
                                          Rcpp::Named("K2") = K2);
    return yList;
  }
  
  double GetProb_SPA_G(const arma::vec MAF_outlier, 
                       const arma::vec residOutlier, 
                       double s, 
                       bool lower_tail,
                       double mean_nonOutlier,
                       double var_nonOutlier)
  {
    double initX = 0;
    
    Rcpp::List rootList = fastgetroot_K1(initX, s, MAF_outlier, mean_nonOutlier, var_nonOutlier, residOutlier);
    double zeta = rootList["root"];
    
    arma::vec k12 = Horg_H2(zeta, residOutlier, MAF_outlier);
    double k1 = k12.at(0) + mean_nonOutlier * zeta + 0.5 * var_nonOutlier * pow(zeta, 2);
    double k2 = k12.at(1) + var_nonOutlier;
    
    double temp1 = zeta * s - k1;
    
    double w = arma::sign(zeta) * sqrt(2 * temp1);
    double v = zeta * sqrt(k2);
    
    double pval = arma::normcdf(arma::sign(lower_tail - 0.5) * (w + log(v/w) / w));
    return pval;
  }

  // Modified helper: returns beta instead of MAFest
  arma::vec logistic_regression_beta(const arma::mat& X, const arma::vec& y) {
    int n = X.n_rows;
    int p = X.n_cols;
    
    arma::mat WX_new(n, p + 1);
    arma::mat X_new = arma::join_horiz(arma::ones(n), X);
    
    arma::vec beta(p+1, arma::fill::zeros);
    double tol = 1e-6;
    int max_iter = 100;
    arma::vec mu(n);
    
    for (int i = 0; i < max_iter; i++) {
        mu = 1.0 / (1.0 + arma::exp(-X_new * beta));
        
        arma::vec W = mu % (1.0 - mu);
        arma::vec z = X_new * beta + (y - mu) / W;
        
        for(int j = 0; j < p+1; j++){
            WX_new.col(j) = X_new.col(j) % W;
        }
        
        arma::vec beta_new = arma::solve(X_new.t() * WX_new, X_new.t() * (W % z));
        
        if (arma::norm(beta_new - beta) < tol) {
          beta = beta_new;
          break;
        }
        beta = beta_new;
    }
    return beta;
  }

// ==========================================
// SPAmixPlusClass Implementation
// ==========================================

SPAmixPlusClass::SPAmixPlusClass(
  const arma::mat& t_resid,
  const arma::mat& t_PCs,
  int t_N,
  double t_SPA_Cutoff,
  const Rcpp::List& t_outlierList,
  const Rcpp::DataFrame& t_sparseGRM,
  const std::string& t_afFilePath,
  const std::string& t_afFilePrecision
) {
  
  // ==== Store AF file info ====
  m_afFilePath = t_afFilePath;
  m_afFilePrecision = t_afFilePrecision;
  
  // ==== Process sparseGRM ====
  // DataFrame with columns: id1_index, id2_index, value
  Rcpp::IntegerVector id1_indices = t_sparseGRM["id1_index"];
  Rcpp::IntegerVector id2_indices = t_sparseGRM["id2_index"];
  Rcpp::NumericVector values = t_sparseGRM["value"];
  
  int n = id1_indices.size();
  m_sparseTriplets.clear();
  m_sparseTriplets.reserve(n);
  
  for(int i = 0; i < n; ++i){
    m_sparseTriplets.emplace_back(id1_indices[i], id2_indices[i], values[i]);
  }
  
  // Standard initialization
  m_resid = t_resid;
  m_PCs = t_PCs;
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  
  m_Npheno = m_resid.n_cols;
  
  m_pvalVec.zeros(m_Npheno);
  m_zScoreVec.zeros(m_Npheno);
  m_BetaVec.zeros(m_Npheno);
  m_SVec.zeros(m_Npheno);
  m_SmeanVec.zeros(m_Npheno);
  m_VarSVec.zeros(m_Npheno);
  
  m_outlierList = t_outlierList;

  // Precompute matrices for fit_lm (from GRAB constructor)
  m_onePlusPCs = arma::join_horiz(arma::ones(t_N), t_PCs);
  arma::mat X_t = m_onePlusPCs.t();
  arma::mat XTX = X_t * m_onePlusPCs;
  arma::mat XTX_inv = arma::inv(XTX); // GRAB uses inv
  m_sqrt_XTX_inv_diag = arma::sqrt(XTX_inv.diag());

  m_diffTime1 = arma::vec(1, arma::fill::zeros);
  m_diffTime2 = arma::vec(1, arma::fill::zeros);
}

// Helper method: fit_lm (Modified to return beta if needed, but keeping signature for now)
// We add a helper to just get beta
arma::vec SPAmixPlusClass::fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues)
{
    int n = m_N;
    int k = m_PCs.n_cols;
    
    // m_onePlusPCs has K+1 columns.
    // coef will be K+1.
    arma::vec coef = arma::solve(m_onePlusPCs, g); 
    
    // We need pvalues for the selection logic, so we must calculate stats
    arma::vec fittedValues = m_onePlusPCs * coef;
    double s2 = arma::sum(arma::square(g - fittedValues)) / (n - k - 1);
    
    // Se calculation relies on diag(inv(XTX)) which is m_sqrt_XTX_inv_diag
    arma::vec se = m_sqrt_XTX_inv_diag * std::sqrt(s2);
    arma::vec t = coef / se;
    
    for(int i = 0; i < k; i++){
      // p-values for PCs (skipping intercept at index 0)
      pvalues[i] = 2 * R::pt(std::abs(t[i+1]), n-k-1, 0, 0); 
    }
    return coef;
}

arma::vec SPAmixPlusClass::fit_lm(const arma::vec& g, arma::vec& pvalues) 
{
    arma::vec coef = fit_lm_get_beta(g, pvalues);
    return m_onePlusPCs * coef;
}

SPAmixPlusClass::AFModelInfo SPAmixPlusClass::computeAFModel(arma::vec t_GVec, double t_altFreq) {
    AFModelInfo model;
    // Defaults
    model.status = 0; // Mean
    // Beta size K+1
    int k = m_PCs.n_cols;
    model.betas = arma::vec(k + 1, arma::fill::zeros);
    
    double MAC_cutoff = 20;
    double PCs_pvalue_cutoff = 0.05;
    double MAF_est_negative_ratio_cutoff = 0.1;
    
    int N = m_N;
    arma::vec g = t_GVec; 
    
    double MAC = t_altFreq * 2.0 * N; 
    arma::vec pvalues(k);
    
    if(MAC <= MAC_cutoff){
       model.status = 0;
       // For status 0, we might just rely on t_altFreq provided at runtime in Step 2, 
       // or we store it in betas[0]. 
       // The original code uses t_altFreq from input arg, not computed from G.
       // However, to be self-contained, let's store it. But t_altFreq varies per marker.
       // The 'model' is specific to a marker. So we can store t_altFreq.
       // But wait, t_altFreq is passed to getAFFromModel(model, t_altFreq).
       // So we don't strictly need to store it if we pass it back.
       // Just set status=0.
    }else{
      arma::vec coef_lm = fit_lm_get_beta(g, pvalues);
      arma::vec fit = m_onePlusPCs * coef_lm; // Compute fit to check bounds
      fit = fit / 2.0; // Scale to 0-1
      
      arma::uvec posZero = arma::find(fit < 0);
      arma::uvec posOne = arma::find(fit > 1);
      
      int nError = posZero.n_elem + posOne.n_elem; 
      double propError = (double)nError / N;
      
      if(propError < MAF_est_negative_ratio_cutoff){
        // Use Linear Model
        model.status = 1;
        model.betas = coef_lm;
      }else{
        arma::uvec posSigPCs = arma::find(pvalues < PCs_pvalue_cutoff);

        if(posSigPCs.n_elem == 0){
           model.status = 0; // Fallback to mean
        }else{
           // Logistic Regression on significant PCs
           arma::mat sigPCs = m_PCs.cols(posSigPCs);
           
           // Impute g for logistic
           arma::vec g0(N, arma::fill::zeros);  
           arma::uvec posg12 = arma::find(g > 0.5);
           g0.elem(posg12).fill(1);
           
           double MAC_after = sum(g0);
           if(MAC_after <= MAC_cutoff){
             model.status = 0; // Fallback
           }else{
             model.status = 2; 
             // Run logistic regression to get coefficients for SUBSET
             arma::vec sub_beta = logistic_regression_beta(sigPCs, g0);
             
             // Map sub_beta back to full beta vector
             // sub_beta has length posSigPCs.n_elem + 1 (intercept)
             model.betas(0) = sub_beta(0); // Intercept
             for(unsigned int i=0; i < posSigPCs.n_elem; ++i){
                // posSigPCs contains indices into m_PCs (0 to K-1).
                // m_onePlusPCs cols are 1 to K.
                // model.betas index matches m_onePlusPCs index.
                int pc_index = posSigPCs(i);
                model.betas(pc_index + 1) = sub_beta(i + 1);
             }
           }
        }
      }
    }
    return model;
}

arma::vec SPAmixPlusClass::getAFFromModel(AFModelInfo t_model, double t_altFreq) {
    if(t_model.status == 0){
        return arma::vec(m_N, arma::fill::value(t_altFreq));
    }
    
    if(t_model.status == 1){ 
       // Linear
       arma::vec fit = m_onePlusPCs * t_model.betas;
       fit = fit / 2.0;
       
       // Apply clamping logic (same as original)
       arma::uvec posZero = arma::find(fit < 0);
       arma::uvec posOne = arma::find(fit > 1);
       fit.elem(posZero).fill(0.0);
       fit.elem(posOne).fill(1.0);
       
       return fit;
    }
    
    if(t_model.status == 2){
       // Logistic
       // We have full betas (with zeros for insignificant PCs)
       // Calculate X * Beta
       arma::vec linear_pred = m_onePlusPCs * t_model.betas;
       
       // Sigmoid
       arma::vec mu = 1.0 / (1.0 + arma::exp(-linear_pred));
       
       // Convert to MAF
       return 1.0 - arma::sqrt(1.0 - mu);
    }
    
    // Default
    return arma::vec(m_N, arma::fill::value(t_altFreq));
}

arma::vec SPAmixPlusClass::getMAFest(arma::vec t_GVec, double t_altFreq) {
    // Legacy wrapper
    AFModelInfo model = computeAFModel(t_GVec, t_altFreq);
    return getAFFromModel(model, t_altFreq);
}

// Overload/Update getMarkerPval
double SPAmixPlusClass::getMarkerPval(arma::vec t_GVec, double t_altFreq) {
    
    arma::vec AFVec = getMAFest(t_GVec, t_altFreq);
    m_MAFVec = AFVec;
    
    // Variance component based on estimated MAF
    arma::vec GVarVec = 2.0 * AFVec % (1.0 - AFVec);

    for(int i = 0; i < m_Npheno; i++){
        Rcpp::List tempOutlierList = m_outlierList[i];
        
        arma::uvec posValue = tempOutlierList["posValue"];
        arma::uvec posOutlier = tempOutlierList["posOutlier"];
        arma::uvec posNonOutlier = tempOutlierList["posNonOutlier"];
        
        arma::vec resid = tempOutlierList["resid"]; // Note: this name might mask method arg if not careful, but ok here
        arma::vec residOutlier = tempOutlierList["residOutlier"];
        arma::vec residNonOutlier = tempOutlierList["residNonOutlier"];
        arma::vec resid2NonOutlier = tempOutlierList["resid2NonOutlier"];

        // resid is subsetted R vector (length = posValue.n_elem)
        arma::vec R_subset = resid; 
        arma::vec GVar_subset = GVarVec.elem(posValue);
        
        // R_new for sparse variance calculation
        arma::vec R_new = R_subset % arma::sqrt(GVar_subset);
        
        double VarS = calculateSparseVariance(R_new, posValue);
        
        double S = arma::sum(t_GVec.elem(posValue) % R_subset);
        double S_mean = 2.0 * arma::sum(R_subset % AFVec.elem(posValue));
        double zScore = (S - S_mean) / std::sqrt(VarS);
        
        m_zScoreVec.at(i) = zScore;
        m_BetaVec.at(i) = (S - S_mean) / VarS;
        m_SVec.at(i) = S;
        m_SmeanVec.at(i) = S_mean;
        m_VarSVec.at(i) = VarS;
        
        if(std::abs(zScore) < m_SPA_Cutoff){
            m_pvalVec.at(i) = arma::normcdf(-1.0*std::abs(zScore))*2.0;
            continue;
        }

        // SPA adjustment
        double S_var_SPAmix = arma::sum(arma::square(R_subset) % GVar_subset);
        double Var_ratio = S_var_SPAmix / VarS;
        
        double S_new = S * std::sqrt(Var_ratio);
        double S_mean_new = S_mean * std::sqrt(Var_ratio);
        
        double S_upper = std::max(S_new, 2.0*S_mean_new - S_new); 
        double S_lower = std::min(S_new, 2.0*S_mean_new - S_new);
        
        arma::vec MAF_outlier = AFVec.elem(posOutlier);
        
        double mean_nonOutlier = arma::sum(residNonOutlier % AFVec.elem(posNonOutlier)) * 2.0; 
        double var_nonOutlier = arma::sum(resid2NonOutlier % AFVec.elem(posNonOutlier) % (1.0 - AFVec.elem(posNonOutlier))) * 2.0;

        double pval1 = GetProb_SPA_G(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
        double pval2 = GetProb_SPA_G(MAF_outlier, residOutlier, S_lower, true, mean_nonOutlier, var_nonOutlier);
        
        m_pvalVec.at(i) = pval1 + pval2;
    }
    return 0.0;
}

double SPAmixPlusClass::getMarkerPvalFromModel(arma::vec t_GVec, AFModelInfo t_model, double t_altFreq) {
    arma::vec AFVec = getAFFromModel(t_model, t_altFreq);
    m_MAFVec = AFVec;
    
    // Variance component based on estimated MAF
    arma::vec GVarVec = 2.0 * AFVec % (1.0 - AFVec);

    for(int i = 0; i < m_Npheno; i++){
        Rcpp::List tempOutlierList = m_outlierList[i];
        
        arma::uvec posValue = tempOutlierList["posValue"];
        arma::uvec posOutlier = tempOutlierList["posOutlier"];
        arma::uvec posNonOutlier = tempOutlierList["posNonOutlier"];
        
        arma::vec resid = tempOutlierList["resid"]; 
        arma::vec residOutlier = tempOutlierList["residOutlier"];
        arma::vec residNonOutlier = tempOutlierList["residNonOutlier"];
        arma::vec resid2NonOutlier = tempOutlierList["resid2NonOutlier"];

        // resid is subsetted R vector (length = posValue.n_elem)
        arma::vec R_subset = resid; 
        arma::vec GVar_subset = GVarVec.elem(posValue);
        
        // R_new for sparse variance calculation
        arma::vec R_new = R_subset % arma::sqrt(GVar_subset);
        
        double VarS = calculateSparseVariance(R_new, posValue);
        
        double S = arma::sum(t_GVec.elem(posValue) % R_subset);
        double S_mean = 2.0 * arma::sum(R_subset % AFVec.elem(posValue));
        double zScore = (S - S_mean) / std::sqrt(VarS);
        
        // Store statistics
        m_zScoreVec.at(i) = zScore;
        m_BetaVec.at(i) = (S - S_mean) / VarS;
        m_SVec.at(i) = S;
        m_SmeanVec.at(i) = S_mean;
        m_VarSVec.at(i) = VarS;
        
        if(std::abs(zScore) < m_SPA_Cutoff){
            // Normal approximation
            m_pvalVec.at(i) = 2.0 * R::pnorm(std::abs(zScore), 0.0, 1.0, 0, 0); 
            continue;
        }
        
        // SPA adjustment
        double S_var_SPAmix = arma::sum(arma::square(R_subset) % GVar_subset);
        double Var_ratio = S_var_SPAmix / VarS;
        
        double S_new = S * std::sqrt(Var_ratio);
        double S_mean_new = S_mean * std::sqrt(Var_ratio);
        
        double S_upper = std::max(S_new, 2.0*S_mean_new - S_new); 
        double S_lower = std::min(S_new, 2.0*S_mean_new - S_new);
        
        arma::vec MAF_outlier = AFVec.elem(posOutlier);
        
        double mean_nonOutlier = arma::sum(residNonOutlier % AFVec.elem(posNonOutlier)) * 2.0; 
        double var_nonOutlier = arma::sum(resid2NonOutlier % AFVec.elem(posNonOutlier) % (1.0 - AFVec.elem(posNonOutlier))) * 2.0;

        double pval1 = GetProb_SPA_G(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
        double pval2 = GetProb_SPA_G(MAF_outlier, residOutlier, S_lower, true, mean_nonOutlier, var_nonOutlier);
        
        m_pvalVec.at(i) = pval1 + pval2;
    }
    return 0.0;
}

double SPAmixPlusClass::calculateSparseVariance(const arma::vec& R_new, const arma::uvec& posValue) {
    double covSum = 0.0;
    
    // Use unordered_set for O(1) lookups
    std::unordered_set<int> validIndices;
    for (auto idx : posValue) {
        validIndices.insert(static_cast<int>(idx));
    }
    
    for (const auto& triplet : m_sparseTriplets) {
      int i = std::get<0>(triplet);
      int j = std::get<1>(triplet);
      
      // Only include if BOTH individuals i and j are in the current subset (posValue)
      if (validIndices.count(i) && validIndices.count(j)) {
        double grmValue = std::get<2>(triplet);
        covSum += grmValue * R_new(i) * R_new(j);
      }
    }
    return 2.0 * covSum - arma::sum(arma::square(R_new));
  }

}


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
