#ifndef WTCOXG_H
#define WTCOXG_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using namespace Rcpp;
using namespace arma;

namespace WtCoxG
{
    // Structure to hold marker information
    struct MarkerInfo
    {
        double AF_ref;
        double AN_ref;
        double TPR;
        double sigma2;
        double pvalue_bat;
        double w_ext;
        double var_ratio_w0;
        double var_ratio_int;
        double var_ratio_ext;
        
        // Default constructor
        MarkerInfo() = default;
        
        // Constructor with all parameters
        MarkerInfo(double af_ref, double an_ref, double tpr, double sig2,
                  double pval_bat, double w_ext_val, double var_ratio_w0_val,
                  double var_ratio_int_val, double var_ratio_ext_val)
            : AF_ref(af_ref), AN_ref(an_ref), TPR(tpr), sigma2(sig2),
              pvalue_bat(pval_bat), w_ext(w_ext_val), var_ratio_w0(var_ratio_w0_val),
              var_ratio_int(var_ratio_int_val), var_ratio_ext(var_ratio_ext_val) {}
    };

    class WtCoxGClass
    {
    private:
        // -------------------- Member Variables --------------------
        arma::vec m_R;                       // Residuals vector
        arma::vec m_w;                       // Weights vector
        DataFrame m_mergeGenoInfo;           // Merged genotype information (updated per chunk)
        std::string m_imputeMethod;          // Imputation method
        double m_cutoff;                     // Cutoff threshold
        std::vector<MarkerInfo> m_markerInfoVec;  // Vector of marker information
        arma::vec m_pvalVec;                 // Vector to store p-values

        // -------------------- Private Helper Methods --------------------
        // Helper method to extract marker info from DataFrame
        void extractMarkerInfo() {
            NumericVector AF_ref = m_mergeGenoInfo["AF_ref"];
            NumericVector AN_ref = m_mergeGenoInfo["AN_ref"];
            NumericVector TPR = m_mergeGenoInfo["TPR"];
            NumericVector sigma2 = m_mergeGenoInfo["sigma2"];
            NumericVector pvalue_bat = m_mergeGenoInfo["pvalue_bat"];
            NumericVector w_ext = m_mergeGenoInfo["w.ext"];
            NumericVector var_ratio_w0 = m_mergeGenoInfo["var.ratio.w0"];
            NumericVector var_ratio_int = m_mergeGenoInfo["var.ratio.int"];
            NumericVector var_ratio_ext = m_mergeGenoInfo["var.ratio.ext"];

            m_markerInfoVec.clear();
            m_markerInfoVec.reserve(AF_ref.size());

            for (int i = 0; i < AF_ref.size(); ++i) {
                m_markerInfoVec.emplace_back(
                    AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
                    w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
                );
            }
        }

    public:
        // -------------------- Constructor --------------------
        WtCoxGClass(const arma::vec& R,
                    const arma::vec& w,
                    const std::string& imputeMethod = "none",
                    double cutoff = 0.1);

        // -------------------- Main Public Interface --------------------
        // Method to calculate p-values for a single marker using marker-specific info
        arma::vec getpvalVec(const arma::vec& GVec, const MarkerInfo& markerInfo);
        
        // Method to calculate p-values for a single marker by index (for chunk-based processing)
        // This method looks up marker info from the stored DataFrame subset
        arma::vec getpvalVec(const arma::vec& GVec, int i);
        
        // Method to calculate p-values for a single marker by genoIndex to avoid indexing mismatch
        arma::vec getpvalVecByGenoIndex(const arma::vec& GVec, uint64_t genoIndex);
        
        // Method to update marker information for current chunk
        void updateMarkerInfo(const DataFrame& mergeGenoInfo_subset);

        // -------------------- Public Getters --------------------
        const arma::vec& getResiduals() const { return m_R; }
        const arma::vec& getWeights() const { return m_w; }
        const DataFrame& getCurrentMarkerInfo() const { return m_mergeGenoInfo; }
        const std::string& getImputeMethod() const { return m_imputeMethod; }
        double getCutoff() const { return m_cutoff; }
        size_t getNumMarkers() const { return m_markerInfoVec.size(); }
    };

    // -------------------- Standalone Functions --------------------
    // Core utility functions - inline implementations for performance
    inline arma::vec impute_missing(const arma::vec& g) {
        arma::vec g_imputed = g;
        arma::uvec missing_idx = arma::find_nonfinite(g);
        
        if (missing_idx.n_elem > 0) {
            arma::uvec non_missing_idx = arma::find_finite(g);
            if (non_missing_idx.n_elem > 0) {
                double mean_val = arma::mean(g.elem(non_missing_idx));
                g_imputed.elem(missing_idx).fill(mean_val);
            }
        }
        
        return g_imputed;
    }

    // Boost-based statistical functions (more accurate than simple approximations)
    inline double pnorm_boost(double x, double mean = 0.0, double sd = 1.0, bool lower_tail = true, bool log_p = false) {
        boost::math::normal dist(mean, sd);
        double result = boost::math::cdf(dist, x);
        if (!lower_tail) result = 1.0 - result;
        if (log_p) result = std::log(result);
        return result;
    }

    inline double qnorm_boost(double p, double mean = 0.0, double sd = 1.0, bool lower_tail = true, bool log_p = false) {
        if (log_p) p = std::exp(p);
        if (!lower_tail) p = 1.0 - p;
        boost::math::normal dist(mean, sd);
        return boost::math::quantile(dist, p);
    }

    inline double qchisq_boost(double p, double df, bool lower_tail = true, bool log_p = false) {
        if (log_p) p = std::exp(p);
        if (!lower_tail) p = 1.0 - p;
        boost::math::chi_squared dist(df);
        return boost::math::quantile(dist, p);
    }

    // MGF and CGF functions - inline for performance in tight loops
    inline double M_G0_cpp(double t, double MAF) {
        return std::pow(1.0 - MAF + MAF * std::exp(t), 2.0);
    }

    inline double M_G1_cpp(double t, double MAF) {
        return 2.0 * (MAF * std::exp(t)) * (1.0 - MAF + MAF * std::exp(t));
    }

    inline double M_G2_cpp(double t, double MAF) {
        double maf_exp_t = MAF * std::exp(t);
        return 2.0 * maf_exp_t * maf_exp_t + 2.0 * maf_exp_t * (1.0 - MAF + maf_exp_t);
    }

    inline double K_G0_cpp(double t, double MAF) {
        return std::log(M_G0_cpp(t, MAF));
    }

    inline double K_G1_cpp(double t, double MAF) {
        return M_G1_cpp(t, MAF) / M_G0_cpp(t, MAF);
    }

    inline double K_G2_cpp(double t, double MAF) {
        double m0 = M_G0_cpp(t, MAF);
        double m1 = M_G1_cpp(t, MAF);
        double m2 = M_G2_cpp(t, MAF);
        return (m0 * m2) / (m0 * m0) - std::pow(m1 / m0, 2.0);
    }

    // Function declarations for implementations in WtCoxG.cpp
    double find_root_brent(std::function<double(double)> f, double a, double b, double tol = 1e-6);
    
    double H_org_cpp(double t, const arma::vec& R, double MAF, double n_ext,
                     double N_all, double sumR, double var_mu_ext,
                     double g_var_est, double meanR, double b);
    
    double H1_adj_cpp(double t, const arma::vec& R, double s, double MAF,
                      double n_ext, double N_all, double sumR, double var_mu_ext,
                      double g_var_est, double meanR, double b);
    
    double H2_cpp(double t, const arma::vec& R, double MAF, double n_ext,
                  double N_all, double sumR, double var_mu_ext,
                  double g_var_est, double meanR, double b);
    
    double GetProb_SPA_G_cpp(double MAF, const arma::vec& R, double s, double n_ext,
                             double N_all, double sumR, double var_mu_ext,
                             double g_var_est, double meanR, double b, bool lower_tail);
    
    // Standalone WtCoxG functions (implementations in WtCoxG.cpp)
    arma::vec SPA_G_one_SNP_homo_cpp(const arma::vec& g_input, const arma::vec& R,
                                     double mu_ext, double n_ext, double b,
                                     double sigma2, double var_ratio, double Cutoff,
                                     double missing_cutoff, double min_mac);

}  // namespace WtCoxG


double WtCoxG_test_cpp(const arma::vec& g_input, const arma::vec& R, const arma::vec& w,
                       double p_bat, double TPR, double sigma2, double b,
                       double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
                       double var_ratio0, double var_ratio1, double mu_ext,
                       double n_ext, double p_cut);

#endif  // WTCOXG_HPP
