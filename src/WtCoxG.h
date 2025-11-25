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

namespace WtCoxG {
class WtCoxGClass
{
private:
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

    std::vector<MarkerInfo> m_markerInfoVec;  // Vector of marker information
    arma::vec m_R;                       // Residuals vector
    arma::vec m_w;                       // Weights vector
    double m_cutoff;                     // batch effect p-value cutoff
    double m_SPA_Cutoff;                 // SPA cutoff

public:

    WtCoxGClass(
        const arma::vec& t_R,
        const arma::vec& t_w,
        const double t_cutoff,
        const double t_SPA_Cutoff
    )
        : m_R(t_R),
        m_w(t_w),
        m_cutoff(t_cutoff),
        m_SPA_Cutoff(t_SPA_Cutoff)
    {
        // Constructor body (if needed)
    }

    // reset m_markerInfoVec with new data from R data.frame
    void updateMarkerInfo(const Rcpp::DataFrame& t_mergeGenoInfo) {
        m_markerInfoVec.clear();

        Rcpp::NumericVector AF_ref = t_mergeGenoInfo["AF_ref"];
        Rcpp::NumericVector AN_ref = t_mergeGenoInfo["AN_ref"];
        Rcpp::NumericVector TPR = t_mergeGenoInfo["TPR"];
        Rcpp::NumericVector sigma2 = t_mergeGenoInfo["sigma2"];
        Rcpp::NumericVector pvalue_bat = t_mergeGenoInfo["pvalue_bat"];
        Rcpp::NumericVector w_ext = t_mergeGenoInfo["w.ext"];
        Rcpp::NumericVector var_ratio_w0 = t_mergeGenoInfo["var.ratio.w0"];
        Rcpp::NumericVector var_ratio_int = t_mergeGenoInfo["var.ratio.int"];
        Rcpp::NumericVector var_ratio_ext = t_mergeGenoInfo["var.ratio.ext"];

        int n = AF_ref.size();
        m_markerInfoVec.reserve(n);

        for (int i = 0; i < n; ++i) {
            m_markerInfoVec.emplace_back(
                AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
                w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
            );
        }
    }
    
    // Method called in Main.cpp to get p-values for a given marker
    arma::vec getpvalVec(const arma::vec& GVec, const int i) {
        arma::vec result(2);
        const MarkerInfo& info = m_markerInfoVec[i];

        // Test with external reference - match R parameter mapping exactly
        result(0) = WtCoxG_test_cpp(
            GVec, m_R, m_w,
            info.pvalue_bat,           // p_bat
            info.TPR,                  // TPR
            info.sigma2,               // sigma2
            info.w_ext,                // b = w.ext
            info.var_ratio_int,        // var_ratio_int (not used in ext test)
            info.var_ratio_w0,         // var_ratio_w0
            info.var_ratio_w0,         // var_ratio_w1 = var_ratio_w0
            info.var_ratio_ext,        // var_ratio0 = var_ratio_ext
            info.var_ratio_ext,        // var_ratio1 = var_ratio_ext
            info.AF_ref,               // mu.ext
            info.AN_ref / 2.0,         // n.ext
            m_cutoff);                 // p_cut
        

        // Test without external reference - match R parameter mapping exactly
        result(1) = WtCoxG_test_cpp(
            GVec, m_R, m_w,
            info.pvalue_bat,           // p_bat
            arma::datum::nan,          // TPR = NA
            arma::datum::nan,          // sigma2 = NA
            0.0,                       // b = 0 (default)
            info.var_ratio_int,        // var_ratio_int
            1.0,                       // var_ratio_w0 = 1 (default)
            1.0,                       // var_ratio_w1 = 1 (default)
            1.0,                       // var_ratio0 = 1 (default)
            1.0,                       // var_ratio1 = 1 (default)
            arma::datum::nan,          // mu.ext = NA
            arma::datum::nan,          // n.ext = NA
            m_cutoff);                 // p_cut

        return result;
    }

private:
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
    
    arma::vec SPA_G_one_SNP_homo_cpp(const arma::vec& g_input, const arma::vec& R,
                                    double mu_ext, double n_ext, double b,
                                    double sigma2, double var_ratio, double Cutoff,
                                    double missing_cutoff, double min_mac);

    double WtCoxG_test_cpp(const arma::vec& g_input, const arma::vec& R, const arma::vec& w,
                        double p_bat, double TPR, double sigma2, double b,
                        double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
                        double var_ratio0, double var_ratio1, double mu_ext,
                        double n_ext, double p_cut);
};

} // namespace WtCoxG


#endif
