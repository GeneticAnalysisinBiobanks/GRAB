// [[Rcpp::depends(RcppArmadillo)]]
#include "WtCoxG.h"
#include <cmath>
#include <algorithm>
#include <functional>

namespace WtCoxG {

double WtCoxGClass::find_root_brent(std::function<double(double)> f, double a, double b, double tol) {

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0) {
        throw std::runtime_error("Root not bracketed");
    }

    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c = a;
    double fc = fa;

    const int max_iter = std::min(50, std::max(15, static_cast<int>(20.0 / tol)));
    tol = std::max(tol, 1e-9);

    double last_improvement = std::abs(b - a);
    int stagnant_iterations = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        double current_gap = std::abs(b - a);
        if (current_gap < tol || std::abs(fb) < tol * 10) {
            return b;
        }

        if (current_gap >= 0.95 * last_improvement) {
            stagnant_iterations++;
            if (stagnant_iterations > 3)
                break;
        } else {
            stagnant_iterations = 0;
        }
        last_improvement = current_gap;

        double s;
        if (fa != fc && fb != fc && std::abs(fa - fc) > 1e-15 && 
            std::abs(fb - fc) > 1e-15) {
            double denom1 = (fa - fb) * (fa - fc);
            double denom2 = (fb - fa) * (fb - fc);
            double denom3 = (fc - fa) * (fc - fb);

            if (std::abs(denom1) > 1e-12 && std::abs(denom2) > 1e-12 && 
                std::abs(denom3) > 1e-12) {
                s = a * fb * fc / denom1 + b * fa * fc / denom2 + 
                    c * fa * fb / denom3;
            } else {
                s = b - fb * (b - a) / (fb - fa);
            }
        } else {
            s = b - fb * (b - a) / (fb - fa);
        }

        if (s <= std::min(a, b) || s >= std::max(a, b)) {
            s = (a + b) / 2.0;
        }

        double fs = f(s);
        c = b;
        fc = fb;

        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    return b;
}

// CGF functions for score test statistic
double WtCoxGClass::H_org_cpp(
    double t, const arma::vec& R, double MAF, double n_ext,
    double N_all, double sumR, double var_mu_ext,
    double g_var_est, double meanR, double b
) {

    double mu_adj = -2.0 * b * sumR * MAF;
    double var_adj = 4.0 * b * b * sumR * sumR * var_mu_ext;
    
    double result = 0.0;
    for (size_t i = 0; i < R.n_elem; ++i) {
        result += K_G0_cpp(t * (R(i) - (1.0 - b) * meanR), MAF);
    }
    
    return result + mu_adj * t + var_adj * t * t / 2.0;
}

double WtCoxGClass::H1_adj_cpp(
    double t, const arma::vec& R, double s, double MAF,
    double n_ext, double N_all, double sumR, double var_mu_ext,
    double g_var_est, double meanR, double b
) {

    double mu_adj = -2.0 * b * sumR * MAF;
    double var_adj = 4.0 * b * b * sumR * sumR * var_mu_ext;
    
    double result = 0.0;
    for (size_t i = 0; i < R.n_elem; ++i) {
        double R_adj = R(i) - (1.0 - b) * meanR;
        result += R_adj * K_G1_cpp(t * R_adj, MAF);
    }
    
    return result + mu_adj + var_adj * t - s;
}

double WtCoxGClass::H2_cpp(
    double t, const arma::vec& R, double MAF, double n_ext,
    double N_all, double sumR, double var_mu_ext,
    double g_var_est, double meanR, double b
) {

    double var_adj = n_ext * std::pow(sumR / N_all, 2.0) * 2.0 * MAF * (1.0 - MAF);
    
    double result = 0.0;
    for (size_t i = 0; i < R.n_elem; ++i) {
        double R_adj = R(i) - (1.0 - b) * meanR;
        result += R_adj * R_adj * K_G2_cpp(t * R_adj, MAF);
    }
    
    return result + var_adj;
}

// SPA probability function
double WtCoxGClass::GetProb_SPA_G_cpp(
    double MAF, const arma::vec& R, double s, double n_ext,
    double N_all, double sumR, double var_mu_ext,
    double g_var_est, double meanR, double b, bool lower_tail
) {
    // Match R logic exactly: use uniroot with extendInt = "yes"
    auto h1_func = [&](double t) {
        return H1_adj_cpp(t, R, s, MAF, n_ext, N_all, sumR, var_mu_ext, 
                            g_var_est, meanR, b);
    };

    // Find root using extended interval approach like R's uniroot with extendInt = "yes"
    double zeta;
    try {
        // Start with initial interval [-1, 1] and extend if needed
        double a = -1.0, b_bound = 1.0;
        double fa = h1_func(a);
        double fb = h1_func(b_bound);
        
        // If root is not bracketed, extend the interval (like extendInt = "yes")
        if (fa * fb > 0) {
            // Try extending the interval
            double factor = 2.0;
            int max_extend = 10;
            for (int i = 0; i < max_extend; ++i) {
                if (std::abs(fa) < std::abs(fb)) {
                    // Extend left
                    a = a * factor;
                    fa = h1_func(a);
                } else {
                    // Extend right
                    b_bound = b_bound * factor;
                    fb = h1_func(b_bound);
                }
                if (fa * fb <= 0) break;
            }
        }
        
        // Find root using Brent's method
        zeta = find_root_brent(h1_func, a, b_bound, 1e-8);
    } catch (...) {
        // If root finding fails, return NaN (don't use fallback like before)
        Rcpp::Rcout << "    Root finding failed, returning NaN" << std::endl;
        return arma::datum::nan;
    }

    // Calculate k1 and k2
    double k1 = H_org_cpp(zeta, R, MAF, n_ext, N_all, sumR, var_mu_ext, 
                            g_var_est, meanR, b);
    double k2 = H2_cpp(zeta, R, MAF, n_ext, N_all, sumR, var_mu_ext, 
                        g_var_est, meanR, b);

    double temp1 = zeta * s - k1;

    // Calculate pval
    // R: w <- sign(zeta) * (2 * temp1)^{1/2}
    double w = (zeta >= 0 ? 1.0 : -1.0) * std::sqrt(2.0 * temp1);
    // R: v <- zeta * (k2)^{1/2}
    double v = zeta * std::sqrt(k2);
    // R: pval <- pnorm(w + 1 / w * log(v / w), lower.tail = lower.tail)
    double pval = pnorm_boost(w + (1.0 / w) * std::log(v / w), 0.0, 1.0, lower_tail, false);

    return pval;
}

// -------------------- Standalone Function Implementations --------------------

arma::vec WtCoxGClass::SPA_G_one_SNP_homo_cpp(
    const arma::vec& g_input, const arma::vec& R,
    double mu_ext, double n_ext, double b,
    double sigma2, double var_ratio, double Cutoff,
    double missing_cutoff, double min_mac
) {
    // Impute missing values
    arma::vec g = impute_missing(g_input);

    // Check missing rate
    arma::uvec missing_idx = arma::find_nonfinite(g_input);
    double missing_rate = static_cast<double>(missing_idx.n_elem) / g_input.n_elem;

    if (std::isnan(mu_ext))
    {
        mu_ext = 0.0;
        n_ext = 0.0;
    }

    double sum_g = arma::sum(g);
    double sum_2_minus_g = arma::sum(2.0 - g);

    if (sum_g < min_mac || sum_2_minus_g < min_mac || missing_rate > missing_cutoff)
    {
        return arma::vec({arma::datum::nan, arma::datum::nan});
    }

    // Score statistic - pre-compute common values
    double N = g.n_elem;
    double mu_int = arma::mean(g) / 2.0;
    double MAF = (1.0 - b) * mu_int + b * mu_ext;
    double sumR = arma::sum(R);
    double N_all = N + n_ext;
    double S = arma::sum(R % (g - 2.0 * MAF));
    S = S / var_ratio;

    // Estimated variance
    double g_var_est = 2.0 * MAF * (1.0 - MAF);
    double var_mu_ext = (n_ext == 0.0) ? 0.0 : (MAF * (1.0 - MAF) / (2.0 * n_ext) + sigma2);

    double meanR = arma::mean(R);
    arma::vec R_adj = R - (1.0 - b) * meanR; // Pre-compute R_adj once
    double S_var = arma::sum(R_adj % R_adj) * g_var_est + 4.0 * b * b * sumR * sumR * var_mu_ext;

    double z = S / std::sqrt(S_var);

    if (std::abs(z) < Cutoff)
    {
        // Use Boost pnorm for accuracy and performance
        double pval_norm = 2.0 * pnorm_boost(-std::abs(z), 0.0, 1.0, true, false);
        pval_norm = std::min(1.0, pval_norm); // Ensure p-value doesn't exceed 1.0
        return arma::vec({pval_norm, pval_norm});
    }
    else
    {
        double pval1 = GetProb_SPA_G_cpp(MAF, R, std::abs(S), n_ext, N_all, sumR,
                                            var_mu_ext, g_var_est, meanR, b, false);
        double pval2 = GetProb_SPA_G_cpp(MAF, R, -std::abs(S), n_ext, N_all, sumR,
                                            var_mu_ext, g_var_est, meanR, b, true);
        double pval_spa = pval1 + pval2;
        pval_spa = std::min(1.0, pval_spa); // Ensure p-value doesn't exceed 1.0
        // Use Boost pnorm for accuracy and performance
        double pval_norm = 2.0 * pnorm_boost(-std::abs(z), 0.0, 1.0, true, false);
        pval_norm = std::min(1.0, pval_norm); // Ensure p-value doesn't exceed 1.0

        return arma::vec({pval_spa, pval_norm});
    }
}


// Function calculates p-value
double WtCoxGClass::WtCoxG_test_cpp(
    const arma::vec& g_input, const arma::vec& R, const arma::vec& w,
    double p_bat, double TPR, double sigma2, double b,
    double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
    double var_ratio0, double var_ratio1, double mu_ext,
    double n_ext, double p_cut
) {
    // Step 1: Imputation of missing SNP
    arma::vec g = g_input;
    arma::uvec na_indices = arma::find_nonfinite(g);
    double missing_rate = static_cast<double>(na_indices.n_elem) / g.n_elem;
    
    if (missing_rate != 0.0) {
        arma::uvec valid_indices = arma::find_finite(g);
        if (valid_indices.n_elem > 0) {
            double mean_val = arma::mean(g.elem(valid_indices));
            g.elem(na_indices).fill(mean_val);
        }
    }
    
    // Step 2: If external MAF is unavailable, return early
    if (std::isnan(mu_ext)) {
        // In R: p.con <- SPA_G.one.SNP_homo(g = g, R = R, mu.ext = NA, n.ext = 0, sigma2 = 0, var.ratio = var.ratio.int)[1]
        // The R SPA function sets mu.ext = 0, n.ext = 0 when mu.ext is NA
        // IMPORTANT: The "ext" case defaults to var.ratio.int = 1, while "noext" case passes the actual var.ratio.int value
        // Distinguish the cases by checking if TPR/sigma2 are NaN (noext) or have values (ext)
        double var_ratio_to_use;
        if (std::isnan(TPR) && std::isnan(sigma2)) {
            // This is the "noext" case - use the passed var_ratio_int
            var_ratio_to_use = var_ratio_int;
        } else {
            // This is the "ext" case - use default value of 1.0 (like R)
            var_ratio_to_use = 1.0;
        }
        arma::vec spa_result = SPA_G_one_SNP_homo_cpp(
            g, R, 0.0, 0.0, 0.0, 0.0, var_ratio_to_use, 2.0, 0.15, 10.0);
        return spa_result(0);
    }
    
    // Step 3: Early return conditions
    double sum_g = arma::sum(g);
    double sum_2_minus_g = arma::sum(2.0 - g);
    
    if (p_bat < p_cut || std::isnan(p_bat) || sum_g < 10 || sum_2_minus_g < 10) {
        return arma::datum::nan;
    }
    
    // Step 4: Main computation starts here
    double meanR = arma::mean(R);
    double sumR = arma::sum(R);
    double mu_int = arma::mean(g) / 2.0;
    // int N = g.n_elem;
    
    // Step 5: Calculate core variables
    double mu = (1.0 - b) * mu_int + b * mu_ext;
    double S = arma::sum(R % (g - 2.0 * mu));
    
    arma::vec w1 = w / (2.0 * arma::sum(w));
    
    double var_mu_ext = mu * (1.0 - mu) / (2.0 * n_ext);
    double var_Sbat = arma::sum(w1 % w1) * 2.0 * mu * (1.0 - mu) + var_mu_ext;
    
    // Step 6: Calculate bounds
    double qnorm_val = qnorm_boost(1.0 - p_cut / 2.0, 0.0, 1.0, true, false);
    double lb = -qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);
    double ub = qnorm_val * std::sqrt(var_Sbat) * std::sqrt(var_ratio_w0);
    
    // Step 7: Calculate denominator
    double c = pnorm_boost(ub / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
    double d = pnorm_boost(lb / std::sqrt(var_ratio_w1), 0.0, std::sqrt(var_Sbat + sigma2), true, true);
    double p_deno = TPR * (std::exp(d) * (std::exp(c - d) - 1.0)) + (1.0 - TPR) * (1.0 - p_cut);
    
    // Step 8: sigma2 = 0 case
    arma::vec spa_result_s0 = SPA_G_one_SNP_homo_cpp(
        g, R, mu_ext, n_ext, b, 0.0, var_ratio0, 2.0, 0.15, 10.0);
    double p_spa_s0 = spa_result_s0(0);
    double qchisq_val = qchisq_boost(p_spa_s0, 1.0, false, false);
    double var_S = S * S / var_ratio0 / qchisq_val;
    
    // Step 9: Calculate covariance matrix components
    arma::vec R_minus_factor = R - (1.0 - b) * meanR;
    double var_int = arma::sum(R_minus_factor % R_minus_factor) * 2.0 * mu * (1.0 - mu);
    
    double cov_Sbat_S = arma::sum(w1 % R_minus_factor) * 2.0 * mu * (1.0 - mu) + 
                        2.0 * b * sumR * var_mu_ext;
    double denominator = var_int + 4.0 * b * b * sumR * sumR * var_mu_ext;
    if (denominator <= 0.0) {
        Rcpp::Rcout << "    Invalid denominator: " << denominator << std::endl;
        return arma::datum::nan;
    }
    cov_Sbat_S = cov_Sbat_S * std::sqrt(var_S / denominator);
    
    // Step 10: Create VAR matrix and calculate p0
    arma::mat VAR(2, 2);
    VAR(0, 0) = var_S;
    VAR(0, 1) = cov_Sbat_S;
    VAR(1, 0) = cov_Sbat_S;
    VAR(1, 1) = var_Sbat;
    
    // Call mvtnorm::pmvnorm
    Rcpp::Function pmvnorm("pmvnorm", Rcpp::Environment::namespace_env("mvtnorm"));
    
    Rcpp::NumericVector lower0 = Rcpp::NumericVector::create(R_NegInf, lb / std::sqrt(var_ratio_w0));
    Rcpp::NumericVector upper0 = Rcpp::NumericVector::create(-std::abs(S / std::sqrt(var_ratio0)), ub / std::sqrt(var_ratio_w0));
    Rcpp::NumericVector mean0 = Rcpp::NumericVector::create(0.0, 0.0);
    
    Rcpp::NumericVector p0_result = pmvnorm(
        Rcpp::Named("lower", lower0),
        Rcpp::Named("upper", upper0), 
        Rcpp::Named("mean", mean0),
        Rcpp::Named("sigma", VAR)
    );
    double p0 = std::max(0.0, std::min(1.0, p0_result[0]));
    
    // Step 11: sigma2 != 0 case
    arma::vec spa_result_s1 = SPA_G_one_SNP_homo_cpp(
        g, R, mu_ext, n_ext, b, sigma2, var_ratio1, 2.0, 0.15, 10.0);
    double p_spa_s1 = spa_result_s1(0);
    
    double var_S1 = S * S / var_ratio1 / qchisq_boost(p_spa_s1, 1.0, false, false);
    
    double cov_Sbat_S1 = arma::sum(w1 % R_minus_factor) * 2.0 * mu * (1.0 - mu) + 
                         2.0 * b * sumR * (var_mu_ext + sigma2);
    cov_Sbat_S1 = cov_Sbat_S1 * std::sqrt(var_S1 / (var_int + 4.0 * b * b * sumR * sumR * (var_mu_ext + sigma2)));
    
    double var_Sbat1 = var_Sbat + sigma2;
    
    // Step 12: Create VAR1 matrix and calculate p1
    arma::mat VAR1(2, 2);
    VAR1(0, 0) = var_S1;
    VAR1(0, 1) = cov_Sbat_S1;
    VAR1(1, 0) = cov_Sbat_S1;
    VAR1(1, 1) = var_Sbat1;
    
    Rcpp::NumericVector lower1 = Rcpp::NumericVector::create(R_NegInf, lb / std::sqrt(var_ratio_w1));
    Rcpp::NumericVector upper1 = Rcpp::NumericVector::create(-std::abs(S / std::sqrt(var_ratio1)), ub / std::sqrt(var_ratio_w1));
    Rcpp::NumericVector mean1 = Rcpp::NumericVector::create(0.0, 0.0);
    
    Rcpp::NumericVector p1_result = pmvnorm(
        Rcpp::Named("lower", lower1),
        Rcpp::Named("upper", upper1),
        Rcpp::Named("mean", mean1),
        Rcpp::Named("sigma", VAR1)
    );
    
    double p1 = std::max(0.0, std::min(1.0, p1_result[0]));
    
    // Step 13: Final result
    double p_con = 2.0 * (TPR * p1 + (1.0 - TPR) * p0) / p_deno;
    return p_con;
}

} // namespace WtCoxG