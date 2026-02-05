#ifndef SPAMIXLOCALPLUSHEADER_HPP
#define SPAMIXLOCALPLUSHEADER_HPP

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <unordered_map>
#include <array>
#include <cstdio>
#include <cstring>  // Modification Date: 2025-09-25 - Added strlen support for header skip logic
#include <functional>
#include <zlib.h>  // Modification Date: 2025-09-11 - Added zlib support for direct .gz reading

using namespace Rcpp;
using namespace arma;

// ==================== Global State Variables ====================
// These are set once via SPAmixPlusLocal_setup and used across all ancestries

namespace SPAmixLocalPlus {
    extern arma::vec g_resid;
    extern std::vector<std::string> g_subjData;
    extern arma::uvec g_outliers;
    extern int g_save_interval;
    extern double g_MAF_cutoff;
    extern double g_MAC_cutoff;
    extern double g_cutoff;
    extern bool g_verbose;
}

// ====================Data Structure Definitions ====================
// Modification Date: 2025-09-03
// Description: Matches data structures from v11 exactly, adapted only for UKB format

struct PhiData {
    arma::uvec i_idx;
    arma::uvec j_idx;
    arma::vec phi_value;
    
    PhiData() {}
    PhiData(const arma::mat& phi_matrix);
    
    // Modification Date: 2025-09-03 - Added SEXP constructor for RCPP support
    PhiData(SEXP sexp_phi) {
        // First try to process as a matrix
        try {
            arma::mat phi_matrix = Rcpp::as<arma::mat>(sexp_phi);
            if (phi_matrix.n_cols >= 3 && phi_matrix.n_rows > 0) {
                i_idx = arma::conv_to<arma::uvec>::from(phi_matrix.col(0));
                j_idx = arma::conv_to<arma::uvec>::from(phi_matrix.col(1));
                phi_value = phi_matrix.col(2);
                return;
            }
        } catch (...) {
            // If matrix conversion fails, try as List
        }
        
        // Process as List
        try {
            Rcpp::List phi_list(sexp_phi);
            if (phi_list.size() >= 3) {
                i_idx = Rcpp::as<arma::uvec>(phi_list["i_idx"]);
                j_idx = Rcpp::as<arma::uvec>(phi_list["j_idx"]);
                phi_value = Rcpp::as<arma::vec>(phi_list["phi_value"]);
            }
        } catch (...) {
            // If all fails, create empty PhiData
        }
    }
    
    bool empty() const { return i_idx.n_elem == 0; }
    size_t size() const { return i_idx.n_elem; }
};

struct NewtonResult {
    double root;
    bool converge;
    int iter;
    
    NewtonResult() : root(R_NaReal), converge(false), iter(0) {}
    NewtonResult(double r, bool c, int i) : root(r), converge(c), iter(i) {}
    
    // Modification Date: 2025-09-03 - Added SEXP conversion support
    operator SEXP() const {
        return Rcpp::List::create(
            Rcpp::Named("root") = root,
            Rcpp::Named("converge") = converge,
            Rcpp::Named("iter") = iter
        );
    }
};

struct SPAResult {
    double pval_spa;
    double pval_normal;
    bool converged;
    
    SPAResult() : pval_spa(R_NaReal), pval_normal(R_NaReal), converged(false) {}
    SPAResult(double spa, double norm, bool conv) : pval_spa(spa), pval_normal(norm), converged(conv) {}
    
    // Modification Date: 2025-09-03 - Added SEXP conversion support
    operator SEXP() const {
        return Rcpp::List::create(
            Rcpp::Named("pval_spa") = pval_spa,
            Rcpp::Named("pval_normal") = pval_normal,
            Rcpp::Named("converged") = converged
        );
    }
};

// ==================== Core Algorithm Function Declarations ====================
// Modification Date: 2025-09-03
// Description: Maintained exact algorithm functions from v11, no calculation logic changes

// CGF Functions (Exact match with v11)
arma::vec K_G0_vec_cpp(const arma::vec& t_vec, double MAF, const arma::vec& h);
arma::vec K_G1_vec_cpp(const arma::vec& t_vec, double MAF, const arma::vec& h);  
arma::vec K_G2_vec_cpp(const arma::vec& t_vec, double MAF, const arma::vec& h);

// Variance Calculation Function (Exact match with v11)
double SPAmixPlus_local_Var_cpp(const arma::vec& g, const arma::vec& haplo_num,
                                const PhiData& phi_A, const PhiData& phi_B,
                                const PhiData& phi_C, const PhiData& phi_D);

// Core Calculation Functions (Exact match with v11)
double Korg_cpp(double t, const arma::vec& mu, const arma::vec& g, double cutoff = 2);
double K1_adj_cpp(double t, const arma::vec& mu, const arma::vec& g, double cutoff = 2);
double K2_adj_cpp(double t, const arma::vec& mu, const arma::vec& g, double cutoff = 2);

// Newton Method Solver (Exact match with v11)
NewtonResult getroot_K1_cpp(double init, const arma::vec& mu, const arma::vec& g, 
                           double cutoff = 2, int maxiter = 100, double tol = 1e-8);

// SPA Method (Exact match with v11)
SPAResult SPA_test_cpp(const arma::vec& g, const arma::vec& mu, const arma::vec& g_tilde,
                      double cutoff = 2, double alpha = 0.05);

// Single SNP Processing Function (Exact match with v11)
// Modification Date: 2025-12-27
// Description: Added support for MAF_cutoff and MAC_cutoff parameters
arma::vec SPAmixPlus_local_related_one_SNP_cpp(
    const arma::vec& g,
    const arma::vec& R,
    const arma::vec& haplo_num,
    const PhiData& phi_A,
    const PhiData& phi_B,
    const PhiData& phi_C,
    const PhiData& phi_D,
    const arma::uvec& posOutlier,
    double Cutoff,
    double MAF_cutoff = 0.0,
    double MAC_cutoff = 0.0
);

// v11 Core Algorithm Functions
double calculate_var_related_optimized_cpp(
    const arma::vec& R,
    const arma::vec& haplo_num,
    double q,
    const PhiData& phi_A,
    const PhiData& phi_B,
    const PhiData& phi_C,
    const PhiData& phi_D
);

NewtonResult get_root_partial_cpp(
    double s,
    const arma::vec& R_outlier,
    const arma::vec& h_outlier,
    double q,
    double mean_normal,
    double var_normal,
    double init_t = 0.0,
    double tol = 1e-6,
    int max_iter = 100
);

SPAResult partial_spa_pval_cpp(
    double s,
    const arma::vec& R,
    const arma::vec& haplo_num,
    double q,
    double var_S,
    const arma::uvec& posOutlier,
    double Cutoff
);

// Debug version of SPA function
SPAResult debug_partial_spa_pval_cpp(
    double s,
    const arma::vec& R,
    const arma::vec& haplo_num,
    double q,
    double var_S,
    const arma::uvec& posOutlier,
    double Cutoff,
    bool verbose = true
);

// ==================== UKB Format Adapter Function Declarations ====================
// Modification Date: 2025-09-03
// Description: Added UKB format processing functions, handling .txt.gz files and transposed matrix format

// UKB Sample ID Reader
std::vector<std::string> read_ukb_sample_ids_cpp(const std::string& geno_file);

// UKB Data Reader Function - Fixed: Reads both genotype and haplotype data
List read_ukb_batch_data_cpp(const std::string& geno_file,
                            const std::string& haplo_file,
                            const std::vector<std::string>& target_sample_ids,
                            int start_snp = 0,
                            int batch_size = 500);

// Modification Date: 2025-09-11 - New: zlib line-by-line reader, no temp files, performance prioritized
class UKBZlibReader {
private:
    gzFile geno_gz, haplo_gz;
    std::vector<std::string> target_sample_ids;  // Target sample ID list
    std::unordered_map<std::string, int> global_sample_to_col;  // Sample ID to column index mapping
    std::string current_geno_line, current_haplo_line;
    int current_line_num;
    bool initialized;
    
public:
    UKBZlibReader();
    ~UKBZlibReader();
    
    bool initialize(const std::string& geno_file, 
                   const std::string& haplo_file,
                   const std::vector<std::string>& target_sample_ids);
    
    bool read_next_snp(std::vector<double>& geno_data, 
                      std::vector<double>& haplo_data,
                      std::string& snp_id,
                      std::array<std::string, 5>& marker_fields);
    
    void close();
    int get_current_line() const { return current_line_num; }
};

UKBZlibReader* create_ukb_zlib_reader(const std::string& geno_file,
                                     const std::string& haplo_file,
                                     const std::vector<std::string>& target_sample_ids);

// ==================================================================
// Phi Estimation
// ==================================================================

// Compute phi ratios for ancestry-specific kinship estimation
Rcpp::List SPAmixLocalPlus_computePhiInCPP(
    const arma::mat& hapcount_matrix,
    const arma::mat& dosage_matrix,
    const arma::uvec& pair_idx1,
    const arma::uvec& pair_idx2,
    const std::string& scenario,
    double phi_threshold,
    double maf_cutoff
);

// ==================================================================
// Global Setup and Helper Functions
// ==================================================================

// Setup global state for SPAmixLocalPlus analysis
void SPAmixPlusLocal_setupInCPP(
    const arma::vec& resid,
    const std::vector<std::string>& subjData,
    Rcpp::Nullable<Rcpp::List> outLierList,
    int save_interval,
    double MAF_cutoff,
    double MAC_cutoff,
    double cutoff,
    bool verbose
);

// Match file samples to global samples and remap outliers
List getSampleMatchIndices_cpp(const std::vector<std::string>& file_sample_ids);

// ==================================================================
// Main Streaming Function
// ==================================================================

// Renamed from SPAmixPlus_local_ukb_high_performance_streaming_cpp
// Now uses global state instead of passing all parameters
void SPAmixPlusLocal_streamInCPP(
    const std::string& geno_file,
    const std::string& haplo_file,
    const std::string& output_file,
    const arma::uvec& file_match_idx,
    const arma::mat& phi_A_mat,
    const arma::mat& phi_B_mat,
    const arma::mat& phi_C_mat,
    const arma::mat& phi_D_mat
);

// ==================== Helper Function Declarations ====================
// Modification Date: 2025-09-03  
// Description: Helper functions maintained consistent with v11

// File Handling
bool file_exists_cpp(const std::string& filename);
bool is_gzipped(const std::string& filename);


#endif // SPAMIXLOCALPLUSHEADER_HPP
