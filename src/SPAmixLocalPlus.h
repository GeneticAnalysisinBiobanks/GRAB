#ifndef SPAMIXLOCALPLUSHEADER_HPP
#define SPAMIXLOCALPLUSHEADER_HPP

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <array>
#include <cmath>
#include <zlib.h> 

using namespace Rcpp;
using namespace arma;

// ==================== Global State Variables ====================
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

// ==================== Data Structures ====================

struct PhiData {
    arma::uvec i_idx;
    arma::uvec j_idx;
    arma::vec phi_value;
    
    PhiData() {}
    PhiData(const arma::mat& phi_matrix);
    
    bool empty() const { return i_idx.n_elem == 0; }
    size_t size() const { return i_idx.n_elem; }
};

struct NewtonResult {
    double root;
    bool converge;
    int iter;
    
    NewtonResult() : root(R_NaReal), converge(false), iter(0) {}
    NewtonResult(double r, bool c, int i) : root(r), converge(c), iter(i) {}
};

struct SPAResult {
    double pval_spa;
    double pval_normal;
    bool converged;
    
    SPAResult() : pval_spa(R_NaReal), pval_normal(R_NaReal), converged(false) {}
    SPAResult(double spa, double norm, bool conv) : pval_spa(spa), pval_normal(norm), converged(conv) {}
};

class UKBZlibReader {
private:
    gzFile geno_gz, haplo_gz;
    std::vector<std::string> target_sample_ids;
    std::unordered_map<std::string, int> global_sample_to_col;
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
};

UKBZlibReader* create_ukb_zlib_reader(const std::string& geno_file,
                                     const std::string& haplo_file,
                                     const std::vector<std::string>& target_sample_ids);

std::vector<std::string> read_ukb_sample_ids_cpp(const std::string& geno_file);

#endif
