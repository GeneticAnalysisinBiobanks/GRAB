// SPAmixLocalPlus UKB v20 - Ultimate Performance Optimization - zlib Streaming Version
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _WIN32
#pragma comment(lib, "zlib")
#endif

#include "SPAmixLocalPlus.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <chrono>

static constexpr double MIN_P_VALUE = std::numeric_limits<double>::min();

namespace SPAmixLocalPlus {
    arma::vec g_resid;
    std::vector<std::string> g_subjData;
    arma::uvec g_outliers;
    int g_save_interval = 100;
    double g_MAF_cutoff = 0.00001;
    double g_MAC_cutoff = 1.0;
    double g_cutoff = 2.0;
    bool g_verbose = true;
}

// ==================== Data Structure Implementation ====================

PhiData::PhiData(const arma::mat& phi_matrix) {
  if (phi_matrix.n_rows == 0 || phi_matrix.n_cols != 3) {
    return;
  }
  i_idx = arma::conv_to<arma::uvec>::from(phi_matrix.col(0)) - 1;
  j_idx = arma::conv_to<arma::uvec>::from(phi_matrix.col(1)) - 1;
  phi_value = phi_matrix.col(2);
}

// ==================== Helper Functions ====================

bool is_gzipped(const std::string& filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f.good()) return false;
    unsigned char byte1, byte2;
    f >> byte1 >> byte2;
    return (byte1 == 0x1f && byte2 == 0x8b);
}

bool file_exists_cpp(const std::string& filename) {
    std::ifstream f(filename.c_str());
    return f.good();
}

// FileReader Class (Simplified for basic reading)
class FileReader {
private:
    gzFile gz_file;
    std::ifstream txt_file;
    bool is_gz;
    bool is_open_flag;
    
public:
    FileReader(const std::string& filename) : gz_file(nullptr), is_open_flag(false) {
        if (!file_exists_cpp(filename)) return;
        
        is_gz = is_gzipped(filename);
        
        if (is_gz) {
            gz_file = gzopen(filename.c_str(), "rb");
            if (gz_file != nullptr) {
                is_open_flag = true;
            }
        } else {
            txt_file.open(filename);
            if (txt_file.is_open()) {
                is_open_flag = true;
            }
        }
    }
    
    ~FileReader() {
        close();
    }
    
    bool is_open() const {
        return is_open_flag;
    }
    
    bool getline(std::string& line) {
        if (!is_open_flag) return false;
        
        if (is_gz) {
            const int BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
            char* buffer = new char[BUFFER_SIZE];
            
            if (gzgets(gz_file, buffer, BUFFER_SIZE) != nullptr) {
                line = buffer;
                if (!line.empty() && line.back() == '\n') line.pop_back();
                if (!line.empty() && line.back() == '\r') line.pop_back();
                delete[] buffer;
                return true;
            }
            delete[] buffer;
            return false;
        } else {
            if (std::getline(txt_file, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            return false;
        }
    }
    
    void close() {
        if (is_gz && gz_file != nullptr) {
            gzclose(gz_file);
            gz_file = nullptr;
        } else if (!is_gz && txt_file.is_open()) {
            txt_file.close();
        }
        is_open_flag = false;
    }
};

int count_file_lines(const std::string& filename) {
    FileReader reader(filename);
    if (!reader.is_open()) {
        return 0;
    }
    int count = 0;
    std::string line;
    while (reader.getline(line)) {
        count++;
    }
    return count;
}


// ==================== KGF Derivatives ====================

arma::vec K_G0_vec_cpp(const arma::vec& t_vec, double MAF, const arma::vec& h) {
  arma::vec exp_t = arma::exp(t_vec);
  arma::vec base = (1.0 - MAF) + MAF * exp_t;
  arma::vec result(t_vec.n_elem);
  
  for (arma::uword i = 0; i < t_vec.n_elem; ++i) {
    if (base(i) > 0) {
      result(i) = h(i) * std::log(base(i));
    } else {
      result(i) = -arma::datum::inf;
    }
  }
  return result;
}

arma::vec K_G1_vec_cpp(const arma::vec& t_vec, double MAF, const arma::vec& h) {
  arma::vec exp_t = arma::exp(t_vec);
  arma::vec denom = (1.0 - MAF) + MAF * exp_t;
  arma::vec result(t_vec.n_elem);
  
  for (arma::uword i = 0; i < t_vec.n_elem; ++i) {
    if (std::abs(denom(i)) > 1e-15) {
      result(i) = h(i) * MAF * exp_t(i) / denom(i);
    } else {
      result(i) = 0.0;
    }
  }
  return result;
}

arma::vec K_G2_vec_cpp(const arma::vec& t_vec, double MAF, const arma::vec& h) {
  arma::vec exp_t = arma::exp(t_vec);
  arma::vec denom = (1.0 - MAF) + MAF * exp_t;
  arma::vec numer = h % (MAF * exp_t * (1.0 - MAF));
  arma::vec result(t_vec.n_elem);
  
  for (arma::uword i = 0; i < t_vec.n_elem; ++i) {
    double denom_sq = denom(i) * denom(i);
    if (denom_sq > 1e-15) {
      result(i) = numer(i) / denom_sq;
    } else {
      result(i) = 0.0;
    }
  }
  return result;
}

// ==================== Core Calculation Functions ====================

double calculate_var_related_optimized_cpp(
    const arma::vec& R,
    const arma::vec& haplo_num,
    double q,
    const PhiData& phi_A,
    const PhiData& phi_B,
    const PhiData& phi_C,
    const PhiData& phi_D
) {
  double q_term = q * (1.0 - q);
  double var_total = 0.0;
  
  // Scenario A: h_i=2, h_j=2
  if (!phi_A.empty()) {
    for (arma::uword k = 0; k < phi_A.size(); ++k) {
      arma::uword i = phi_A.i_idx(k);
      arma::uword j = phi_A.j_idx(k);
      
      if (i < haplo_num.n_elem && j < haplo_num.n_elem &&
          i < R.n_elem && j < R.n_elem) {
        if (std::abs(haplo_num(i) - 2.0) < 1e-10 && std::abs(haplo_num(j) - 2.0) < 1e-10) {
          var_total += 4.0 * q_term * phi_A.phi_value(k) * R(i) * R(j);
        }
      }
    }
  }
  
  // Scenario B: h_i=2, h_j=1
  if (!phi_B.empty()) {
    for (arma::uword k = 0; k < phi_B.size(); ++k) {
      arma::uword i = phi_B.i_idx(k);
      arma::uword j = phi_B.j_idx(k);
      
      if (i < haplo_num.n_elem && j < haplo_num.n_elem &&
          i < R.n_elem && j < R.n_elem) {
        if (std::abs(haplo_num(i) - 2.0) < 1e-10 && std::abs(haplo_num(j) - 1.0) < 1e-10) {
          var_total += 2.0 * q_term * phi_B.phi_value(k) * R(i) * R(j);
        }
      }
    }
  }
  
  // Scenario C: h_i=1, h_j=2
  if (!phi_C.empty()) {
    for (arma::uword k = 0; k < phi_C.size(); ++k) {
      arma::uword i = phi_C.i_idx(k);
      arma::uword j = phi_C.j_idx(k);
      
      if (i < haplo_num.n_elem && j < haplo_num.n_elem &&
          i < R.n_elem && j < R.n_elem) {
        if (std::abs(haplo_num(i) - 1.0) < 1e-10 && std::abs(haplo_num(j) - 2.0) < 1e-10) {
          var_total += 2.0 * q_term * phi_C.phi_value(k) * R(i) * R(j);
        }
      }
    }
  }
  
  // Scenario D: h_i=1, h_j=1
  if (!phi_D.empty()) {
    for (arma::uword k = 0; k < phi_D.size(); ++k) {
      arma::uword i = phi_D.i_idx(k);
      arma::uword j = phi_D.j_idx(k);
      
      if (i < haplo_num.n_elem && j < haplo_num.n_elem &&
          i < R.n_elem && j < R.n_elem) {
        if (std::abs(haplo_num(i) - 1.0) < 1e-10 && std::abs(haplo_num(j) - 1.0) < 1e-10) {
          var_total += 1.0 * q_term * phi_D.phi_value(k) * R(i) * R(j);
        }
      }
    }
  }
  
  // Add all i=j autocorrelation variance terms
  for (arma::uword i = 0; i < R.n_elem; ++i) {
    var_total += R(i) * R(i) * haplo_num(i) * q_term;
  }
  
  return var_total;
}

NewtonResult get_root_partial_cpp(
    double s,
    const arma::vec& R_outlier,
    const arma::vec& h_outlier,
    double q,
    double mean_normal,
    double var_normal,
    double init_t = 0.0,
    double tol = 0.001,
    int max_iter = 100
) {
  arma::vec init_values = {0.0, -1.0, 1.0, -2.0, 2.0};
  
  for (arma::uword init_idx = 0; init_idx < init_values.n_elem; ++init_idx) {
    double t_cur = init_values(init_idx);
    bool converged = false;
    int iter = 0;
    
    for (iter = 1; iter <= max_iter; ++iter) {
      double K1_out = 0.0, K2_out = 0.0;
      
      if (R_outlier.n_elem > 0) {
        arma::vec tR = t_cur * R_outlier;
        tR = arma::clamp(tR, -700.0, 700.0);
        
        arma::vec K1_terms = K_G1_vec_cpp(tR, q, h_outlier);
        arma::vec K2_terms = K_G2_vec_cpp(tR, q, h_outlier);
        
        K1_out = arma::sum(R_outlier % K1_terms);
        K2_out = arma::sum((R_outlier % R_outlier) % K2_terms);
      }
      
      double K1_total = K1_out + mean_normal + var_normal * t_cur - s;
      double K2_total = K2_out + var_normal;
      
      if (std::abs(K1_total) < tol) {
        converged = true;
        break;
      }
      if (K2_total <= 1e-10) break;
      
      double delta_t = -K1_total / K2_total;
      delta_t = std::max(std::min(delta_t, 2.0), -2.0);
      t_cur = t_cur + delta_t;
      t_cur = std::max(std::min(t_cur, 20.0), -20.0);
      
      if (std::abs(delta_t) < tol) {
        converged = true;
        break;
      }
    }
    
    if (converged) {
      return NewtonResult(t_cur, true, iter);
    }
  }
  
  return NewtonResult(0.0, false, max_iter);
}

SPAResult partial_spa_pval_cpp(
    double s,
    const arma::vec& R,
    const arma::vec& haplo_num,
    double q,
    double var_S,
    const arma::uvec& posOutlier,
    double Cutoff
) {
  double S_mean = arma::sum(q * haplo_num % R);
  double z = (s - S_mean) / std::sqrt(var_S);
  double p_norm = 2.0 * R::pnorm(-std::abs(z), 0.0, 1.0, true, false);
  
  if (std::abs(z) < Cutoff) {
    return SPAResult(p_norm, p_norm, false);
  }
  
  arma::vec MAF_est_Vec(R.n_elem);
  MAF_est_Vec.fill(q);
  arma::vec g_var_est_Vec = haplo_num % MAF_est_Vec % (1.0 - MAF_est_Vec);
  double S_var_SPAmix_local = arma::sum(R % R % g_var_est_Vec);
  double Var_ratio = S_var_SPAmix_local / var_S;
  
  double S_new = s * std::sqrt(Var_ratio);
  double S_mean_new = S_mean * std::sqrt(Var_ratio);
  
  if (posOutlier.n_elem == 0) {
    return SPAResult(p_norm, p_norm, false);
  }
  
  arma::uvec posOutlier_0based = posOutlier - 1;
  arma::uvec all_indices = arma::regspace<arma::uvec>(0, R.n_elem - 1);
  arma::uvec posNormal = all_indices;
  for (arma::uword i = 0; i < posOutlier_0based.n_elem; ++i) {
    arma::uword outlier_idx = posOutlier_0based(i);
    if (outlier_idx < R.n_elem) {
      posNormal = posNormal.elem(arma::find(posNormal != outlier_idx));
    }
  }
  
  arma::vec R_out = R.elem(posOutlier_0based);
  arma::vec h_out = haplo_num.elem(posOutlier_0based);
  arma::vec R_norm = R.elem(posNormal);
  arma::vec h_norm = haplo_num.elem(posNormal);
  
  double mean_normal = arma::sum(R_norm % (q * h_norm));
  double var_normal = arma::sum(R_norm % R_norm % (q * (1.0 - q) * h_norm));
  
  double s_upper = std::max(S_new, 2.0 * S_mean_new - S_new);
  double s_lower = std::min(S_new, 2.0 * S_mean_new - S_new);
  
  // Upper tail probability
  NewtonResult root_info_upper = get_root_partial_cpp(s_upper, R_out, h_out, q, mean_normal, var_normal);
  double zeta_upper = root_info_upper.converge ? root_info_upper.root : 0.0;
  
  double K_out_upper = 0.0, K2_out_upper = 0.0;
  if (R_out.n_elem > 0) {
    arma::vec tR_out_upper = zeta_upper * R_out;
    K_out_upper = arma::sum(K_G0_vec_cpp(tR_out_upper, q, h_out));
    K2_out_upper = arma::sum((R_out % R_out) % K_G2_vec_cpp(tR_out_upper, q, h_out));
  }
  
  double k1_upper = K_out_upper + mean_normal * zeta_upper + 0.5 * var_normal * zeta_upper * zeta_upper;
  double k2_upper = K2_out_upper + var_normal;
  double temp1_upper = zeta_upper * s_upper - k1_upper;
  
  double log_pval1 = 0.0;
  bool use_log_pval1 = false;
  double pval1 = 0.0;
  
  if (temp1_upper <= 0 || k2_upper <= 0) {
    if (std::abs(z) > 8.0) {
      log_pval1 = R::pnorm(std::abs(z), 0.0, 1.0, false, true);
      if (!std::isfinite(log_pval1) || log_pval1 < -700.0) {
        double z_abs = std::abs(z);
        log_pval1 = -0.5 * z_abs * z_abs - std::log(z_abs) - 0.5 * std::log(2.0 * M_PI);
      }
      use_log_pval1 = true;
    } else {
        pval1 = R::pnorm(std::abs(z), 0.0, 1.0, false, false);
        if (pval1 <= 0.0) pval1 = MIN_P_VALUE;
    }
  } else {
    double w_upper = std::copysign(std::sqrt(2.0 * temp1_upper), zeta_upper);
    double v_upper = zeta_upper * std::sqrt(k2_upper);
    
    if (std::abs(w_upper) < 1e-12 || std::abs(v_upper) < 1e-12 || 
        !std::isfinite(w_upper) || !std::isfinite(v_upper)) {
        if (std::abs(z) >= Cutoff) {
          pval1 = MIN_P_VALUE;
        } else {
          if (std::abs(z) > 8.0) {
            log_pval1 = R::pnorm(std::abs(z), 0.0, 1.0, false, true);
            if (!std::isfinite(log_pval1) || log_pval1 < -700.0) {
              double z_abs = std::abs(z);
              log_pval1 = -0.5 * z_abs * z_abs - std::log(z_abs) - 0.5 * std::log(2.0 * M_PI);
            }
            use_log_pval1 = true;
          } else {
            pval1 = R::pnorm(std::abs(z), 0.0, 1.0, false, false);
            if (pval1 <= 0.0) pval1 = MIN_P_VALUE;
          }
        }
    } else {
      double lr_arg_upper = w_upper + (1.0 / w_upper) * std::log(v_upper / w_upper);
      if (!std::isfinite(lr_arg_upper)) {
        if (std::abs(z) >= Cutoff) {
          pval1 = MIN_P_VALUE;
        } else {
          if (std::abs(z) > 8.0) {
            log_pval1 = R::pnorm(std::abs(z), 0.0, 1.0, false, true);
            use_log_pval1 = true;
          } else {
            pval1 = R::pnorm(std::abs(z), 0.0, 1.0, false, false);
            if (pval1 <= 0.0) pval1 = MIN_P_VALUE;
          }
        }
      } else {
        if (std::abs(lr_arg_upper) > 8.0) {
          log_pval1 = R::pnorm(lr_arg_upper, 0.0, 1.0, false, true);
          if (!std::isfinite(log_pval1) || log_pval1 < -700.0) {
            double lr_abs = std::abs(lr_arg_upper);
            log_pval1 = -0.5 * lr_abs * lr_abs - std::log(lr_abs) - 0.5 * std::log(2.0 * M_PI);
          }
          use_log_pval1 = true;
        } else {
          pval1 = R::pnorm(lr_arg_upper, 0.0, 1.0, false, false);
          if (pval1 <= 0.0) pval1 = MIN_P_VALUE;
        }
      }
    }
  }
  
  NewtonResult root_info_lower = get_root_partial_cpp(s_lower, R_out, h_out, q, mean_normal, var_normal);
  double zeta_lower = root_info_lower.converge ? root_info_lower.root : 0.0;
  
  double K_out_lower = 0.0, K2_out_lower = 0.0;
  if (R_out.n_elem > 0) {
    arma::vec tR_out_lower = zeta_lower * R_out;
    K_out_lower = arma::sum(K_G0_vec_cpp(tR_out_lower, q, h_out));
    K2_out_lower = arma::sum((R_out % R_out) % K_G2_vec_cpp(tR_out_lower, q, h_out));
  }
  
  double k1_lower = K_out_lower + mean_normal * zeta_lower + 0.5 * var_normal * zeta_lower * zeta_lower;
  double k2_lower = K2_out_lower + var_normal;
  double temp1_lower = zeta_lower * s_lower - k1_lower;
  
  double log_pval2 = 0.0;
  bool use_log_pval2 = false;
  double pval2 = 0.0;
  
  if (temp1_lower <= 0 || k2_lower <= 0) {
    if (std::abs(z) > 8.0) {
      log_pval2 = R::pnorm(-std::abs(z), 0.0, 1.0, true, true);
      if (!std::isfinite(log_pval2) || log_pval2 < -700.0) {
        double z_abs = std::abs(z);
        log_pval2 = -0.5 * z_abs * z_abs - std::log(z_abs) - 0.5 * std::log(2.0 * M_PI);
      }
      use_log_pval2 = true;
    } else {
        pval2 = R::pnorm(-std::abs(z), 0.0, 1.0, true, false);
        if (pval2 <= 0.0) pval2 = MIN_P_VALUE;
    }
  } else {
    double w_lower = std::copysign(std::sqrt(2.0 * temp1_lower), zeta_lower);
    double v_lower = zeta_lower * std::sqrt(k2_lower);
    
    if (std::abs(w_lower) < 1e-12 || std::abs(v_lower) < 1e-12 || 
        !std::isfinite(w_lower) || !std::isfinite(v_lower)) {
        if (std::abs(z) >= Cutoff) {
          pval2 = MIN_P_VALUE;
        } else {
          if (std::abs(z) > 8.0) {
            log_pval2 = R::pnorm(-std::abs(z), 0.0, 1.0, true, true);
            if (!std::isfinite(log_pval2) || log_pval2 < -700.0) {
              double z_abs = std::abs(z);
              log_pval2 = -0.5 * z_abs * z_abs - std::log(z_abs) - 0.5 * std::log(2.0 * M_PI);
            }
            use_log_pval2 = true;
          } else {
            pval2 = R::pnorm(-std::abs(z), 0.0, 1.0, true, false);
            if (pval2 <= 0.0) pval2 = MIN_P_VALUE;
          }
        }
    } else {
      double lr_arg_lower = w_lower + (1.0 / w_lower) * std::log(v_lower / w_lower);
      if (!std::isfinite(lr_arg_lower)) {
        if (std::abs(z) >= Cutoff) {
          pval2 = MIN_P_VALUE;
        } else {
          if (std::abs(z) > 8.0) {
            log_pval2 = R::pnorm(-std::abs(z), 0.0, 1.0, true, true);
            use_log_pval2 = true;
          } else {
            pval2 = R::pnorm(-std::abs(z), 0.0, 1.0, true, false);
            if (pval2 <= 0.0) pval2 = MIN_P_VALUE;
          }
        }
      } else {
        if (std::abs(lr_arg_lower) > 8.0) {
          log_pval2 = R::pnorm(lr_arg_lower, 0.0, 1.0, true, true);
          if (!std::isfinite(log_pval2) || log_pval2 < -700.0) {
            double lr_abs = std::abs(lr_arg_lower);
            log_pval2 = -0.5 * lr_abs * lr_abs - std::log(lr_abs) - 0.5 * std::log(2.0 * M_PI);
          }
          use_log_pval2 = true;
        } else {
          pval2 = R::pnorm(lr_arg_lower, 0.0, 1.0, true, false);
          if (pval2 <= 0.0) pval2 = MIN_P_VALUE;
        }
      }
    }
  }
  
  double pval_spa_total;
  bool upper_newton_success = root_info_upper.converge;
  bool lower_newton_success = root_info_lower.converge;
  
  if (!upper_newton_success && !lower_newton_success) {
    pval_spa_total = R_NaN;
  } else if (!upper_newton_success && lower_newton_success) {
    if (use_log_pval2) {
      if (log_pval2 > -600.0) {
        pval_spa_total = std::exp(log_pval2);
      } else {
    pval_spa_total = MIN_P_VALUE;
      }
    } else {
      pval_spa_total = pval2;
    }
  } else if (upper_newton_success && !lower_newton_success) {
    if (use_log_pval1) {
      if (log_pval1 > -600.0) {
        pval_spa_total = std::exp(log_pval1);
      } else {
    pval_spa_total = MIN_P_VALUE;
      }
    } else {
      pval_spa_total = pval1;
    }
  } else {
    if (use_log_pval1 && use_log_pval2) {
      double max_log = std::max(log_pval1, log_pval2);
      double log_sum = max_log + std::log(std::exp(log_pval1 - max_log) + std::exp(log_pval2 - max_log));
      
      if (log_sum > -600.0) {  
        pval_spa_total = std::exp(log_sum);
      } else {
        pval_spa_total = MIN_P_VALUE;
      }
    } else if (use_log_pval1) {
      if (log_pval1 > -600.0) {
        pval_spa_total = std::exp(log_pval1) + pval2;
      } else {
    pval_spa_total = pval2 + MIN_P_VALUE;
      }
    } else if (use_log_pval2) {
      if (log_pval2 > -600.0) {
        pval_spa_total = pval1 + std::exp(log_pval2);
      } else {
    pval_spa_total = pval1 + MIN_P_VALUE;
      }
    } else {
      pval_spa_total = pval1 + pval2;
    }
  }
  
  if (ISNAN(pval_spa_total)) {
  } else if (!std::isfinite(pval_spa_total)) {
  } else if (pval_spa_total <= 0) {
    pval_spa_total = MIN_P_VALUE;
  } else if (pval_spa_total > 1.0) {
    pval_spa_total = 1.0;
  }
  
  if (!ISNAN(pval_spa_total) && pval_spa_total == 0.0) {
    pval_spa_total = MIN_P_VALUE;
  }
  if (p_norm == 0.0) {
    p_norm = MIN_P_VALUE;
  }
  
  return SPAResult(pval_spa_total, p_norm, true);
}


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
    double MAF_cutoff,
    double MAC_cutoff
) {
    double total_alleles = arma::sum(haplo_num);
    
    if (total_alleles <= 0.0) {
        arma::vec result(9);
        result(0) = 0.0;  
        result(1) = arma::sum(arma::conv_to<arma::vec>::from(arma::find_nonfinite(g))) / static_cast<double>(g.n_elem);  
        result(2) = arma::datum::nan;  
        result(3) = arma::datum::nan;  
        result(4) = arma::datum::nan;  
        result(5) = arma::datum::nan;  
        result(6) = arma::datum::nan;  
        result(7) = arma::datum::nan;   
        result(8) = 0.0;  
        return result;
    }
    
    double q = arma::sum(g) / total_alleles;
    arma::uword N = g.n_elem;
    
    arma::vec g_processed = g;  
    
    double AltAlleleCount = arma::sum(g_processed);
    
    double minor_allele_count = (q > 0.5) ? (total_alleles - AltAlleleCount) : AltAlleleCount;
    
    if (q < MAF_cutoff || q > (1.0 - MAF_cutoff) || minor_allele_count < MAC_cutoff) {
        arma::vec result(9);
        result(0) = q;  
        result(1) = arma::sum(arma::conv_to<arma::vec>::from(arma::find_nonfinite(g))) / static_cast<double>(N);  
        result(2) = arma::datum::nan;  
        result(3) = arma::datum::nan;  
        result(4) = arma::datum::nan;  
        result(5) = arma::datum::nan;  
        result(6) = arma::datum::nan;  
        result(7) = arma::datum::nan;  
        result(8) = AltAlleleCount;  
        return result;
    }
    
    double S = arma::sum(g_processed % R);
    double S_mean = arma::sum(q * haplo_num % R);
    
    double var_S = calculate_var_related_optimized_cpp(
        R, haplo_num, q, phi_A, phi_B, phi_C, phi_D
    );
    
    arma::vec result(9);  
    
    if (var_S <= 0.0) {
        result(0) = q;  
        result(1) = arma::sum(arma::conv_to<arma::vec>::from(arma::find_nonfinite(g))) / static_cast<double>(N);  
        result(2) = arma::datum::nan;  
        result(3) = arma::datum::nan;  
        result(4) = S;  
        result(5) = S_mean;  
        result(6) = var_S;  
        result(7) = 0.0;  
        result(8) = AltAlleleCount;  
    } else {
        double z = (S - S_mean) / std::sqrt(var_S);
        
        double pval_norm;
        if (std::abs(z) > 8.0) {
            double log_pval_one_tail = R::pnorm(-std::abs(z), 0.0, 1.0, true, true);
            
            if (!std::isfinite(log_pval_one_tail) || log_pval_one_tail < -700.0) {
                double z_abs = std::abs(z);
                log_pval_one_tail = -0.5 * z_abs * z_abs - std::log(z_abs) - 0.5 * std::log(2.0 * M_PI);
            }
            
            double log_pval_two_tail = std::log(2.0) + log_pval_one_tail;
            
            if (log_pval_two_tail > -600.0) {  
                pval_norm = std::exp(log_pval_two_tail);
            } else {
                pval_norm = MIN_P_VALUE;
            }
        } else {
            pval_norm = 2.0 * R::pnorm(-std::abs(z), 0.0, 1.0, true, false);
            if (pval_norm <= 0.0) pval_norm = MIN_P_VALUE;
        }
        
        SPAResult spa_result = partial_spa_pval_cpp(S, R, haplo_num, q, var_S, posOutlier, Cutoff);
        
        double pval_spa_fixed = spa_result.pval_spa;
        
        if (!ISNAN(pval_spa_fixed) && std::isfinite(pval_spa_fixed) && pval_spa_fixed > 1.0) {
            pval_spa_fixed = 1.0;
        }
        
        result(0) = q;  
        result(1) = arma::sum(arma::conv_to<arma::vec>::from(arma::find_nonfinite(g))) / static_cast<double>(N);  
        result(2) = pval_spa_fixed;  
        result(3) = pval_norm;  
        result(4) = S;  
        result(5) = S_mean;  
        result(6) = var_S;  
        result(7) = z;  
        result(8) = AltAlleleCount;  
    }
    
    return result;
}

// ==================== UKBZlibReader Implementation ====================

UKBZlibReader::UKBZlibReader() : geno_gz(nullptr), haplo_gz(nullptr), 
                                 current_line_num(0), initialized(false) {}

UKBZlibReader::~UKBZlibReader() {
    close();
}

bool UKBZlibReader::initialize(const std::string& geno_file, 
                              const std::string& haplo_file,
                              const std::vector<std::string>& target_sample_ids) {
    close();  
    
    this->target_sample_ids = target_sample_ids;
    
    geno_gz = gzopen(geno_file.c_str(), "rb");
    if (geno_gz == nullptr) {
        Rcpp::Rcout << "    Unable to open geno file: " << geno_file << std::endl;
        return false;
    }
    
    haplo_gz = gzopen(haplo_file.c_str(), "rb");
    if (haplo_gz == nullptr) {
        Rcpp::Rcout << "    Unable to open haplo file: " << haplo_file << std::endl;
        gzclose(geno_gz);
        geno_gz = nullptr;
        return false;
    }
    
    gzclose(geno_gz);
    geno_gz = nullptr;
    
    FileReader reader(geno_file);
    if (!reader.is_open()) {
        Rcpp::Rcout << "    Unable to open file to read header: " << geno_file << std::endl;
        return false;
    }
    
    std::string header_line;
    if (!reader.getline(header_line)) {
        Rcpp::Rcout << "    Unable to read file header: " << geno_file << std::endl;
        return false;
    }
    reader.close();
    
    geno_gz = gzopen(geno_file.c_str(), "rb");
    if (geno_gz == nullptr) {
        Rcpp::Rcout << "    Unable to re-open geno file: " << geno_file << std::endl;
        return false;
    }
    
    char skip_buffer[1024];
    if (gzgets(geno_gz, skip_buffer, sizeof(skip_buffer)) == nullptr) {
        Rcpp::Rcout << "    Unable to skip geno file header line" << std::endl;
        return false;
    }
    size_t skipped_len = std::strlen(skip_buffer);
    if (skipped_len == 0 || skip_buffer[skipped_len - 1] != '\n') {
        int ch;
        while ((ch = gzgetc(geno_gz)) != -1 && ch != '\n') {
        }
    }
    
    std::istringstream iss(header_line);
    std::string field;
    std::vector<std::string> header_fields;
    while (std::getline(iss, field, '\t')) {
        header_fields.push_back(field);
    }
    
    global_sample_to_col.clear();
    for (size_t i = 5; i < header_fields.size(); ++i) {  
        global_sample_to_col[header_fields[i]] = static_cast<int>(i);  
    }
    
    if (SPAmixLocalPlus::g_verbose) {
        Rcpp::Rcout << "    Sample mapping construction complete: " << global_sample_to_col.size() 
                   << " UKB samples, processing " << target_sample_ids.size() << " target samples" << std::endl;
    }
    
    char skip_haplo_buffer[1024];
    if (gzgets(haplo_gz, skip_haplo_buffer, sizeof(skip_haplo_buffer)) == nullptr) {
        Rcpp::Rcout << "    Unable to skip haplo file header line" << std::endl;
        return false;
    }
    size_t skipped_haplo_len = std::strlen(skip_haplo_buffer);
    if (skipped_haplo_len == 0 || skip_haplo_buffer[skipped_haplo_len - 1] != '\n') {
        int ch;
        while ((ch = gzgetc(haplo_gz)) != -1 && ch != '\n') {
        }
    }
    
    current_line_num = 1;  
    initialized = true;
    
    if (SPAmixLocalPlus::g_verbose) {
        Rcpp::Rcout << "    zlib reader initialized successfully" << std::endl;
    }
    return true;
}

bool UKBZlibReader::read_next_snp(std::vector<double>& geno_data, 
                                  std::vector<double>& haplo_data,
                                  std::string& snp_id,
                                  std::array<std::string, 5>& marker_fields) {
    if (!initialized || geno_gz == nullptr || haplo_gz == nullptr) {
        return false;
    }
    
    const int buffer_size = 10 * 1024 * 1024;  
    char* buffer = new char[buffer_size];
    
    try {
        if (gzgets(geno_gz, buffer, buffer_size) == nullptr) {
            delete[] buffer;
            return false;  
        }
        current_geno_line = std::string(buffer);
        if (current_geno_line.back() == '\n') current_geno_line.pop_back();
        
        if (gzgets(haplo_gz, buffer, buffer_size) == nullptr) {
            delete[] buffer;
            Rcpp::Rcout << "    Error geno and haplo file line count mismatch" << std::endl;
            return false;
        }
        current_haplo_line = std::string(buffer);
        if (current_haplo_line.back() == '\n') current_haplo_line.pop_back();
        
        delete[] buffer;
        
        std::istringstream geno_iss(current_geno_line);
        std::string field;
        std::vector<std::string> geno_fields;
        while (std::getline(geno_iss, field, '\t')) {
            geno_fields.push_back(field);
        }
        
        std::istringstream haplo_iss(current_haplo_line);
        std::vector<std::string> haplo_fields;
        while (std::getline(haplo_iss, field, '\t')) {
            haplo_fields.push_back(field);
        }
        
        marker_fields = {"", "", "", "", ""};
        for (int meta_idx = 0; meta_idx < 5 && meta_idx < static_cast<int>(geno_fields.size()); ++meta_idx) {
            marker_fields[meta_idx] = geno_fields[meta_idx];
        }

        if (geno_fields.size() > 2) {
            snp_id = geno_fields[2];
        } else {
            snp_id = "SNP_" + std::to_string(current_line_num);
        }
        
        geno_data.clear();
        haplo_data.clear();
        geno_data.reserve(target_sample_ids.size());
        haplo_data.reserve(target_sample_ids.size());
        
        for (size_t i = 0; i < target_sample_ids.size(); ++i) {
            const std::string& sample_id = target_sample_ids[i];
            auto it = global_sample_to_col.find(sample_id);
            if (it != global_sample_to_col.end()) {
                int col_idx = it->second;
                if (col_idx < static_cast<int>(geno_fields.size()) && 
                    col_idx < static_cast<int>(haplo_fields.size())) {
                    
                    const std::string& geno_field_val = geno_fields[col_idx];
                    const std::string& haplo_field_val = haplo_fields[col_idx];
                    
                    if (geno_field_val.empty() || geno_field_val == "NA" || geno_field_val == "." || geno_field_val == "-9" ||
                        haplo_field_val.empty() || haplo_field_val == "NA" || haplo_field_val == "." || haplo_field_val == "-9") {
                        geno_data.push_back(R_NaReal);
                        haplo_data.push_back(R_NaReal);
                    } else {
                        try {
                            double geno_val = std::stod(geno_field_val);
                            double haplo_val = std::stod(haplo_field_val);
                            geno_data.push_back(geno_val);
                            haplo_data.push_back(haplo_val);
                        } catch (...) {
                            geno_data.push_back(R_NaReal);
                            haplo_data.push_back(R_NaReal);
                        }
                    }
                } else {
                    geno_data.push_back(R_NaReal);
                    haplo_data.push_back(R_NaReal);
                }
            } else {
                geno_data.push_back(R_NaReal);
                haplo_data.push_back(R_NaReal);
            }
        }
        
        current_line_num++;
        return true;
        
    } catch (...) {
        delete[] buffer;
        Rcpp::Rcout << "    Error Exception occurred while reading SNP data" << std::endl;
        return false;
    }
}

void UKBZlibReader::close() {
    if (geno_gz != nullptr) {
        gzclose(geno_gz);
        geno_gz = nullptr;
    }
    if (haplo_gz != nullptr) {
        gzclose(haplo_gz);
        haplo_gz = nullptr;
    }
    initialized = false;
}

UKBZlibReader* create_ukb_zlib_reader(const std::string& geno_file,
                                     const std::string& haplo_file,
                                     const std::vector<std::string>& target_sample_ids) {
    UKBZlibReader* reader = new UKBZlibReader();
    if (reader->initialize(geno_file, haplo_file, target_sample_ids)) {
        return reader;
    } else {
        delete reader;
        return nullptr;
    }
}

// ==================== Exported Functions ====================

// [[Rcpp::export]]
void SPAmixPlusLocal_streamInCPP(
    const std::string& geno_file,
    const std::string& haplo_file,
    const std::string& output_file,
    const arma::uvec& file_match_idx,
    const arma::mat& phi_A_mat,
    const arma::mat& phi_B_mat,
    const arma::mat& phi_C_mat,
    const arma::mat& phi_D_mat
) {
    using namespace SPAmixLocalPlus;
    
    auto overall_start = std::chrono::high_resolution_clock::now();
    
    arma::vec R_matched(file_match_idx.n_elem);
    for (arma::uword i = 0; i < file_match_idx.n_elem; i++) {
        R_matched(i) = g_resid(file_match_idx(i));
    }
    
    std::vector<arma::uword> outlier_indices;
    for (arma::uword i = 0; i < g_outliers.n_elem; i++) {
        arma::uword global_idx = g_outliers(i);
        for (arma::uword j = 0; j < file_match_idx.n_elem; j++) {
            if (file_match_idx(j) == global_idx) {
                outlier_indices.push_back(j);
                break;
            }
        }
    }
    arma::uvec posOutlier = arma::conv_to<arma::uvec>::from(outlier_indices);
    
    if (g_verbose) {
        Rcpp::Rcout << "\n    SPAmixPlusLocal Streaming Analysis" << std::endl;
        Rcpp::Rcout << "    Genotype File: " << geno_file << std::endl;
        Rcpp::Rcout << "    Haplotype File: " << haplo_file << std::endl;
        Rcpp::Rcout << "    Output File: " << output_file << std::endl;
        Rcpp::Rcout << "    Matched Samples: " << R_matched.n_elem << std::endl;
        Rcpp::Rcout << "    Outliers (remapped): " << posOutlier.n_elem << std::endl;
        Rcpp::Rcout << "    Save Interval: Every " << g_save_interval << " SNPs" << std::endl;
        Rcpp::Rcout << "    Filter Params: MAF >= " << g_MAF_cutoff << ", MAC >= " << g_MAC_cutoff << std::endl;
    }
    
    int total_snps = count_file_lines(geno_file);
    if (total_snps > 0) {
        total_snps = total_snps - 1;  
    }
    if (total_snps <= 0) total_snps = 10000000;  
    
    std::vector<std::string> target_sample_ids;
    target_sample_ids.reserve(file_match_idx.n_elem);
    
    std::vector<std::string> ukb_sample_ids = read_ukb_sample_ids_cpp(geno_file);
    
    if (g_verbose) {
        Rcpp::Rcout << "    Stats UKB File Total Samples: " << ukb_sample_ids.size() << std::endl;
        Rcpp::Rcout << "    Target file_match_idx Range: [" << file_match_idx.min() << ", " << file_match_idx.max() << "]" << std::endl;
    }
    
    for (arma::uword i = 0; i < file_match_idx.n_elem; ++i) {
        arma::uword ukb_idx = file_match_idx(i);  
        target_sample_ids.push_back(ukb_sample_ids[ukb_idx]);
        if (g_verbose && i < 5) {
            Rcpp::Rcout << "      Mapping[" << i << "]: file_match_idx=" << ukb_idx 
                       << " -> sample_id=" << ukb_sample_ids[ukb_idx] << std::endl;
        }
    }
    
    auto init_start = std::chrono::high_resolution_clock::now();
    UKBZlibReader* reader = create_ukb_zlib_reader(geno_file, haplo_file, target_sample_ids);
    if (reader == nullptr) {
        Rcpp::stop("Error Unable to initialize zlib reader");
    }
    
    // Test read
    std::vector<double> test_geno, test_haplo;
    std::string test_snp_id;
    std::array<std::string, 5> test_marker_fields;
    arma::uword expected_sample_count = file_match_idx.n_elem;
    
    if (reader->read_next_snp(test_geno, test_haplo, test_snp_id, test_marker_fields)) {
        if (test_geno.size() == expected_sample_count && 
            test_haplo.size() == expected_sample_count &&
            test_geno.size() == R_matched.n_elem) {
            if (g_verbose) {
                Rcpp::Rcout << "    Success Data format validated: " << expected_sample_count 
                           << " samples, Test SNP: " << test_snp_id << std::endl;
            }
        } else {
            delete reader;
            Rcpp::stop("Error Data format mismatch: Expected " + std::to_string(expected_sample_count) + 
                      " samples, Actual geno=" + std::to_string(test_geno.size()) + 
                      ", haplo=" + std::to_string(test_haplo.size()));
        }
    } else {
        delete reader;
        Rcpp::stop("Error Unable to read test SNP for format validation");
    }
    
    delete reader;
    reader = create_ukb_zlib_reader(geno_file, haplo_file, target_sample_ids);
    if (reader == nullptr) {
        Rcpp::stop("Error Unable to re-initialize zlib reader");
    }
    
    auto init_end = std::chrono::high_resolution_clock::now();
    auto init_duration = std::chrono::duration_cast<std::chrono::milliseconds>(init_end - init_start);
    
    if (g_verbose) {
        Rcpp::Rcout << "    zlib reader initialized, data format validated (" << init_duration.count() << "ms)" << std::endl;
    }
    
    std::ofstream output(output_file);
    if (!output.is_open()) {
        delete reader;
        Rcpp::stop("Error Unable to create output file: " + output_file);
    }
    
    output << "rsID\tCHROM\tPOS\tREF\tALT\tAltFreq\tMissingRate\tPvalue\tPvalueNormal\tStat\tStatMean\tStatVar\tzScore\tAltCounts\tBetaG\n";

    std::vector<std::string> result_buffer;
    result_buffer.reserve(g_save_interval);
    
    int processed_snps = 0;
    int saved_snps = 0;
    std::vector<double> geno_data, haplo_data;
    std::string snp_id;
    std::array<std::string, 5> marker_fields;
    
    auto process_start = std::chrono::high_resolution_clock::now();
    
    while (processed_snps < total_snps && reader->read_next_snp(geno_data, haplo_data, snp_id, marker_fields)) {
        
        arma::vec g(geno_data.data(), geno_data.size(), false, true);
        arma::vec haplo_num(haplo_data.data(), haplo_data.size(), false, true);
        
        try {
            arma::vec result = SPAmixPlus_local_related_one_SNP_cpp(
                g, R_matched, haplo_num,
                PhiData(phi_A_mat), PhiData(phi_B_mat), 
                PhiData(phi_C_mat), PhiData(phi_D_mat),
                posOutlier, g_cutoff, g_MAF_cutoff, g_MAC_cutoff
            );
            
            if (result.n_elem >= 9) {
                double maf = result(0);          
                double missing_rate = result(1); 
                double pval_spa = result(2);     
                double pval_norm = result(3);    
                double stat = result(4);         
                double mean = result(5);         
                double var = result(6);          
                double z = result(7);            
                double mac = result(8);          
                
                double beta_g = std::numeric_limits<double>::quiet_NaN();
                if (std::isfinite(stat) && std::isfinite(var) && var != 0.0) {
                    beta_g = stat / var;
                }

                auto extract_meta = [&](int idx) -> std::string {
                    if (idx < 0 || idx >= static_cast<int>(marker_fields.size())) {
                        return "NA";
                    }
                    const std::string& value = marker_fields[idx];
                    if (value.empty()) {
                        return "NA";
                    }
                    return value;
                };

                std::string chrom = extract_meta(0);
                std::string pos = extract_meta(1);
                std::string ref = extract_meta(3);
                std::string alt = extract_meta(4);

                std::ostringstream result_line;
                result_line << snp_id << "\t"
                           << chrom << "\t"
                           << pos << "\t"
                           << ref << "\t"
                           << alt << "\t"
                           << std::fixed << std::setprecision(8) << maf << "\t"
                           << std::fixed << std::setprecision(0) << missing_rate << "\t"
                           << std::scientific << std::setprecision(8) << pval_spa << "\t"
                           << std::scientific << std::setprecision(8) << pval_norm << "\t"
                           << std::fixed << std::setprecision(8) << stat << "\t"
                           << std::fixed << std::setprecision(8) << mean << "\t"
                           << std::fixed << std::setprecision(8) << var << "\t"
                           << std::fixed << std::setprecision(8) << z << "\t"
                           << std::fixed << std::setprecision(0) << mac << "\t";

                if (std::isfinite(beta_g)) {
                    result_line << std::fixed << std::setprecision(8) << beta_g;
                } else {
                    result_line << "NA";
                }
                
                result_buffer.push_back(result_line.str());
            } else {
                std::ostringstream result_line;
                result_line << snp_id
                           << "\tNA\tNA\tNA\tNA"   
                           << "\tNA\tNA\tNA\tNA"   
                           << "\tNA\tNA\tNA\tNA\tNA\tNA"; 
                result_buffer.push_back(result_line.str());
            }
            
        } catch (...) {
            std::ostringstream result_line;
            result_line << snp_id
                       << "\tNA\tNA\tNA\tNA"
                       << "\tNA\tNA\tNA\tNA"
                       << "\tNA\tNA\tNA\tNA\tNA\tNA";
            result_buffer.push_back(result_line.str());
        }
        
        processed_snps++;
        
        if (result_buffer.size() >= static_cast<size_t>(g_save_interval) || processed_snps >= total_snps) {
            for (const std::string& line : result_buffer) {
                output << line << "\n";
            }
            output.flush();  
            
            saved_snps += result_buffer.size();
            result_buffer.clear();
            
            if (g_verbose) {
                auto current_time = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - process_start);
                Rcpp::Rcout << "    Processed " << processed_snps << "/" << total_snps 
                           << " SNPs, Saved " << saved_snps << " results ("
                           << elapsed.count() << "s)" << std::endl;
            }
        }
    }
    
    output.close();
    delete reader;
    
    auto overall_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(overall_end - overall_start);
    
    if (g_verbose) {
        Rcpp::Rcout << "\n    SPAmixPlusLocal Analysis complete" << std::endl;
        Rcpp::Rcout << "    Total Processed: " << processed_snps << " SNPs" << std::endl;
        Rcpp::Rcout << "    Results saved to: " << output_file << std::endl;
        Rcpp::Rcout << "    Total time: " << total_duration.count() << " seconds" << std::endl;
    }
    
}

// [[Rcpp::export]]
std::vector<std::string> read_ukb_sample_ids_cpp(const std::string& geno_file) {
    std::vector<std::string> sample_ids;
    
    FileReader reader(geno_file);
    if (!reader.is_open()) {
        Rcpp::stop("Unable to read file: " + geno_file);
    }
    
    std::string header_line;
    if (!reader.getline(header_line)) {
        Rcpp::stop("Unable to read file header: " + geno_file);
    }
    
    std::istringstream iss(header_line);
    std::string token;
    int col_idx = 0;
    
    while (std::getline(iss, token, '\t')) {
        if (col_idx >= 5) {  
            sample_ids.push_back(token);
        }
        col_idx++;
    }
    
    if (SPAmixLocalPlus::g_verbose) {
        Rcpp::Rcout << "    Extracted " << sample_ids.size() << " sample IDs from UKB file header" << std::endl;
    }
    
    return sample_ids;
}

// [[Rcpp::export]]
Rcpp::List SPAmixLocalPlus_computePhiInCPP(
    const arma::mat& hapcount_matrix,
    const arma::mat& dosage_matrix,
    const arma::uvec& pair_idx1,
    const arma::uvec& pair_idx2,
    const std::string& scenario,
    double phi_threshold = 0.0,
    double maf_cutoff = 0.01
) {
    int n_snps = hapcount_matrix.n_rows;
    int n_samples = hapcount_matrix.n_cols;
    int n_pairs = pair_idx1.n_elem;
    
    if ((int)pair_idx2.n_elem != n_pairs) {
        Rcpp::stop("Pair index vectors must have same length");
    }
    
    std::vector<arma::uword> directed_pair_idx1, directed_pair_idx2;
    for (int k = 0; k < n_pairs; k++) {
        if (pair_idx1(k) != pair_idx2(k)) {
            directed_pair_idx1.push_back(pair_idx1(k));
            directed_pair_idx2.push_back(pair_idx2(k));
            directed_pair_idx1.push_back(pair_idx2(k));
            directed_pair_idx2.push_back(pair_idx1(k));
        }
    }
    
    arma::uvec final_pair_idx1 = arma::conv_to<arma::uvec>::from(directed_pair_idx1);
    arma::uvec final_pair_idx2 = arma::conv_to<arma::uvec>::from(directed_pair_idx2);
    int n_directed_pairs = final_pair_idx1.n_elem;
    
    if ((int)dosage_matrix.n_rows != n_snps || (int)dosage_matrix.n_cols != n_samples) {
        Rcpp::stop("Hapcount and dosage matrix dimensions mismatch");
    }
    
    arma::vec ratio_sums(n_directed_pairs, arma::fill::zeros);
    arma::uvec valid_counts(n_directed_pairs, arma::fill::zeros);
    
    int snps_passed_maf = 0;
    double maf_high_threshold = 1.0 - maf_cutoff;
    
    for (int s = 0; s < n_snps; s++) {
        double total_dosage = 0.0;
        double total_hapcount = 0.0;
        
        for (int i = 0; i < n_samples; i++) {
            total_dosage += dosage_matrix(s, i);
            total_hapcount += hapcount_matrix(s, i);
        }
        
        double global_maf = (total_hapcount > 1e-10) ? (total_dosage / total_hapcount) : 0.0;
        
        if (global_maf <= maf_cutoff || global_maf >= maf_high_threshold) {
            continue;
        }
        
        snps_passed_maf++;
        
        for (int p = 0; p < n_directed_pairs; p++) {
            int idx1 = final_pair_idx1(p);
            int idx2 = final_pair_idx2(p);
            
            if (idx1 >= (int)n_samples || idx2 >= (int)n_samples) continue;
            
            double h_i = hapcount_matrix(s, idx1);
            double h_j = hapcount_matrix(s, idx2);
            double g_i = dosage_matrix(s, idx1);
            double g_j = dosage_matrix(s, idx2);
            
            double q = global_maf;
            
            bool scenario_match = false;
            int hi_int = (int)std::round(h_i);
            int hj_int = (int)std::round(h_j);
            
            if (scenario == "A" && hi_int == 2 && hj_int == 2) scenario_match = true;
            else if (scenario == "B" && hi_int == 2 && hj_int == 1) scenario_match = true;
            else if (scenario == "C" && hi_int == 1 && hj_int == 2) scenario_match = true;
            else if (scenario == "D" && hi_int == 1 && hj_int == 1) scenario_match = true;
            
            if (scenario_match) {
                double numerator = (g_i - h_i * q) * (g_j - h_j * q);
                double denominator = h_i * h_j * q * (1.0 - q);
                
                if (std::abs(denominator) > 1e-15) {
                    ratio_sums(p) += numerator / denominator;
                    valid_counts(p) += 1;
                }
            }
        }
    }
    
    return Rcpp::List::create(
        Rcpp::Named("ratio_sums") = ratio_sums,
        Rcpp::Named("valid_counts") = valid_counts,
        Rcpp::Named("directed_pair_idx1") = final_pair_idx1,
        Rcpp::Named("directed_pair_idx2") = final_pair_idx2,
        Rcpp::Named("n_snps_processed") = snps_passed_maf,
        Rcpp::Named("n_pairs") = n_directed_pairs
    );
}

// [[Rcpp::export]]
void SPAmixPlusLocal_setupInCPP(
    const arma::vec& resid,
    const std::vector<std::string>& subjData,
    Rcpp::Nullable<Rcpp::List> outLierList,
    int save_interval,
    double MAF_cutoff,
    double MAC_cutoff,
    double cutoff,
    bool verbose
) {
    using namespace SPAmixLocalPlus;
    
    g_resid = resid;
    g_subjData = subjData;
    g_save_interval = save_interval;
    g_MAF_cutoff = MAF_cutoff;
    g_MAC_cutoff = MAC_cutoff;
    g_cutoff = cutoff;
    g_verbose = verbose;
    
    if (outLierList.isNotNull()) {
        Rcpp::List outList(outLierList);
        if (outList.containsElementNamed("posOutlier")) {
            Rcpp::IntegerVector outlier_vec = outList["posOutlier"];
            g_outliers = arma::conv_to<arma::uvec>::from(Rcpp::as<arma::vec>(outlier_vec));
        } else {
            g_outliers = arma::uvec();  
        }
    } else {
        g_outliers = arma::uvec();  
    }
    
    if (g_verbose) {
        Rcpp::Rcout << "    SPAmixPlusLocal Global Setup" << std::endl;
        Rcpp::Rcout << "    Residuals: " << g_resid.n_elem << " samples" << std::endl;
        Rcpp::Rcout << "    Outliers: " << g_outliers.n_elem << std::endl;
        Rcpp::Rcout << "    MAF cutoff: " << g_MAF_cutoff << std::endl;
        Rcpp::Rcout << "    MAC cutoff: " << g_MAC_cutoff << std::endl;
        Rcpp::Rcout << "    SPA cutoff: " << g_cutoff << std::endl;
    }
}

// [[Rcpp::export]]
List getSampleMatchIndices_cpp(const std::vector<std::string>& file_sample_ids) {
    using namespace SPAmixLocalPlus;
    
    std::unordered_map<std::string, int> global_map;
    for (size_t i = 0; i < g_subjData.size(); i++) {
        global_map[g_subjData[i]] = i;
    }
    
    std::vector<int> file_indices;  
    std::vector<std::string> matched_sample_ids;  
    std::vector<int> global_to_matched;  
    global_to_matched.resize(g_subjData.size(), -1);
    
    int matched_idx = 0;
    for (size_t file_idx = 0; file_idx < file_sample_ids.size(); file_idx++) {
        auto it = global_map.find(file_sample_ids[file_idx]);
        if (it != global_map.end()) {
            int global_idx = it->second;
            file_indices.push_back(file_idx);
            matched_sample_ids.push_back(file_sample_ids[file_idx]);
            global_to_matched[global_idx] = matched_idx;
            matched_idx++;
        }
    }
    
    std::vector<int> valid_outliers;
    for (size_t i = 0; i < g_outliers.n_elem; i++) {
        int global_outlier_idx = g_outliers(i);
        if (global_outlier_idx >= 0 && global_outlier_idx < (int)global_to_matched.size()) {
            int matched_outlier_idx = global_to_matched[global_outlier_idx];
            if (matched_outlier_idx >= 0) {
                valid_outliers.push_back(matched_outlier_idx);
            }
        }
    }
    
    if (g_verbose) {
        Rcpp::Rcout << "    Sample matching: " << matched_sample_ids.size() 
                    << "/" << g_subjData.size() << " global samples found in file" << std::endl;
        Rcpp::Rcout << "    Remapped outliers: " << valid_outliers.size() 
                    << "/" << g_outliers.n_elem << std::endl;
    }
    
    return List::create(
        Named("file_indices") = file_indices,
        Named("global_to_matched") = global_to_matched, 
        Named("valid_outliers") = valid_outliers
    );
}

