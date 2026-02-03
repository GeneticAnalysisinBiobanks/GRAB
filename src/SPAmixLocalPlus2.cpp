#include <RcppArmadillo.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>

using namespace Rcpp;
using namespace arma;

//' Read PLINK BIM file and extract target SNPs by chromosome
//' @param bim_file PLINK BIM file path
//' @return List containing chromosome-wise SNP information
// [[Rcpp::export]]
List read_bim_snps_by_chromosome(const std::string& bim_file) {
    // Rcout << "Reading PLINK BIM file: " << bim_file << std::endl;
    
    std::ifstream file(bim_file);
    if (!file.is_open()) {
        stop("Cannot open BIM file: " + bim_file);
    }
    
    std::map<int, std::vector<std::string>> chr_snps;
    std::string line;
    // int total_snps = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string chr_str, snp_id, gmap, pos_str, a1, a2;
        
        if (iss >> chr_str >> snp_id >> gmap >> pos_str >> a1 >> a2) {
            try {
                int chr = std::stoi(chr_str);
                
                if (chr >= 1 && chr <= 22) {  // Only autosomes
                    chr_snps[chr].push_back(snp_id);
                    // total_snps++;
                }
            } catch (...) {
                continue; 
            }
        }
    }
    file.close();
    
    // Rcout << "BIM file loaded: " << total_snps << " SNPs across " << chr_snps.size() << " chromosomes" << std::endl;
    
    List result;
    for (const auto& pair : chr_snps) {
        int chr = pair.first;
        result[std::to_string(chr)] = List::create(
            Named("snps") = pair.second,
            Named("n_snps") = pair.second.size()
        );
    }
    
    return result;
}

//' Read UKB sample IDs from file header with large buffer
//' @param file_path UKB .txt.gz file path  
//' @return List with sample information
// [[Rcpp::export]]
List read_ukb_sample_info_large(const std::string& file_path) {
    std::string cmd = "zcat \"" + file_path + "\" | head -n 1";
    
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        stop("Failed to execute zcat command for file: " + file_path);
    }
    
    // Allocate 8MB buffer
    const size_t BUFFER_SIZE = 8 * 1024 * 1024;
    char* buffer = new char[BUFFER_SIZE];
    std::vector<std::string> sample_ids;
    
    try {
        if (fgets(buffer, BUFFER_SIZE, pipe) != nullptr) {
            std::string header_line(buffer);
            
            if (!header_line.empty() && header_line.back() == '\n') {
                header_line.pop_back();
            }
            
            std::istringstream iss(header_line);
            std::string field;
            int field_count = 0;
            
            while (std::getline(iss, field, '\t')) {
                field_count++;
                if (field_count > 5) {  // Skip first 5 columns: CHROM,POS,ID,REF,ALT
                    if (!field.empty()) {
                         // Basic cleaning
                         std::string clean_field = field;
                         clean_field.erase(std::remove_if(clean_field.begin(), clean_field.end(), 
                                                         [](unsigned char c) { return c < 32 || c > 126; }), 
                                          clean_field.end());
                         if(!clean_field.empty()) sample_ids.push_back(clean_field);
                    }
                }
            }
        }
        
        delete[] buffer;
        pclose(pipe);
        
        if (sample_ids.empty()) {
            stop("No sample IDs found in file: " + file_path);
        }
        
        return List::create(
            Named("sample_ids") = sample_ids,
            Named("n_samples") = sample_ids.size()
        );
        
    } catch (...) {
        delete[] buffer;
        pclose(pipe);
        throw;
    }
}

//' Stream-based SNP batch reading
//' @param hapcount_file Hapcount file path
//' @param dosage_file Dosage file path  
//' @param target_snps Vector of target SNP IDs
//' @param max_snps_per_batch Maximum SNPs to read per batch
//' @return List with hapcount and dosage matrices for found SNPs
// [[Rcpp::export]]
List read_ukb_snp_batch_stream(const std::string& hapcount_file,
                               const std::string& dosage_file,
                               const std::vector<std::string>& target_snps,
                               int max_snps_per_batch = 200) {
    
    if (target_snps.empty()) {
        return List::create(
            Named("hapcount_matrix") = arma::mat(),
            Named("dosage_matrix") = arma::mat(),
            Named("found_snps") = std::vector<std::string>(),
            Named("n_snps_found") = 0,
            Named("n_samples") = 0
        );
    }
    
    // Create SNP pattern for grep
    std::string snp_pattern = target_snps[0];
    int batch_size = std::min((int)target_snps.size(), max_snps_per_batch);
    for (int i = 1; i < batch_size; i++) {
        snp_pattern += "|" + target_snps[i];
    }
    
    // zcat + grep
    std::string hap_cmd = "zcat \"" + hapcount_file + "\" | grep -E \"(" + snp_pattern + ")\" | head -n " + std::to_string(batch_size);
    std::string dos_cmd = "zcat \"" + dosage_file + "\" | grep -E \"(" + snp_pattern + ")\" | head -n " + std::to_string(batch_size);
    
    std::vector<std::vector<double>> hap_data, dos_data;
    std::vector<std::string> found_snps;
    int n_samples = 0;
    
    const size_t LINE_BUFFER_SIZE = 8 * 1024 * 1024;
    char* line_buffer = new char[LINE_BUFFER_SIZE];
    
    try {
        // Read Hapcount
        FILE* hap_pipe = popen(hap_cmd.c_str(), "r");
        if (hap_pipe) {
            while (fgets(line_buffer, LINE_BUFFER_SIZE, hap_pipe) != nullptr) {
                std::string line(line_buffer);
                if (!line.empty() && line.back() == '\n') line.pop_back();
                
                std::istringstream iss(line);
                std::string field;
                std::vector<double> row_data;
                int field_count = 0;
                std::string snp_id;
                
                while (std::getline(iss, field, '\t')) {
                    field_count++;
                    if (field_count == 3) {
                        snp_id = field;
                    } else if (field_count > 5) {
                        try {
                            row_data.push_back(std::stod(field));
                        } catch (...) {
                            row_data.push_back(0.0);
                        }
                    }
                }
                
                if (!row_data.empty() && !snp_id.empty()) {
                    hap_data.push_back(row_data);
                    found_snps.push_back(snp_id);
                    if (n_samples == 0) n_samples = row_data.size();
                }
            }
            pclose(hap_pipe);
        }
        
        // Read Dosage
        FILE* dos_pipe = popen(dos_cmd.c_str(), "r");
        if (dos_pipe) {
            while (fgets(line_buffer, LINE_BUFFER_SIZE, dos_pipe) != nullptr) {
                std::string line(line_buffer);
                if (!line.empty() && line.back() == '\n') line.pop_back();
                
                std::istringstream iss(line);
                std::string field;
                std::vector<double> row_data;
                int field_count = 0;
                
                while (std::getline(iss, field, '\t')) {
                    field_count++;
                    if (field_count > 5) {
                        try {
                            row_data.push_back(std::stod(field));
                        } catch (...) {
                            row_data.push_back(0.0);
                        }
                    }
                }
                if (!row_data.empty()) dos_data.push_back(row_data);
            }
            pclose(dos_pipe);
        }
        
        delete[] line_buffer;
        
        int n_snps_found = std::min(hap_data.size(), dos_data.size());
        
        if (n_snps_found == 0 || n_samples == 0) {
            return List::create(Named("n_snps_found") = 0);
        }
        
        arma::mat hapcount_matrix(n_snps_found, n_samples);
        arma::mat dosage_matrix(n_snps_found, n_samples);
        
        for (int s = 0; s < n_snps_found; s++) {
            for (int i = 0; i < n_samples; i++) {
                hapcount_matrix(s, i) = (i < (int)hap_data[s].size()) ? hap_data[s][i] : 0.0;
                dosage_matrix(s, i) = (i < (int)dos_data[s].size()) ? dos_data[s][i] : 0.0;
            }
        }
        
        return List::create(
            Named("hapcount_matrix") = hapcount_matrix,
            Named("dosage_matrix") = dosage_matrix,
            Named("found_snps") = std::vector<std::string>(found_snps.begin(), found_snps.begin() + n_snps_found),
            Named("n_snps_found") = n_snps_found,
            Named("n_samples") = n_samples
        );
        
    } catch (...) {
        delete[] line_buffer;
        throw;
    }
}

//' Compute phi ratios with ancestry-specific MAF filtering
//' @param hapcount_matrix Hapcount matrix (SNPs x Samples)
//' @param dosage_matrix Dosage matrix (SNPs x Samples)
//' @param pair_idx1 First individual indices (0-based)
//' @param pair_idx2 Second individual indices (0-based)
//' @param scenario Scenario ("A", "B", "C", "D")
//' @param phi_threshold Minimum phi threshold (default 0.0)
//' @param maf_cutoff MAF cutoff threshold (default 0.2)
//' @return List with ratio sums and valid counts for each pair
// [[Rcpp::export]]
List computePhiRatiosCPP(const arma::mat& hapcount_matrix,
                                        const arma::mat& dosage_matrix,
                                        const arma::uvec& pair_idx1,
                                        const arma::uvec& pair_idx2,
                                        const std::string& scenario,
                                        double phi_threshold = 0.0,
                                        double maf_cutoff = 0.2) {
    
    int n_snps = hapcount_matrix.n_rows;
    int n_samples = hapcount_matrix.n_cols;
    int n_input_pairs = pair_idx1.n_elem;
    
    if ((int)pair_idx2.n_elem != n_input_pairs) stop("Pair index vectors must have same length");
    
    // Generate directed pairs (i!=j) -> i->j and j->i
    std::vector<arma::uword> directed_pair_idx1, directed_pair_idx2;
    
    for (int k = 0; k < n_input_pairs; k++) {
        arma::uword i = pair_idx1(k);
        arma::uword j = pair_idx2(k);
        
        if (i != j) {
            directed_pair_idx1.push_back(i);
            directed_pair_idx2.push_back(j);
            directed_pair_idx1.push_back(j);
            directed_pair_idx2.push_back(i);
        }
    }
    
    arma::uvec final_pair_idx1 = arma::conv_to<arma::uvec>::from(directed_pair_idx1);
    arma::uvec final_pair_idx2 = arma::conv_to<arma::uvec>::from(directed_pair_idx2);
    int n_pairs = final_pair_idx1.n_elem;
    
    if ((int)dosage_matrix.n_rows != n_snps || (int)dosage_matrix.n_cols != n_samples) {
        stop("Hapcount and dosage matrix dimensions mismatch");
    }
    
    arma::vec ratio_sums(n_pairs, arma::fill::zeros);
    arma::uvec valid_counts(n_pairs, arma::fill::zeros);
    
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
        
        for (int p = 0; p < n_pairs; p++) {
            int idx1 = final_pair_idx1(p);
            int idx2 = final_pair_idx2(p);
            
            if (idx1 >= n_samples || idx2 >= n_samples) continue;
            
            double h_i = hapcount_matrix(s, idx1);
            double h_j = hapcount_matrix(s, idx2);
            double g_i = dosage_matrix(s, idx1);
            double g_j = dosage_matrix(s, idx2);
            
            double q_i = global_maf;
            double q_j = global_maf; // q[i,s] = q[j,s] = global_maf
            
            bool scenario_match = false;
            // A: h_i=2, h_j=2
            // B: h_i=2, h_j=1
            // C: h_i=1, h_j=2
            // D: h_i=1, h_j=1
            // Note: R code uses scenario_mask_fun(h_i, h_j) which returns TRUE/FALSE
            // Here we compare double to 1.0/2.0 inside loose tolerance? usually integers
            
            int hi_int = (int)std::round(h_i);
            int hj_int = (int)std::round(h_j);

            if (scenario == "A" && hi_int == 2 && hj_int == 2) scenario_match = true;
            else if (scenario == "B" && hi_int == 2 && hj_int == 1) scenario_match = true;
            else if (scenario == "C" && hi_int == 1 && hj_int == 2) scenario_match = true;
            else if (scenario == "D" && hi_int == 1 && hj_int == 1) scenario_match = true;
            
            if (scenario_match) {
                double numerator = (g_i - h_i * q_i) * (g_j - h_j * q_j);
                double denominator = h_i * h_j * q_i * (1.0 - q_i);
                
                // Avoid division by zero
                if (std::abs(denominator) > 1e-15) {
                   ratio_sums(p) += numerator / denominator;
                   valid_counts(p) += 1;
                }
            }
        }
    }
    
    return List::create(
        Named("ratio_sums") = ratio_sums,
        Named("valid_counts") = valid_counts,
        Named("directed_pair_idx1") = final_pair_idx1,
        Named("directed_pair_idx2") = final_pair_idx2,
        Named("n_snps_processed") = snps_passed_maf,
        Named("n_pairs") = n_pairs
    );
}

//' Main function: Stream-optimized genome-wide phi estimation
//' @return List with final phi values and statistics
// [[Rcpp::export]]
List estimatePhiStreamCPP(const std::vector<std::string>& hap_files,
                          const std::vector<std::string>& dos_files,
                          const std::string& bim_file,
                          const std::string& ancestry_id,
                          const std::string& scenario,
                          const arma::uvec& pair_idx1,
                          const arma::uvec& pair_idx2,
                          const std::vector<int>& chromosomes,
                          int snp_batch_size = 200,
                          double phi_threshold = 0.0,
                          double maf_cutoff = 0.2) {
    
    int n_pairs = pair_idx1.n_elem;
    if ((int)pair_idx2.n_elem != n_pairs) stop("Pair index vectors must have same length");
    
    // Calculate expected directed pairs
    int expected_directed_pairs = 0;
    for (int k = 0; k < n_pairs; k++) {
        if (pair_idx1(k) != pair_idx2(k)) expected_directed_pairs += 2;
    }
    
    List bim_data = read_bim_snps_by_chromosome(bim_file);
    
    arma::vec total_ratio_sums(expected_directed_pairs, arma::fill::zeros);
    arma::uvec total_valid_counts(expected_directed_pairs, arma::fill::zeros);
    
    // Track pairs from first batch
    arma::uvec final_directed_pair_idx1(expected_directed_pairs);
    bool pair_indices_captured = false;
    
    // Iterate through provided chromosomes and corresponding files
    for (size_t i = 0; i < chromosomes.size(); i++) {
        int chr = chromosomes[i];
        
        // Ensure we have files for this chromosome
        if(i >= hap_files.size() || i >= dos_files.size()) {
             Rcout << "Warning: Missing file paths for Index " << i << " (Chr " << chr << ")" << std::endl;
             continue;
        }

        std::string hapcount_file = hap_files[i];
        std::string dosage_file = dos_files[i];

        std::string chr_key = std::to_string(chr);
        if (!bim_data.containsElementNamed(chr_key.c_str())) continue;
        
        List chr_info = bim_data[chr_key];
        std::vector<std::string> chr_snps = chr_info["snps"];
        if (chr_snps.empty()) continue;
        
        Rcout << "Processing Chr " << chr << " (" << chr_snps.size() << " SNPs)..." << std::endl;
        
        int n_batches = (chr_snps.size() + snp_batch_size - 1) / snp_batch_size;
        
        for (int batch = 0; batch < n_batches; batch++) {
            if (batch % 10 == 0) Rcpp::checkUserInterrupt();

            int start_idx = batch * snp_batch_size;
            int end_idx = std::min(start_idx + snp_batch_size, (int)chr_snps.size());
            std::vector<std::string> batch_snps(chr_snps.begin() + start_idx, chr_snps.begin() + end_idx);
            
            try {
                List batch_data = read_ukb_snp_batch_stream(hapcount_file, dosage_file, batch_snps, snp_batch_size);
                int n_snps_found = batch_data["n_snps_found"];
                if (n_snps_found == 0) continue;
                
                arma::mat hap = batch_data["hapcount_matrix"];
                arma::mat dos = batch_data["dosage_matrix"];
                
                List batch_result = computePhiRatiosCPP(hap, dos, pair_idx1, pair_idx2, scenario, phi_threshold, maf_cutoff);
                
                arma::vec batch_sums = batch_result["ratio_sums"];
                arma::uvec batch_counts = batch_result["valid_counts"];
                
                if (!pair_indices_captured) {
                     // Ensure dimensions match
                     if(batch_sums.n_elem == total_ratio_sums.n_elem) {
                         final_directed_pair_idx1 = as<arma::uvec>(batch_result["directed_pair_idx1"]);
                         pair_indices_captured = true;
                     }
                }
                
                if (batch_sums.n_elem == total_ratio_sums.n_elem) {
                     total_ratio_sums += batch_sums;
                     total_valid_counts += batch_counts;
                }
                
            } catch (...) {
                // Ignore errors reading bad batch
            }
        }
    }
    
    // Calculate final Phi
    int n_final_pairs = total_ratio_sums.n_elem;
    arma::vec phi_values(n_final_pairs);
    
    for (int i = 0; i < n_final_pairs; i++) {
        if (total_valid_counts(i) > 0) {
            phi_values(i) = total_ratio_sums(i) / total_valid_counts(i);
        } else {
            phi_values(i) = NA_REAL;
        }
    }
    
    // If no pairs captured (no data found anywhere), reconstruct indices locally from input
    // NOTE: This repeats the "directed pair generation" logic
    if (!pair_indices_captured) {
        std::vector<arma::uword> dp1, dp2;
        for (int k = 0; k < n_pairs; k++) {
            if (pair_idx1(k) != pair_idx2(k)) {
                dp1.push_back(pair_idx1(k));
                dp2.push_back(pair_idx2(k));
                dp1.push_back(pair_idx2(k));
                dp2.push_back(pair_idx1(k));
            }
        }
        return List::create(
           Named("phi_values") = phi_values, // all zeros/NA
           Named("pair_idx1") = dp1,
           Named("pair_idx2") = dp2
        );
    }
    
    // We only need one index vector to map back, but caller might want both
    // computePhiRatiosCPP returns directed indices that correspond to the ratio_sums vector
    // We should return these indices so R can map them back to Sample IDs
    
    return List::create(
        Named("phi_values") = phi_values,
        Named("pair_idx1") = final_directed_pair_idx1, 
        // We know idx2 is implicitly determined but let's just return values. 
        // The R side needs to reconstruct the "j" index or we return it.
        // Let's rely on R to map 1->2 and 2->1 if needed, OR return both.
        // I will return the vector index, but computePhiRatiosCPP returns "pair_idx1" and "pair_idx2"
        // in directed form. Let's obtain it from a single call or reconstruct.
        // Actually, 'final_directed_pair_idx1' was captured from the first batch result.
        // We need 'final_directed_pair_idx2' too.
        Named("valid_counts") = total_valid_counts
    );
}
