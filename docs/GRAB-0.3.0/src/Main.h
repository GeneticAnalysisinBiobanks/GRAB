#include <RcppArmadillo.h>


// Unified interface for extracting genotype data from different file formats.
// Handles missing data tracking and provides comprehensive marker information.
arma::vec Unified_getOneMarker(
  std::string t_genoType,                    // Genotype file format ("PLINK", "BGEN")
  uint64_t t_gIndex,                         // Marker index in genotype file
  std::string &t_ref,                        // Reference allele (output)
  std::string &t_alt,                        // Alternate allele (output)
  std::string &t_marker,                     // Marker ID (output)
  uint32_t &t_pd,                            // Physical position (output)
  std::string &t_chr,                        // Chromosome (output)
  double &t_altFreq,                         // Alternate allele frequency (output)
  double &t_altCounts,                       // Alternate allele count (output)
  double &t_missingRate,                     // Missing genotype rate (output)
  double &t_imputeInfo,                      // Imputation quality score (output)
  bool t_isOutputIndexForMissing,            // Whether to output missing data indices
  std::vector<uint32_t> &t_indexForMissing,  // Indices of missing genotypes (output)
  bool t_isOnlyOutputNonZero,                // Whether to output only non-zero genotypes
  std::vector<uint32_t> &t_indexForNonZero   // Indices of non-zero genotypes (output)
);
