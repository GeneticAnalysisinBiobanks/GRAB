#include <RcppArmadillo.h>

namespace PLINK { class PlinkClass; }
namespace BGEN { class BgenClass; }
namespace DenseGRM { class DenseGRMClass; }
namespace POLMM { class POLMMClass; }
namespace SPACox { class SPACoxClass; }
namespace SPAmix { class SPAmixClass; }
namespace SPAGRM { class SPAGRMClass; }
namespace SAGELD { class SAGELDClass; }
namespace WtCoxG { class WtCoxGClass; }
namespace SPAsqr { class SPAsqrClass; }
namespace LEAF { class LEAFClass; }
namespace SPAmixPlus { class SPAmixPlusClass; }

extern PLINK::PlinkClass* ptr_gPLINKobj;
extern BGEN::BgenClass* ptr_gBGENobj;
extern DenseGRM::DenseGRMClass* ptr_gDenseGRMobj;

extern POLMM::POLMMClass* ptr_gPOLMMobj;
extern SPACox::SPACoxClass* ptr_gSPACoxobj;
extern SPAmix::SPAmixClass* ptr_gSPAmixobj;
extern SPAGRM::SPAGRMClass* ptr_gSPAGRMobj;
extern SAGELD::SAGELDClass* ptr_gSAGELDobj;
extern WtCoxG::WtCoxGClass* ptr_gWtCoxGobj;
extern SPAsqr::SPAsqrClass* ptr_gSPAsqrobj;
extern LEAF::LEAFClass* ptr_gLEAFobj;
extern SPAmixPlus::SPAmixPlusClass* ptr_gSPAmixPlusobj;

extern std::string g_impute_method;
extern double g_missingRate_cutoff;
extern unsigned int g_omp_num_threads;
extern double g_marker_minMAF_cutoff;
extern double g_marker_minMAC_cutoff;


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
