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
extern arma::uvec g_group;
extern bool g_ifOutGroup;
extern unsigned int g_nGroup;
extern arma::sp_mat g_SparseGRM;

void setPOLMMobjInCPP(
  arma::mat t_muMat,
  arma::mat t_iRMat,
  arma::mat t_Cova,
  arma::uvec t_yVec,
  double t_tau,
  bool t_printPCGInfo,
  double t_tolPCG,
  int t_maxiterPCG,
  double t_varRatio,
  double t_SPA_cutoff,
  bool t_flagSparseGRM,
  arma::uvec t_group,
  bool t_ifOutGroup,
  unsigned int t_nGroup
);

void setSPAGRMobjInCPP(
  arma::vec t_resid,
  arma::vec t_resid_unrelated_outliers,
  double t_sum_R_nonOutlier,
  double t_R_GRM_R_nonOutlier,
  double t_R_GRM_R_TwoSubjOutlier,
  double t_R_GRM_R,
  arma::vec t_MAF_interval,
  Rcpp::List t_TwoSubj_list,
  Rcpp::List t_ThreeSubj_list,
  double t_SPA_Cutoff,
  double t_zeta,
  double t_tol
);

void setSAGELDobjInCPP(
  std::string t_Method,
  arma::mat t_XTs,
  arma::mat t_SS,
  arma::mat t_AtS,
  arma::mat t_Q,
  arma::mat t_A21,
  arma::mat t_TTs,
  arma::mat t_Tys,
  arma::vec t_sol,
  arma::vec t_blups,
  double t_sig,
  arma::vec t_resid,
  arma::vec t_resid_G,
  arma::vec t_resid_GxE,
  arma::vec t_resid_E,
  arma::vec t_resid_unrelated_outliers,
  arma::vec t_resid_unrelated_outliers_G,
  arma::vec t_resid_unrelated_outliers_GxE,
  double t_sum_R_nonOutlier,
  double t_sum_R_nonOutlier_G,
  double t_sum_R_nonOutlier_GxE,
  double t_R_GRM_R,
  double t_R_GRM_R_G,
  double t_R_GRM_R_GxE,
  double t_R_GRM_R_G_GxE,
  double t_R_GRM_R_E,
  double t_R_GRM_R_nonOutlier,
  double t_R_GRM_R_nonOutlier_G,
  double t_R_GRM_R_nonOutlier_GxE,
  double t_R_GRM_R_nonOutlier_G_GxE,
  double t_R_GRM_R_TwoSubjOutlier,
  double t_R_GRM_R_TwoSubjOutlier_G,
  double t_R_GRM_R_TwoSubjOutlier_GxE,
  double t_R_GRM_R_TwoSubjOutlier_G_GxE,
  Rcpp::List t_TwoSubj_list,
  Rcpp::List t_ThreeSubj_list,
  arma::vec t_MAF_interval,
  double t_zScoreE_cutoff,
  double t_SPA_Cutoff,
  double t_zeta,
  double t_tol
);

void setSPAmixobjInCPP(
  arma::mat t_resid,
  arma::mat t_PCs,
  int t_N,
  double t_SPA_Cutoff,
  Rcpp::List t_outlierList
);

void setSPACoxobjInCPP(
  arma::mat t_cumul,
  arma::vec t_mresid,
  arma::mat t_XinvXX,
  arma::mat t_tX,
  int t_N,
  double t_pVal_covaAdj_Cutoff,
  double t_SPA_Cutoff
);

void setWtCoxGobjInCPP(
  arma::vec t_mresid,
  arma::vec t_weight,
  double t_cutoff,
  double t_SPA_Cutoff
);

void setSPAsqrobjInCPP(
  arma::vec t_taus,
  arma::mat t_Resid_mat,
  Rcpp::List t_Resid_unrelated_outliers_lst,
  arma::vec t_sum_R_nonOutlier_vec,
  arma::vec t_R_GRM_R_nonOutlier_vec,
  arma::vec t_R_GRM_R_TwoSubjOutlier_vec,
  arma::vec t_R_GRM_R_vec,
  arma::vec t_MAF_interval,
  Rcpp::List t_TwoSubj_list_lst,
  Rcpp::List t_CLT_union_lst,
  Rcpp::List t_ThreeSubj_family_idx_lst,
  Rcpp::List t_ThreeSubj_stand_S_lst,
  double t_SPA_Cutoff,
  double t_zeta,
  double t_tol
);

void setLEAFobjInCPP(
  Rcpp::List t_residuals,
  Rcpp::List t_weight,
  Rcpp::List t_clusterIdx,
  double t_cutoff,
  double t_SPA_Cutoff
);

void setSPAmixPlusobjInCPP(
  arma::mat t_resid,
  arma::mat t_PCs,
  int t_N,
  double t_SPA_Cutoff,
  Rcpp::List t_outlierList,
  Rcpp::DataFrame t_sparseGRM,
  std::string t_afFilePath,
  std::string t_afFilePrecision
);


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
