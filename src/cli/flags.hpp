// flags.hpp — Flag and method definitions for composable CLI help
//
// Each FlagDef describes one command-line flag.  File-accepting flags carry
// a fileInfo string with detailed format documentation.
//
// Each MethodDef references null-terminated arrays of FlagDef pointers
// (required / optional).  Help text is generated from these at runtime,
// so no per-method help strings need to be maintained by hand.
#pragma once

namespace cli {

// ── Flag definition ────────────────────────────────────────────────

struct FlagDef {
    const char* flag;       // "--threads", or combined display string for groups
    const char* metavar;    // "INT"/"FLOAT"/"FILE"/"PREFIX", nullptr for bool/group
    const char* brief;      // one-line description
    const char* fileInfo;   // extended file-format help, or nullptr
};

// ── Method / utility-mode definition ───────────────────────────────

struct MethodDef {
    const char* name;                   // "SPACox", "cal-ind-af-coef", ...
    const char* desc;                   // one-line description
    const FlagDef* const* required;     // null-terminated
    const FlagDef* const* optional;     // null-terminated
    const char* residNote;              // per-method --null-resid columns, or nullptr
    const char* outputCols;             // output column description
    const char* extra;                  // additional notes, or nullptr
};

// ════════════════════════════════════════════════════════════════════
//  Flag instances
// ════════════════════════════════════════════════════════════════════

// ── Shared required ────────────────────────────────────────────────

inline const FlagDef kBfile = {
    "--bfile", "PREFIX", "PLINK binary genotype prefix (.bed/.bim/.fam)", nullptr};

inline const FlagDef kPfile = {
    "--pfile", "PREFIX", "PLINK2 genotype prefix (.pgen/.pvar/.psam)", nullptr};

inline const FlagDef kVcf = {
    "--vcf", "FILE", "VCF or BCF genotype file", nullptr};

inline const FlagDef kBgen = {
    "--bgen", "FILE", "BGEN genotype file", nullptr};

// Combined display entry for genotype input (exactly one required)
inline const FlagDef kGeno_input = {
    "--bfile PREFIX | --pfile PREFIX | --vcf FILE | --bgen FILE", nullptr,
    "Genotype input (exactly one)", nullptr};

inline const FlagDef kNullResid = {
    "--null-resid", "FILE", "Null model residual file",
    "Formats: (A) header with IID/#IID column, or (B) pure numeric matrix.\n"
    "  Format A: FID/IID detected from header; remaining columns are values.\n"
    "            Header names (beyond IID/FID) do not affect parsing;\n"
    "            column order determines which values are residuals/weights/etc.\n"
    "  Format B: row count must equal .fam subject count; .fam order assumed.\n"
    "\n"
    "Column layout by method:\n"
    "  SPACox / SPAGRM / SPAmix / SPAmixPlus : RESID [RESID2 ...]\n"
    "  WtCoxG / LEAF                          : RESID  WEIGHT  INDICATOR\n"
    "  SPAsqr                                 : R_tau1  R_tau2  ...  R_tauK\n"
    "\n"
    "Multi-column: each value column runs a separate GWAS.\n"
    "  Single column  -> output to --out directly.\n"
    "  Multiple cols  -> PREFIX.1.gz, PREFIX.2.gz, ..."};

inline const FlagDef kOut = {
    "--out", "FILE", "Output file path", nullptr};

inline const FlagDef kCovar = {
    "--covar", "FILE", "Covariate file (plink2 format, columns selected by --covar-name)",
    "Plink2-compatible IID + named columns.\n"
    "Use --covar-name to select specific columns as covariates.\n"
    "If --covar-name is omitted, all non-metadata columns are used."};

inline const FlagDef kPheno = {
    "--pheno", "FILE", "Phenotype file (plink2 format, IID + named columns)",
    "Plink2-compatible IID + named columns.\n"
    "Use --pheno-binary / --pheno-surv / --pheno-quant to designate the response.\n"
    "Use --covar-name to select covariate columns.\n"
    "Use --pc-cols to select PC columns for LEAF / SPAmix."};

inline const FlagDef kCovarName = {
    "--covar-name", "COL_IDS", "Comma-separated covariate column names",
    "Selects columns from --covar (preferred) or --pheno as covariates.\n"
    "An intercept is added automatically."};

inline const FlagDef kBinaryPheno = {
    "--pheno-binary", "COL", "Column name for binary phenotype (case/control)",
    "Selects a 0/1 column from --pheno for logistic regression."};

inline const FlagDef kSurvPheno = {
    "--pheno-surv", "TIME:EVENT", "Colon-separated survival time and event column names",
    "Selects survival time and event columns from --pheno for Cox regression.\n"
    "Example: --pheno-surv TIME:EVENT"};

inline const FlagDef kPcCols = {
    "--pc-cols", "COL_IDS", "Comma-separated PC column names (default: PC1,PC2,PC3,PC4)",
    "Selects columns from --covar or --pheno as principal components.\n"
    "Used for K-means clustering in LEAF and AF estimation in SPAmix."};

inline const FlagDef kRefAf = {
    "--ref-af", "FILE", "Reference allele frequency (plink2 .afreq or numeric)",
    "Formats: (A) plink2 .afreq: #CHROM  ID  REF  ALT  ALT_FREQS  OBS_CT\n"
    "           Matched to .bim by (CHROM, ID) with allele flip detection.\n"
    "         (B) Two-column numeric: ALT_FREQS  OBS_CT (rows in .bim order).\n"
    "For LEAF: comma-separated, one file per reference population."};

inline const FlagDef kSpGrmGrab = {
    "--sp-grm-grab", "FILE", "Sparse GRM, GRAB format (ID1 ID2 VALUE)",
    "Whitespace-delimited; #-lines skipped.\n"
    "Subjects not in .fam are silently dropped.\n"
    "Mutually exclusive with --sp-grm-plink2."};

inline const FlagDef kSpGrmPlink2 = {
    "--sp-grm-plink2", "FILE", "Sparse GRM, plink2 .grm.sp format",
    "Companion .grm.id (FID  IID) auto-detected from file path.\n"
    "If .grm.id is absent, 0-based indices assumed to match .fam order.\n"
    "Mutually exclusive with --sp-grm-grab."};

// Combined display entry used in method required/optional arrays
inline const FlagDef kSpGrm = {
    "--sp-grm-grab FILE | --sp-grm-plink2 FILE", nullptr,
    "Sparse GRM (exactly one)", nullptr};

inline const FlagDef kIndAfCoef = {
    "--ind-af-coef", "FILE", "Pre-computed individual AF model",
    "Produced by --cal-ind-af-coef.\n"
    "Text/gz: #CHROM  ID  STATUS  BETA0  BETA1 ...\n"
    "Binary (.bin): native format."};

inline const FlagDef kPairwiseIbd = {
    "--pairwise-ibd", "FILE", "Pairwise IBD file",
    "Tab-separated: #ID1  ID2  pa  pb  pc\n"
    "Produced by grab --cal-pairwise-ibd."};

// ── Scalar options ─────────────────────────────────────────────────

inline const FlagDef kPrevalence   = {"--prevalence",               "FLOAT",
    "Disease prevalence in reference population", nullptr};
inline const FlagDef kBatchPThresh = {"--batch-effect-p-threshold", "FLOAT",
    "Batch-effect p-value cutoff (default: 0.05)", nullptr};
inline const FlagDef kCovarPThresh = {"--covar-p-threshold",        "FLOAT",
    "Covariate p-value cutoff (default: 5e-5)", nullptr};
inline const FlagDef kSpaZThresh   = {"--spa-z-threshold",          "FLOAT",
    "SPA z-score cutoff (default: 2.0)", nullptr};
inline const FlagDef kOutlierIqr   = {"--outlier-iqr-threshold",    "FLOAT",
    "IQR outlier multiplier (default: 1.5)", nullptr};
inline const FlagDef kOutlierAbs   = {"--outlier-abs-bound",        "FLOAT",
    "Absolute outlier cutoff (default: 0.55)", nullptr};
inline const FlagDef kThreads      = {"--threads",    "INT",
    "Worker threads (default: 1)", nullptr};
inline const FlagDef kChunkSize    = {"--chunk-size", "INT",
    "Markers per chunk (default: 8192, min: 256)", nullptr};
inline const FlagDef kGeno         = {"--geno", "FLOAT",
    "Per-marker missing rate cutoff (default: 0.1)", nullptr};
inline const FlagDef kMaf          = {"--maf",  "FLOAT",
    "Min minor allele frequency (default: 1e-4)", nullptr};
inline const FlagDef kMac          = {"--mac",  "FLOAT",
    "Min minor allele count (default: 10)", nullptr};
inline const FlagDef kMinMafIbd    = {"--min-maf-ibd", "FLOAT",
    "Min MAF for IBD calculation (default: 0.01)", nullptr};
inline const FlagDef kNClusters    = {"--leaf-nclusters", "INT",
    "Number of K-means clusters for LEAF (default: from --ref-af count)", nullptr};
inline const FlagDef kSeed          = {"--seed", "INT",
    "Random seed for reproducibility (default: 0 = use random device)", nullptr};

inline const FlagDef kQuantPheno   = {"--pheno-quant", "COL",
    "Column name for quantitative phenotype (SPAsqr)", nullptr};
inline const FlagDef kSpasqrTaus   = {"--spasqr-taus", "LIST",
    "Comma-separated tau levels for SPAsqr (default: 0.1,0.3,0.5,0.7,0.9)", nullptr};

// ════════════════════════════════════════════════════════════════════
//  Method definitions
// ════════════════════════════════════════════════════════════════════

// ── SPACox ─────────────────────────────────────────────────────────
inline const FlagDef* const kSPACoxReq[] = {
    &kGeno_input, &kNullResid, &kOut, nullptr};
inline const FlagDef* const kSPACoxOpt[] = {
    &kCovar, &kCovarName, &kCovarPThresh,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kSPACox = {
    "SPACox",
    "Saddlepoint Approximation for Cox proportional hazards",
    kSPACoxReq, kSPACoxOpt,
    "IID  RESID [RESID2 ...]",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  SPACox_P  SPACox_Z",
    nullptr};

// ── SPAGRM ─────────────────────────────────────────────────────────
inline const FlagDef* const kSPAGRMReq[] = {
    &kGeno_input, &kNullResid, &kOut, &kSpGrm, &kPairwiseIbd, nullptr};
inline const FlagDef* const kSPAGRMOpt[] = {
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kSPAGRM = {
    "SPAGRM",
    "SPA with sparse GRM relatedness correction",
    kSPAGRMReq, kSPAGRMOpt,
    "IID  RESID [RESID2 ...]",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  SPAGRM_P  SPAGRM_Z",
    "Generate pairwise IBD file with: grab --cal-pairwise-ibd"};

// ── SPAmix ─────────────────────────────────────────────────────────
inline const FlagDef* const kSPAmixReq[] = {
    &kGeno_input, &kNullResid, &kOut, &kPcCols, nullptr};
inline const FlagDef* const kSPAmixOpt[] = {
    &kPheno, &kCovar, &kIndAfCoef, &kSpGrm, &kOutlierIqr,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kSPAmix = {
    "SPAmix",
    "SPA with individual ancestry-based allele frequencies",
    kSPAmixReq, kSPAmixOpt,
    "IID  RESID [RESID2 ...]",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  SPAmix_P  SPAmix_Z",
    "Pre-compute AF model for speed: grab --cal-ind-af-coef\n"
    "Optional --sp-grm-grab / --sp-grm-plink2 enables GRM-based variance (same as SPAmixPlus)."};

// ── SPAmixPlus ─────────────────────────────────────────────────────
inline const FlagDef* const kSPAmixPlusReq[] = {
    &kGeno_input, &kNullResid, &kOut, &kPcCols, &kSpGrm, nullptr};
inline const FlagDef* const kSPAmixPlusOpt[] = {
    &kPheno, &kCovar, &kIndAfCoef, &kOutlierIqr,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kSPAmixPlus = {
    "SPAmixPlus",
    "SPAmix with additional sparse GRM relatedness correction",
    kSPAmixPlusReq, kSPAmixPlusOpt,
    "IID  RESID [RESID2 ...]",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  SPAmixPlus_P  SPAmixPlus_Z",
    nullptr};

// ── SPAsqr ─────────────────────────────────────────────────────────
inline const FlagDef* const kSPAsqrReq[] = {
    &kGeno_input, &kOut, &kSpGrm, nullptr};
inline const FlagDef* const kSPAsqrOpt[] = {
    &kNullResid, &kPheno, &kCovar, &kCovarName, &kQuantPheno, &kSpasqrTaus,
    &kOutlierIqr, &kOutlierAbs,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kSPAsqr = {
    "SPAsqr",
    "SPA for quantile regression residuals (multiple tau levels)",
    kSPAsqrReq, kSPAsqrOpt,
    "IID  R_tau1  R_tau2  ...  R_tauK  (via --null-resid)\n"
    "  or --pheno + --pheno-quant + --spasqr-taus for built-in quantile regression",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P\n"
    "  P_CCT  P_tau{val}... Z_tau{val}...  (--pheno path)\n"
    "  Z_tau1 ... Z_tauK   P_tau1 ... P_tauK   P_CCT  (--null-resid path)",
    nullptr};

// ── WtCoxG ─────────────────────────────────────────────────────────
inline const FlagDef* const kWtCoxGReq[] = {
    &kGeno_input, &kOut, &kRefAf, &kPrevalence, nullptr};
inline const FlagDef* const kWtCoxGOpt[] = {
    &kNullResid, &kPheno, &kCovar, &kCovarName, &kBinaryPheno, &kSurvPheno,
    &kSpGrm, &kBatchPThresh,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kWtCoxG = {
    "WtCoxG",
    "Weighted Cox regression for time-to-event GWAS",
    kWtCoxGReq, kWtCoxGOpt,
    "IID  RESID  WEIGHT  INDICATOR  (via --null-resid) or --pheno + regression",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P\n"
    "  p_ext  p_noext  z_ext  z_noext  p_batch",
    nullptr};

// ── LEAF ───────────────────────────────────────────────────────────
inline const FlagDef* const kLEAFReq[] = {
    &kGeno_input, &kOut, &kRefAf, &kPrevalence, nullptr};
inline const FlagDef* const kLEAFOpt[] = {
    &kNullResid, &kPheno, &kCovar, &kCovarName, &kBinaryPheno, &kSurvPheno, &kPcCols,
    &kNClusters, &kSeed, &kSpGrm, &kBatchPThresh,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kLEAF = {
    "LEAF",
    "Local Ethnicity-Aware GWAS, WtCoxG per ancestry cluster",
    kLEAFReq, kLEAFOpt,
    "IID  RESID  WEIGHT  INDICATOR  (per cluster, via --null-resid) or --pheno + regression",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P\n"
    "  meta.p_ext  meta.p_noext  cl1.p_ext  cl1.p_noext  cl1.p_batch  [cl2 ...]",
    "--null-resid path: comma-sep per cluster; --pheno path: auto K-means on --pc-cols.\n"
    "The number of clusters and reference populations are independent.\n"
    "Summix estimates per-cluster ancestry proportions from the reference populations."};

// ── Utility mode: cal-ind-af-coef ──────────────────────────────────
inline const FlagDef* const kCalAfReq[] = {
    &kGeno_input, &kPcCols, &kOut, nullptr};
inline const FlagDef* const kCalAfOpt[] = {
    &kPheno, &kCovar, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kCalIndAfCoef = {
    "cal-ind-af-coef",
    "Compute per-marker individual-ancestry AF model",
    kCalAfReq, kCalAfOpt,
    nullptr,
    "Binary .bin file (pass to --ind-af-coef for SPAmix/SPAmixPlus).",
    nullptr};

// ── Utility mode: cal-pairwise-ibd ────────────────────────────────
inline const FlagDef* const kCalIbdReq[] = {
    &kGeno_input, &kOut, &kSpGrm, nullptr};
inline const FlagDef* const kCalIbdOpt[] = {
    &kMinMafIbd, nullptr};
inline const MethodDef kCalPairwiseIbd = {
    "cal-pairwise-ibd",
    "Compute pairwise IBD probabilities from sparse GRM + genotypes",
    kCalIbdReq, kCalIbdOpt,
    nullptr,
    "Tab-separated: #ID1  ID2  pa  pb  pc\n"
    "Pass to --pairwise-ibd for SPAGRM.",
    nullptr};

// ════════════════════════════════════════════════════════════════════
//  Lookup tables (null-terminated)
// ════════════════════════════════════════════════════════════════════

inline const MethodDef* const kAllMethods[] = {
    &kSPACox, &kSPAGRM, &kSPAmix, &kSPAmixPlus, &kSPAsqr, &kWtCoxG, &kLEAF,
    nullptr};

inline const MethodDef* const kAllUtilModes[] = {
    &kCalIndAfCoef, &kCalPairwiseIbd, nullptr};

// File-accepting flags (for --help <flag-topic>)
inline const FlagDef* const kFileFlags[] = {
    &kNullResid, &kPheno, &kCovar, &kRefAf,
    &kSpGrmGrab, &kSpGrmPlink2, &kIndAfCoef, &kPairwiseIbd,
    nullptr};

// All flags grouped for --help options
inline const FlagDef* const kInputFlags[] = {
    &kBfile, &kPfile, &kVcf, &kBgen,
    &kNullResid, &kOut,
    &kPheno, &kCovar, &kCovarName, &kBinaryPheno, &kSurvPheno, &kPcCols,
    &kRefAf,
    &kSpGrmGrab, &kSpGrmPlink2, &kIndAfCoef, &kPairwiseIbd,
    nullptr};

inline const FlagDef* const kNumericFlags[] = {
    &kPrevalence, &kBatchPThresh, &kCovarPThresh,
    &kSpaZThresh, &kOutlierIqr, &kOutlierAbs,
    &kThreads, &kChunkSize, &kNClusters, &kSeed, &kGeno, &kMaf, &kMac,
    &kMinMafIbd,
    nullptr};

} // namespace cli
