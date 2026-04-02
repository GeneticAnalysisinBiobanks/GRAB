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
    "--bfile", "PREFIX", "PLINK binary genotype prefix", nullptr};

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

// ── Method-specific file flags ─────────────────────────────────────

inline const FlagDef kCovar = {
    "--covar", "FILE", "Covariate file (plink2 .cov format)",
    "Formats: (A) header with IID/#IID, or (B) pure numeric matrix.\n"
    "  Format A: FID/SID/PAT/MAT/SEX/PHENO* columns auto-skipped.\n"
    "            Remaining column names do not matter; all used as covariates.\n"
    "  Format B: all columns are covariates; .fam order assumed.\n"
    "An intercept column is added automatically."};

inline const FlagDef kEigenvec = {
    "--eigenvec", "FILE", "Eigenvector file (plink2 .eigenvec)",
    "Formats: (A) header with IID/#IID, or (B) pure numeric matrix.\n"
    "  Format A: FID/SID skipped; PC columns identified by name (^[Pp][Cc]).\n"
    "            Unlike --null-resid/--covar, column names matter here.\n"
    "  Format B: all columns treated as PCs; .fam order assumed.\n"
    "Produced by plink2 --pca."};

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

// ════════════════════════════════════════════════════════════════════
//  Method definitions
// ════════════════════════════════════════════════════════════════════

// ── SPACox ─────────────────────────────────────────────────────────
inline const FlagDef* const kSPACoxReq[] = {
    &kBfile, &kNullResid, &kOut, nullptr};
inline const FlagDef* const kSPACoxOpt[] = {
    &kCovar, &kCovarPThresh,
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
    &kBfile, &kNullResid, &kOut, &kSpGrm, &kPairwiseIbd, nullptr};
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
    &kBfile, &kNullResid, &kOut, &kEigenvec, nullptr};
inline const FlagDef* const kSPAmixOpt[] = {
    &kIndAfCoef, &kOutlierIqr,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kSPAmix = {
    "SPAmix",
    "SPA with individual ancestry-based allele frequencies",
    kSPAmixReq, kSPAmixOpt,
    "IID  RESID [RESID2 ...]",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  SPAmix_P  SPAmix_Z",
    "Pre-compute AF model for speed: grab --cal-ind-af-coef"};

// ── SPAmixPlus ─────────────────────────────────────────────────────
inline const FlagDef* const kSPAmixPlusReq[] = {
    &kBfile, &kNullResid, &kOut, &kEigenvec, &kSpGrm, nullptr};
inline const FlagDef* const kSPAmixPlusOpt[] = {
    &kIndAfCoef, &kOutlierIqr,
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
    &kBfile, &kNullResid, &kOut, &kSpGrm, nullptr};
inline const FlagDef* const kSPAsqrOpt[] = {
    &kOutlierIqr, &kOutlierAbs,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kSPAsqr = {
    "SPAsqr",
    "SPA for quantile regression residuals (multiple tau levels)",
    kSPAsqrReq, kSPAsqrOpt,
    "IID  R_tau1  R_tau2  ...  R_tauK",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P\n"
    "  Z_tau1 ... Z_tauK   P_tau1 ... P_tauK   P_CCT",
    nullptr};

// ── WtCoxG ─────────────────────────────────────────────────────────
inline const FlagDef* const kWtCoxGReq[] = {
    &kBfile, &kNullResid, &kOut, &kRefAf, &kPrevalence, nullptr};
inline const FlagDef* const kWtCoxGOpt[] = {
    &kSpGrm, &kBatchPThresh,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kWtCoxG = {
    "WtCoxG",
    "Weighted Cox regression for time-to-event GWAS",
    kWtCoxGReq, kWtCoxGOpt,
    "IID  RESID  WEIGHT  INDICATOR",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P\n"
    "  p_ext  p_noext  z_ext  z_noext  p_batch",
    nullptr};

// ── LEAF ───────────────────────────────────────────────────────────
inline const FlagDef* const kLEAFReq[] = {
    &kBfile, &kNullResid, &kOut, &kRefAf, &kPrevalence, nullptr};
inline const FlagDef* const kLEAFOpt[] = {
    &kSpGrm, &kBatchPThresh,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kLEAF = {
    "LEAF",
    "Local Ethnicity-Aware GWAS, WtCoxG per ancestry cluster",
    kLEAFReq, kLEAFOpt,
    "IID  RESID  WEIGHT  INDICATOR  (per cluster file)",
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P\n"
    "  meta.p_ext  meta.p_noext  cl1.p_ext  cl1.p_noext  cl1.p_batch  [cl2 ...]",
    "Comma-separated --null-resid per ancestry cluster, --ref-af per ref-pop.\n"
    "The number of clusters and reference populations are independent.\n"
    "Summix estimates per-cluster ancestry proportions from the reference populations."};

// ── Utility mode: cal-ind-af-coef ──────────────────────────────────
inline const FlagDef* const kCalAfReq[] = {
    &kBfile, &kEigenvec, &kOut, nullptr};
inline const FlagDef* const kCalAfOpt[] = {
    &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, nullptr};
inline const MethodDef kCalIndAfCoef = {
    "cal-ind-af-coef",
    "Compute per-marker individual-ancestry AF model",
    kCalAfReq, kCalAfOpt,
    nullptr,
    "Binary .bin file (pass to --ind-af-coef for SPAmix/SPAmixPlus).",
    nullptr};

// ── Utility mode: cal-pairwise-ibd ────────────────────────────────
inline const FlagDef* const kCalIbdReq[] = {
    &kBfile, &kOut, &kSpGrm, nullptr};
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
    &kNullResid, &kCovar, &kEigenvec, &kRefAf,
    &kSpGrmGrab, &kSpGrmPlink2, &kIndAfCoef, &kPairwiseIbd,
    nullptr};

// All flags grouped for --help options
inline const FlagDef* const kInputFlags[] = {
    &kBfile, &kNullResid, &kOut,
    &kCovar, &kEigenvec, &kRefAf,
    &kSpGrmGrab, &kSpGrmPlink2, &kIndAfCoef, &kPairwiseIbd,
    nullptr};

inline const FlagDef* const kNumericFlags[] = {
    &kPrevalence, &kBatchPThresh, &kCovarPThresh,
    &kSpaZThresh, &kOutlierIqr, &kOutlierAbs,
    &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac,
    &kMinMafIbd,
    nullptr};

} // namespace cli
