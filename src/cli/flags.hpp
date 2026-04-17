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
    const char *flag;     // "--threads", or combined display string for groups
    const char *metavar;  // "INT"/"FLOAT"/"FILE"/"PREFIX", nullptr for bool/group
    const char *brief;    // one-line description
    const char *fileInfo; // extended file-format help, or nullptr
};

// ── Method / utility-mode definition ───────────────────────────────

struct MethodDef {
    const char *name;               // "SPACox", "cal-ind-af-coef", ...
    const char *desc;               // one-line description
    const FlagDef *const *required; // null-terminated
    const FlagDef *const *optional; // null-terminated
    const char *residNote;          // per-method residual columns, or nullptr
    const char *outputCols;         // output column description
    const char *extra;              // additional notes, or nullptr
};

// ════════════════════════════════════════════════════════════════════
//  Flag instances
// ════════════════════════════════════════════════════════════════════

// ── Shared required ────────────────────────────────────────────────

inline const FlagDef kBfile = {
    "--bfile", "PREFIX", "PLINK binary genotype prefix (.bed/.bim/.fam)",
    nullptr
};

inline const FlagDef kPfile = {
    "--pfile", "PREFIX", "PLINK2 genotype prefix (.pgen/.pvar/.psam)",
    nullptr
};

inline const FlagDef kVcf = {
    "--vcf", "FILE", "VCF or BCF genotype file",
    nullptr
};

inline const FlagDef kBgen = {
    "--bgen", "FILE", "BGEN genotype file",
    nullptr
};

// Combined display entry for genotype input (exactly one required)
inline const FlagDef kGeno_input = {
    "--bfile PREFIX | --pfile PREFIX | --vcf FILE | --bgen FILE", nullptr,
    "Genotype input (exactly one)",
    nullptr
};

inline const FlagDef kOut = {
    "--out", "PREFIX", "Output prefix (suffix added per mode/method)",
    nullptr
};

inline const FlagDef kCovar = {
    "--covar", "FILE", "Covariate file (strict format, columns selected by --covar-name)",
    R"(Strict format: mandatory header line, first column = subject ID.
  All column names must match [0-9A-Za-z_\-.]+.
Use --covar-name to select specific columns as covariates.
If --covar-name is omitted, all data columns are used.
When --covar is absent, --covar-name selects columns from --pheno instead.)"
};

inline const FlagDef kPheno = {
    "--pheno", "FILE", "Phenotype file (strict format, first column = subject ID)",
    R"(Strict format: mandatory header line, first column = subject ID.
  All column names must match [0-9A-Za-z_\-.]+.
Use --pheno-name to select specific phenotype columns.
Use --resid-name to select pre-computed residual columns.)"
};

inline const FlagDef kCovarName = {
    "--covar-name", "COL_IDS", "Comma-separated covariate column names",
    R"(Selects columns from --covar as covariates.
When --covar is absent, selects from --pheno instead.
An intercept is added automatically.)"
};

inline const FlagDef kPhenoName = {
    "--pheno-name", "COL_IDS", "Comma-separated phenotype column names",
    R"(Selects columns from --pheno for analysis.
Required by POLMM (ordinal), WtCoxG/LEAF (TIME:EVENT), SPAsqr (quantitative).)"
};

inline const FlagDef kResidName = {
    "--resid-name", "COL_IDS", "Comma-separated residual column names",
    R"(Selects columns from --pheno for use as residuals.
Required by SPACox, SPAGRM, SAGELD, SPAmix, SPAmixPlus, SPAmixLocalPlus.
Default: all residual columns.)"
};

inline const FlagDef kPcCols = {
    "--pc-cols", "COL_IDS", "Comma-separated PC column names (default: PC1,PC2,PC3,PC4)",
    R"(Selects columns from --covar or --pheno as principal components.
Used for K-means clustering in LEAF and AF estimation in SPAmix.)"
};

inline const FlagDef kRefAf = {
    "--ref-af", "FILE", "Reference allele frequency (plink2 .afreq or numeric)",
    R"(Formats: (A) plink2 .afreq: #CHROM  ID  REF  ALT  ALT_FREQS  OBS_CT
           Matched to .bim by (CHROM, ID) with allele flip detection.
         (B) Two-column numeric: ALT_FREQS  OBS_CT (rows in .bim order).
For LEAF: comma-separated, one file per reference population.)"
};

inline const FlagDef kSpGrmGrab = {
    "--sp-grm-grab", "FILE", "Sparse GRM, GRAB format (ID1 ID2 VALUE)",
    R"(Whitespace-delimited with mandatory header: ID1  ID2  VALUE
(header may have # prefix).  #-comment lines skipped.
Subjects not in .fam are silently dropped.
Mutually exclusive with --sp-grm-plink2.)"
};

inline const FlagDef kSpGrmPlink2 = {
    "--sp-grm-plink2", "FILE", "Sparse GRM, plink2 .grm.sp format",
    R"(Companion .grm.id (FID  IID) auto-detected from file path.
If .grm.id is absent, 0-based indices assumed to match .fam order.
Mutually exclusive with --sp-grm-grab.)"
};

// Combined display entry used in method required/optional arrays
inline const FlagDef kSpGrm = {
    "--sp-grm-grab FILE | --sp-grm-plink2 FILE", nullptr, "Sparse GRM (exactly one)",
    nullptr
};

inline const FlagDef kIndAfCoef = {
    "--ind-af-coef", "FILE", "Pre-computed individual AF model",
    R"(Produced by --cal-ind-af-coef.
Text/gz: #STATUS  BETA0  BETA1 ...
Rows in filtered-marker order (positional matching).)"
};

inline const FlagDef kPairwiseIbd = {
    "--pairwise-ibd", "FILE", "Pairwise IBD file",
    R"(Tab-separated: #ID1  ID2  pa  pb  pc
Produced by grab --cal-pairwise-ibd.)"
};

// ── Scalar options ─────────────────────────────────────────────────

inline const FlagDef kPrevalence = {
    "--prevalence", "FLOAT", "Disease prevalence in reference population",
    nullptr
};

inline const FlagDef kBatchPThresh = {
    "--batch-effect-p-threshold", "FLOAT",
    "Batch-effect p-value cutoff (default: 0.05)",
    nullptr
};

inline const FlagDef kCovarPThresh = {
    "--covar-p-threshold", "FLOAT", "Covariate p-value cutoff (default: 5e-5)",
    nullptr
};
inline const FlagDef kSpaZThresh = {
    "--spa-z-threshold", "FLOAT", "SPA z-score cutoff (default: 2.0)",
    nullptr
};

inline const FlagDef kOutlierIqr = {
    "--outlier-iqr-threshold", "FLOAT", "IQR outlier multiplier (default: 1.5)",
    nullptr
};

inline const FlagDef kOutlierAbs = {
    "--outlier-abs-bound", "FLOAT", "Absolute outlier cutoff (default: 0.55)",
    nullptr
};

inline const FlagDef kThreads = {
    "--threads", "INT", "Worker threads (default: 1)",
    nullptr
};

inline const FlagDef kChunkSize = {
    "--chunk-size", "INT", "Markers per chunk (default: 8192, min: 256)",
    nullptr
};

inline const FlagDef kCompression = {
    "--compression", "gz|zst", "Compress output files (default: plain text)",
    nullptr
};

inline const FlagDef kCompressionLevel = {
    "--compression-level", "INT", "Compression level (gz: 1-9, zst: 1-22; default: 0 = library default)",
    nullptr
};

inline const FlagDef kGeno = {
    "--geno", "FLOAT", "Per-marker missing rate cutoff (default: 0.1)",
    nullptr
};

inline const FlagDef kMaf = {
    "--maf", "FLOAT", "Min minor allele frequency (default: 1e-5)",
    nullptr
};

inline const FlagDef kMac = {
    "--mac", "FLOAT", "Min minor allele count (default: 10)",
    nullptr
};

inline const FlagDef kHwe = {
    "--hwe", "FLOAT", "Exclude markers with HWE p < threshold (default: 0, disabled)",
    nullptr
};

inline const FlagDef kMinMafIbd = {
    "--min-maf-ibd", "FLOAT", "Min MAF for IBD calculation (default: 0.01)",
    nullptr
};

inline const FlagDef kNClusters = {
    "--leaf-nclusters", "INT",
    "Number of K-means clusters for LEAF (default: from --ref-af count)",
    nullptr
};

inline const FlagDef kKeep = {
    "--keep", "FILE", "Restrict analysis to subjects listed in file (FID IID per line)",
    nullptr
};

inline const FlagDef kRemove = {
    "--remove", "FILE", "Exclude subjects listed in file (FID IID per line)",
    nullptr
};

inline const FlagDef kPredList = {
    "--pred-list", "FILE", "Regenie step 1 pred.list file for LOCO analysis (phenoName locoPath per line)",
    R"(Space-separated file with two columns: phenotype name and path to
the corresponding .loco file from Regenie step 1.)"
};

inline const FlagDef kExtract = {
    "--extract", "FILE", "Restrict analysis to SNPs listed in file (one ID per line)",
    R"(Single-column file of SNP IDs (matching .bim column 2).
Applied to both --bfile and --admix-bfile inputs.)"
};

inline const FlagDef kExclude = {
    "--exclude", "FILE", "Exclude SNPs listed in file from analysis (one ID per line)",
    R"(Single-column file of SNP IDs (matching .bim column 2).
Applied to both --bfile and --admix-bfile inputs.)"
};

inline const FlagDef kChr = {
    "--chr", "NUMS", "Restrict analysis to specified chromosomes",
    R"(Comma-separated chromosome numbers or ranges.
Examples: --chr 5   --chr 2,3   --chr 1-4,6-8,22)"
};

inline const FlagDef kAdmixBfile = {
    "--admix-bfile", "PREFIX",
    "Admixed ancestry binary genotype prefix (.abed/.bim/.fam)",
    R"(Binary format storing 2K tracks (dosage + hapcount per ancestry).
Shares standard PLINK .fam and .bim files.
The .abed file has a 16-byte header with magic 0xAD4D, version,
number of ancestries K, number of subjects N, and markers M.)"
};

inline const FlagDef kAdmixPhi = {
    "--admix-phi", "FILE", "Pre-computed phi file from --cal-phi (wide format)",
    R"(Tab-separated: idx1  idx2  anc0_A  anc0_B  anc0_C  anc0_D  [anc1_A ...]
Indices are 0-based into .fam row order. One row per related pair.)"
};

inline const FlagDef kMsp = {
    "--rfmix-msp", "FILE", "MSP local-ancestry file from rfmix2 (for --make-abed)",
    R"(rfmix2 output format:
  Line 1: #Subpopulation order/codes: 0=POP0\t1=POP1\t...  (K inferred)
  Line 2: #chm\tspos\tepos\tsgpos\tegpos\tn snps\tIID0.0\tIID0.1\t...
  Data:   chrom\tspos(0-based)\tepos(excl)\t...\tancestry_calls...
Used with --vcf to produce admixed .abed ancestry tracks.)"
};

inline const FlagDef kAdmixTextPrefix = {
    "--admix-text-prefix", "PREFIX", "extract_tracts text output prefix (for --make-abed)",
    R"({PREFIX}.anc{k}.dosage[.gz]   and   {PREFIX}.anc{k}.hapcount[.gz]  (k=0,1,...)
Header: CHROM  POS  ID  REF  ALT  SAMPLE1  SAMPLE2  ...
Values: integer 0-2 per subject; K auto-detected from file presence.)"};
;

inline const FlagDef kSeed = {
    "--seed", "INT", "Random seed for reproducibility (default: 0 = use random device)",
    nullptr
};

inline const FlagDef kSpasqrTaus = {
    "--spasqr-taus", "LIST",
    "Comma-separated tau levels for SPAsqr, max 20 (default: 0.1,0.3,0.5,0.7,0.9)",
    nullptr
};

inline const FlagDef kSpasqrTol = {
    "--spasqr-tol", "FLOAT",
    "Convergence tolerance for conquer quantile regression (default: 1e-7)",
    nullptr
};

inline const FlagDef kSpasqrH = {
    "--spasqr-h", "FLOAT",
    "Explicit bandwidth for conquer (mutually exclusive with --spasqr-h-scale)",
    nullptr
};

inline const FlagDef kSpasqrHScale = {
    "--spasqr-h-scale", "FLOAT",
    "Divisor for IQR-based bandwidth: h = IQR(Y) / SCALE (default: 3; mutually exclusive with --spasqr-h)",
    nullptr
};


// ════════════════════════════════════════════════════════════════════
//  Method definitions
// ════════════════════════════════════════════════════════════════════

// ── SPACox ─────────────────────────────────────────────────────────
inline const FlagDef *const kSPACoxReq[] = {
    &kGeno_input, &kPheno, &kOut,
    nullptr
};
inline const FlagDef *const kSPACoxOpt[] = {
    &kCovar,       &kCovarName,        &kResidName,    &kCovarPThresh, &kSpaZThresh, &kThreads, &kChunkSize,
    &kCompression, &kCompressionLevel, &kGeno,         &kMaf,         &kMac,        &kHwe,     &kChr,
    nullptr
};
inline const MethodDef kSPACox = {
    "SPACox",
    "Saddlepoint Approximation for Cox proportional hazards",
    kSPACoxReq,
    kSPACoxOpt,
    "#IID  RESID  [RESID2  ...] (from --pheno, selected by --resid-name)",
    R"(Per residual column: PREFIX.PHENO.SPACox[.gz|.zst]
  CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z)",
    nullptr,
};

// ── SPAGRM ─────────────────────────────────────────────────────────
inline const FlagDef *const kSPAGRMReq[] = {
    &kGeno_input, &kPheno, &kOut, &kSpGrm, &kPairwiseIbd,
    nullptr
};
inline const FlagDef *const kSPAGRMOpt[] = {
    &kSpaZThresh, &kThreads, &kChunkSize, &kCompression, &kCompressionLevel,
    &kGeno,       &kMaf,     &kMac,       &kHwe,         &kChr,
    nullptr
};
inline const MethodDef kSPAGRM = {
    "SPAGRM",
    "SPA with sparse GRM relatedness correction",
    kSPAGRMReq,
    kSPAGRMOpt,
    "#IID  RESID  [RESID2  ...] (from --pheno, selected by --resid-name)",
    R"(Per residual column: PREFIX.PHENO.SPAGRM[.gz|.zst]
  CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z)",
    "Generate pairwise IBD file with: grab --cal-pairwise-ibd",
};

// ── SAGELD ─────────────────────────────────────────────────────────
inline const FlagDef *const kSAGELDReq[] = {
    &kGeno_input, &kPheno, &kOut, &kSpGrm, &kPairwiseIbd,
    nullptr
};
inline const FlagDef *const kSAGELDOpt[] = {
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf, &kMac, &kHwe, &kChr,
    nullptr
};
inline const MethodDef kSAGELD = {
    "SAGELD",
    "G x E interaction analysis for longitudinal data with GRM correction",
    kSAGELDReq,
    kSAGELDOpt,
    "#IID  R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...] (from --pheno)",
    R"(CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
  P_G  P_Gx<E1>  [...]  Z_G  Z_Gx<E1>  [...])",
    "Residuals from lme4::lmer(); generate IBD with: grab --cal-pairwise-ibd",
};

// ── SPAmix ─────────────────────────────────────────────────────────
inline const FlagDef *const kSPAmixReq[] = {
    &kGeno_input, &kPheno, &kOut, &kPcCols,
    nullptr
};
inline const FlagDef *const kSPAmixOpt[] = {
    &kPheno,      &kCovar,   &kIndAfCoef, &kSpGrm,       &kOutlierIqr,
    &kSpaZThresh, &kThreads, &kChunkSize, &kCompression, &kCompressionLevel,
    &kGeno,       &kMaf,     &kMac,       &kHwe,         &kChr,
    nullptr
};
inline const MethodDef kSPAmix = {
    "SPAmix",
    "SPA with individual ancestry-based allele frequencies",
    kSPAmixReq,
    kSPAmixOpt,
    "#IID  RESID  [RESID2  ...] (from --pheno, selected by --resid-name)",
    R"(Per residual column: PREFIX.PHENO.SPAmix[.gz|.zst]
  CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z  BETA  SE)",
    R"(Pre-compute AF model for speed: grab --cal-af-coef
Optional --sp-grm-grab / --sp-grm-plink2 enables GRM-based variance (same as SPAmixPlus).)",
};

// ── SPAmixPlus ─────────────────────────────────────────────────────
inline const FlagDef *const kSPAmixPlusReq[] = {
    &kGeno_input, &kPheno, &kOut, &kPcCols, &kSpGrm,
    nullptr
};
inline const FlagDef *const kSPAmixPlusOpt[] = {
    &kCovar,     &kIndAfCoef,        &kOutlierIqr, &kSpaZThresh, &kThreads,
    &kChunkSize, &kCompression, &kCompressionLevel, &kGeno,       &kMaf,        &kMac,
    &kHwe,       &kChr,
    nullptr
};
inline const MethodDef kSPAmixPlus = {
    "SPAmixPlus",
    "SPAmix with additional sparse GRM relatedness correction",
    kSPAmixPlusReq,
    kSPAmixPlusOpt,
    "#IID  RESID  [RESID2  ...] (from --pheno, selected by --resid-name)",
    R"(Per residual column: PREFIX.PHENO.SPAmixP[.gz|.zst]
  CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z  BETA  SE)",
    nullptr,
};

// ── POLMM ──────────────────────────────────────────────────────────
inline const FlagDef *const kPOLMMReq[] = {
    &kGeno_input, &kPheno, &kPhenoName, &kOut, &kSpGrm,
    nullptr
};

inline const FlagDef *const kPOLMMOpt[] = {
    &kCovar,   &kCovarName, &kSpaZThresh,
    &kThreads, &kChunkSize, &kGeno,     &kMaf,
    &kMac,     &kHwe,       &kChr,
    nullptr
};

inline const MethodDef kPOLMM = {
    "POLMM",
    "Proportional Odds Logistic Mixed Model for ordinal categorical GWAS",
    kPOLMMReq,
    kPOLMMOpt,
    nullptr,
    "CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  POLMM_P  POLMM_Z  POLMM_BETA  POLMM_SE",
    R"(Fits null model internally (all-in-one). Ordinal levels auto-detected.
Uses sparse GRM for random effects via PCG solver.)",
};

// ── SPAsqr ─────────────────────────────────────────────────────────
inline const FlagDef *const kSPAsqrReq[] = {
    &kGeno_input, &kOut, &kSpGrm,
    nullptr
};

inline const FlagDef *const kSPAsqrOpt[] = {
    &kPheno,        &kCovar,      &kCovarName,  &kResidName,
    &kPhenoName,    &kSpasqrTaus, &kSpasqrTol,  &kSpasqrH,
    &kSpasqrHScale, &kOutlierIqr, &kOutlierAbs,
    &kSpaZThresh,   &kThreads,    &kChunkSize,  &kGeno, &kMaf,
    &kMac,          &kHwe,        &kChr,        &kPredList,
    nullptr
};

inline const MethodDef kSPAsqr = {
    "SPAsqr",
    "SPA for quantile regression residuals (multiple tau levels)",
    kSPAsqrReq,
    kSPAsqrOpt,
    "--pheno + --pheno-name + --spasqr-taus for built-in quantile regression",
    R"(CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
  P_CCT  P_tau{val}... Z_tau{val}...)",
    nullptr,
};

// ── WtCoxG ─────────────────────────────────────────────────────────
inline const FlagDef *const kWtCoxGReq[] = {
    &kGeno_input, &kOut, &kRefAf, &kPrevalence,
    nullptr
};

inline const FlagDef *const kWtCoxGOpt[] = {
    &kPheno,      &kCovar,  &kCovarName, &kResidName,
    &kPhenoName,  &kSpGrm,  &kBatchPThresh,
    &kSpaZThresh, &kThreads, &kChunkSize, &kGeno, &kMaf,
    &kMac,        &kHwe,    &kChr,
    nullptr
};

inline const MethodDef kWtCoxG = {
    "WtCoxG",
    "Weighted Cox regression for time-to-event GWAS",
    kWtCoxGReq,
    kWtCoxGOpt,
    "--pheno + --pheno-name (TIME:EVENT pair or binary)",
    R"(CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
  p_ext  p_noext  z_ext  z_noext  p_batch)",
    nullptr,
};

// ── LEAF ───────────────────────────────────────────────────────────
inline const FlagDef *const kLEAFReq[] = {
    &kGeno_input, &kOut, &kRefAf, &kPrevalence,
    nullptr
};

inline const FlagDef *const kLEAFOpt[] = {
    &kPheno,     &kCovar,     &kCovarName, &kResidName, &kPhenoName,
    &kPcCols,    &kNClusters, &kSeed,      &kSpGrm,     &kBatchPThresh, &kSpaZThresh,
    &kThreads,   &kChunkSize, &kGeno,      &kMaf,       &kMac,          &kHwe,        &kChr,
    nullptr
};

inline const MethodDef kLEAF = {
    "LEAF",
    "Local Ethnicity-Aware GWAS, WtCoxG per ancestry cluster",
    kLEAFReq,
    kLEAFOpt,
    "--pheno + --pheno-name + --pc-cols for auto K-means clustering",
    R"(CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
  meta.p_ext  meta.p_noext  cl1.p_ext  cl1.p_noext  cl1.p_batch  [cl2 ...])",
    R"(--pheno path: auto K-means on --pc-cols.
The number of clusters and reference populations are independent.
Summix estimates per-cluster ancestry proportions from the reference populations.)",
};

// ── Utility mode: cal-phi ───────────────────────────────────────────
inline const FlagDef *const kCalPhiReq[] = {
    &kAdmixBfile, &kSpGrm, &kOut,
    nullptr
};

inline const FlagDef *const kCalPhiOpt[] = {
    &kKeep,        &kRemove,           &kExtract, &kExclude, &kChr,
    &kCompression, &kCompressionLevel,
    nullptr
};

inline const MethodDef kCalPhi = {
    "cal-phi",
    "Compute phi kinship matrices from admixed genotypes + sparse GRM",
    kCalPhiReq,
    kCalPhiOpt,
    nullptr,
    R"(Wide tab-separated: idx1  idx2  anc0_A  anc0_B  anc0_C  anc0_D  [anc1_A ...]
Output: PREFIX.phi[.gz|.zst])",
    "Pass output to --admix-phi for SPAmixLocalPlus GWAS.",
};

// ── SPAmixLocalPlus (── --pheno + --admix-bfile + --admix-phi) ──
// Not in kAllMethods—dispatched automatically when the three flags are present.
inline const FlagDef *const kSPAmixLocalPlusReq[] = {
    &kAdmixBfile, &kPheno, &kAdmixPhi, &kOut,
    nullptr
};

inline const FlagDef *const kSPAmixLocalPlusOpt[] = {
    &kExtract,          &kExclude, &kOutlierIqr, &kSpaZThresh, &kThreads, &kChunkSize, &kCompression,
    &kCompressionLevel, &kGeno,    &kMaf,        &kMac,        &kHwe,     &kChr,
    nullptr
};

inline const MethodDef kSPAmixLocalPlus = {
    "SPAmixLocalPlus",
    "Local-ancestry-specific GWAS with SPA + phi-based variance correction",
    kSPAmixLocalPlusReq,
    kSPAmixLocalPlusOpt,
    "IID  RESID [RESID2 ...] (from --pheno, selected by --resid-name)",
    R"(CHROM  POS  ID  REF  ALT  P_CCT
  anc0_AltFreq  anc0_MissingRate  anc0_P  anc0_Pnorm  anc0_Stat  anc0_Var  anc0_zScore  anc0_AltCounts  anc0_BetaG
  anc1_AltFreq  ...  (repeated for each ancestry k))",
    R"(Two-phase workflow:
  1. grab --cal-phi --admix-bfile PREFIX --sp-grm-plink2 FILE --out OUTPUT_PREFIX
  2. grab --method SPAmixLocalPlus --admix-bfile PREFIX --admix-phi OUTPUT_PREFIX.phi --pheno FILE --out PREFIX
Output: PREFIX.PHENO.LocalP[.gz|.zst] per residual column)",
};

// ── Utility mode: make-abed ────────────────────────────────────────
inline const FlagDef *const kMakeAbedReq[] = {
    &kOut,
    nullptr
};

inline const FlagDef *const kMakeAbedOpt[] = {
    &kVcf, &kMsp, &kAdmixTextPrefix,
    nullptr
};

inline const MethodDef kMakeAbed = {
    "make-abed",
    "Build .abed admixed ancestry binary from VCF+MSP or extract_tracts output",
    kMakeAbedReq,
    kMakeAbedOpt,
    nullptr,
    "{prefix}.abed  {prefix}.bim  {prefix}.fam",
    R"(Two modes (mutually exclusive):
  --vcf FILE --rfmix-msp FILE  phased VCF/BCF + rfmix2 MSP  ->  .abed
  --admix-text-prefix PREFIX   extract_tracts text output   ->  .abed
--out PREFIX writes PREFIX.abed, PREFIX.bim, PREFIX.fam.
Pass PREFIX as --admix-bfile to SPAmixLocalPlus or --cal-phi.)",
};

// ── Utility mode: cal-af-coef ──────────────────────────────────────
inline const FlagDef *const kCalAfReq[] = {
    &kGeno_input, &kPcCols, &kOut,
    nullptr
};

inline const FlagDef *const kCalAfOpt[] = {
    &kPheno,   &kCovar,     &kKeep, &kRemove, &kCompression, &kCompressionLevel,
    &kThreads, &kChunkSize, &kGeno, &kMaf,    &kMac,         &kHwe,             &kChr,
    nullptr
};

inline const MethodDef kCalAfCoef = {
    "cal-af-coef",
    "Compute per-marker individual-ancestry AF model",
    kCalAfReq,
    kCalAfOpt,
    nullptr,
    "Output: PREFIX.afc[.gz|.zst] (pass to --ind-af-coef for SPAmix/SPAmixPlus)",
    nullptr,
};

// ── Utility mode: cal-pairwise-ibd ────────────────────────────────
inline const FlagDef *const kCalIbdReq[] = {
    &kGeno_input, &kOut, &kSpGrm,
    nullptr
};

inline const FlagDef *const kCalIbdOpt[] = {
    &kKeep, &kRemove, &kMinMafIbd, &kCompression, &kCompressionLevel,
    nullptr
};

inline const MethodDef kCalPairwiseIbd = {
    "cal-pairwise-ibd",
    "Compute pairwise IBD probabilities from sparse GRM + genotypes",
    kCalIbdReq,
    kCalIbdOpt,
    nullptr,
    R"(Output: PREFIX.ibd[.gz|.zst]
Tab-separated: #ID1  ID2  pa  pb  pc
Pass to --pairwise-ibd for SPAGRM.)",
    nullptr,
};

// ════════════════════════════════════════════════════════════════════
//  Lookup tables (null-terminated)
// ════════════════════════════════════════════════════════════════════

inline const MethodDef *const kAllMethods[] = {
    &kSPACox, &kSPAGRM, &kSAGELD, &kSPAmix, &kSPAmixPlus, &kSPAmixLocalPlus,
    &kPOLMM,  &kSPAsqr, &kWtCoxG, &kLEAF,
    nullptr
};

inline const MethodDef *const kAllUtilModes[] = {
    &kCalAfCoef, &kCalPairwiseIbd, &kCalPhi, &kMakeAbed,
    nullptr
};

// File-accepting flags (for --help <flag-topic>)
inline const FlagDef *const kFileFlags[] = {
    &kPheno,        &kCovar,       &kRefAf,    &kSpGrmGrab,
    &kSpGrmPlink2, &kIndAfCoef, &kPairwiseIbd, &kAdmixPhi,
    nullptr
};

// All flags grouped for --help options
inline const FlagDef *const kInputFlags[] = {
    &kBfile,       &kPfile,       &kVcf,
    &kBgen,        &kAdmixBfile,
    &kOut,         &kCompression, &kCompressionLevel,
    &kPheno,       &kCovar,       &kCovarName,
    &kPhenoName,   &kResidName,
    &kPcCols,      &kRefAf,       &kSpGrmGrab,
    &kSpGrmPlink2, &kIndAfCoef,   &kPairwiseIbd,
    &kAdmixPhi,    &kMsp,         &kAdmixTextPrefix,
    &kExtract,     &kExclude,     &kChr,
    nullptr
};

inline const FlagDef *const kNumericFlags[] = {
    &kPrevalence, &kBatchPThresh, &kCovarPThresh,     &kSpaZThresh, &kOutlierIqr, &kOutlierAbs,
    &kThreads,    &kChunkSize,    &kCompressionLevel, &kNClusters,  &kSeed,       &kGeno,
    &kMaf,        &kMac,          &kMinMafIbd,
    nullptr
};

inline const FlagDef *const kFilterFlags[] = {
    &kExtract, &kExclude, &kChr,
    nullptr
};

} // namespace cli
