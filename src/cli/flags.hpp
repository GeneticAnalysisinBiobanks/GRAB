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
    "--vcf", "FILE", "VCF genotype file (.vcf or BGZF-compressed .vcf.gz)",
    R"(Specifies a VCF (Variant Call Format) input file, either as plain text
or BGZF-compressed (typically produced by `bgzip` and conventionally named
`.vcf.gz`).  BGZF is the VCF specification's standard compression format;
plain gzip is also tolerated by the htslib backend but is not part of the
VCF standard.  zstd-compressed VCF is not supported on this path.
Matches plink2's --vcf flag: a BCF2 file passed here is rejected with a
message redirecting to --bcf.  Mutually exclusive with --bcf.)"
};

inline const FlagDef kBcf = {
    "--bcf", "FILE", "BCF2 genotype file",
    R"(Specifies a BCF2 (binary VCF) input file.  Matches plink2's --bcf flag:
a VCF text file passed here is rejected with a message redirecting to --vcf.
Mutually exclusive with --vcf.)"
};

inline const FlagDef kBgen = {
    "--bgen", "FILE <REF/ALT mode>",
    "BGEN genotype file (REF/ALT mode required, matches plink2 --bgen)",
    R"(BGEN itself does not encode which of the two listed alleles is REF, so
a REF/ALT mode must follow the filename.  Syntax follows plink2 --bgen:
  ref-first    The first allele for each variant is REF
               (raw UK Biobank / IMPUTE / qctool convention).
  ref-last     The last allele for each variant is REF
               (default for plink2-exported .bgen files; alleles[0] is ALT).
  ref-unknown  The last allele is treated as provisional-REF.  GRAB has no
               provisional-REF output column, so it behaves identically to
               ref-last and emits a [WARN] noting that REF is provisional.)"
};

// Combined display entry for genotype input (exactly one required)
inline const FlagDef kGeno_input = {
    "--bfile PREFIX | --pfile PREFIX | --vcf FILE | --bcf FILE | --bgen FILE <REF/ALT mode>",
    nullptr,
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
    R"(Selects columns from --covar as null-model covariates.
When --covar is absent, selects from --pheno instead.
An intercept is added automatically.
For SPAmix / SPAmixPlus: --pc-cols is *not* automatically added to
the null model; list any PCs you want adjusted in --covar-name too.)"
};

inline const FlagDef kPhenoName = {
    "--pheno-name", "COL_IDS", "Comma-separated phenotype column names",
    R"(Selects columns from --pheno for analysis.
Required by WtCoxG/LEAF (TIME:EVENT), SPAsqr (quantitative).)"
};

inline const FlagDef kResidName = {
    "--resid-name", "COL_IDS", "Comma-separated residual column names",
    R"(Selects columns from --pheno for use as residuals.
Required by SPACox, SPAGRM, SAGELD, SPAmix, SPAmixPlus, SPAmixLocalPlus.
Default: all residual columns.)"
};

inline const FlagDef kRegressionModel = {
    "--regression-model", "MODEL",
    "Null-model regression family: auto | linear | logistic | cox | ordinal",
    R"(Used with --pheno-name to fit a null model in-process and feed the
fitted residuals into the downstream score test (SPACox / SPAGRM /
SPAmix / SPAmixPlus).  Mutually exclusive with --resid-name.
  auto     — per-token inference from column values + colon syntax (default)
  linear   — regression::linearResiduals (intercept added automatically)
  logistic — regression::logisticResiduals (0/1 response; intercept added)
  cox      — regression::coxResiduals; --pheno-name uses TIME:EVENT pairs
             (e.g. SURVTIME1:STATUS1,SURVTIME2:STATUS2)
  ordinal  — regression::cumulativeLogitFit (integer-coded levels 0..J−1))"
};

inline const FlagDef kSaveResid = {
    "--save-resid", nullptr,
    "Write fitted residuals to PREFIX.null.resid (re-loadable via --resid-name)",
    R"(Only valid together with --pheno-name + --regression-model.  Writes a
plain-text strict-format file: IID column followed by one column per
fitted phenotype (NaN entries as 'NA').  The file can be reloaded in
a later invocation via --pheno PREFIX.null.resid --resid-name ....)"
};

inline const FlagDef kPcCols = {
    "--pc-cols", "COL_IDS", "Comma-separated PC column names (default: PC1,PC2,PC3,PC4)",
    R"(Selects columns from --covar or --pheno as principal components.
Used for K-means clustering in LEAF and for the per-individual AF
model in SPAmix / SPAmixPlus / --cal-af-coef.
Does NOT enter the null-model design.  To adjust the null model
for PCs as well, list them in --covar-name explicitly.)"
};

inline const FlagDef kRefAf = {
    "--ref-af", "FILE", "Reference allele frequency (plink2 .afreq or two-column numeric)",
    R"(Two accepted file formats:

(A) plink2 --freq output (.afreq).  Header line starts with '#' and lists
    the columns CHROM, ID, REF, ALT (or ALT1), ALT_FREQS (or ALT1_FREQ),
    OBS_CT in any order; column positions are inferred from the header.
    Records are matched to .bim by (CHROM, ID).  When the afreq REF/ALT
    pair is reversed relative to the .bim record, the loader flips the
    frequency to 1 − ALT_FREQS automatically.

(B) Two-column numeric fallback: ALT_FREQS  OBS_CT, one row per .bim
    record in .bim order, with no header.  The number of rows must equal
    the .bim marker count.

For --method LEAF: pass a comma-separated list of files, one per
reference population (e.g. --ref-af pop1.afreq,pop2.afreq).)"
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
    R"(plink2-style sparse GRM.  The file lists one related pair per line as
'ID1  ID2  VALUE', where VALUE is the kinship coefficient.  The companion
.grm.id file ('FID  IID' per subject) is auto-detected from the file path;
if .grm.id is absent, the loader assumes 0-based indices that match .fam
order.  Subjects not present in the genotype .fam are silently dropped.
Mutually exclusive with --sp-grm-grab.)"
};

// Combined display entry used in method required/optional arrays.
// --sp-grm-grab is also accepted by the parser but intentionally not
// advertised in --help; per-method usage prints only --sp-grm-plink2.
inline const FlagDef kSpGrm = {
    "--sp-grm-plink2", "FILE", "Sparse GRM (plink2 .grm.sp format)",
    nullptr
};

inline const FlagDef kIndAfCoef = {
    "--ind-af-coef", "FILE", "Pre-computed individual AF model",
    R"(Per-marker individual-ancestry allele-frequency model produced by
'grab2 --cal-af-coef'.  Plain text or gzip/zstd compressed; the loader
infers the codec from the file suffix (.gz, .zst).  Header line begins
with '#STATUS' and is followed by one BETA column per principal
component: '#STATUS  BETA0  BETA1  BETA2  ...'.  Rows are stored in
filtered-marker order (the order produced by --cal-af-coef on the same
genotype file with the same QC thresholds) and are matched
positionally — no marker ID is stored in the file, so the AF model
must be regenerated whenever the input genotype set or its QC filters
change.

AF coefficient scope.  The file stores a single AF model per marker,
fit on the subject set used during --cal-af-coef.  When passed to
SPAmix / SPAmixPlus with multiple phenotypes, every phenotype reuses
that one model regardless of its own missingness mask; this is the
intended cost of pre-computation.  By contrast, on-the-fly AF (no
--ind-af-coef) re-fits the AF model from each unique missingness
pattern present among the phenotypes in that run, so the model used
for a given phenotype depends only on that phenotype's non-missing
subjects.  In both modes, running a single phenotype alone and
running it jointly with other phenotypes produces byte-identical
per-marker results for that phenotype.)"
};

inline const FlagDef kPairwiseIbd = {
    "--pairwise-ibd", "FILE", "Pairwise IBD file",
    R"(Tab-separated: ID1  ID2  pa  pb  pc
Produced by grab2 --cal-pairwise-ibd.)"
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
    "--outlier-iqr-multiplier", "FLOAT",
    "IQR outlier multiplier for SPAmix/SPAmixPlus/SPAmixLocalPlus/WtCoxG/LEAF/SPAsqr/SPAGRM (default: 1.5)",
    R"(Tukey rule: a residual is an outlier iff r < q25 − k·IQR or r > q75 + k·IQR,
where k is the value of this flag.  Larger k → wider band → fewer outliers;
smaller k → narrower band → more outliers.  Outlier residuals enter the SPA
through their exact CGF; non-outlier residuals are folded into a second-
order Taylor (Gaussian) approximation.

  k = 1.5   (default)  ≈ 1 % outliers under a Gaussian residual distribution
  k = 1.0              ≈ 4 %
  k = 0                ≈ 50 % (band collapses to [q25, q75])
  k ≤ −0.5             100 % outliers (band collapses to a point at k = −0.5
                       and inverts for k < −0.5); equivalent to "all residuals
                       go through the exact CGF, no Taylor approximation".
                       This matches the R reference WtCoxG/LEAF implementation
                       byte-for-byte at the algorithmic level.
)"
};

inline const FlagDef kOutlierAbs = {
    "--spasqr-outlier-abs-bound", "FLOAT",
    "SPAsqr-only absolute outlier cutoff (default: 0.55)",
    nullptr
};

inline const FlagDef kSpagrmControlOutlier = {
    "--spagrm-control-outlier", nullptr,
    "Enable iterative IQR-ratio adjustment so the outlier share stays in (0, 5%] (SPAGRM, default: off)",
    R"(Mirrors the ControlOutlier argument of SPAGRM.NullModel() in the R reference.
Flag is parameterless: present  → enabled; absent → disabled (default).

  Enabled  — if 0 outliers are detected, shrink the IQR ratio by 0.8 until at
             least one outlier appears; if more than 5%% of subjects are
             outliers, increase the ratio by 0.5 until the share is at most 5%%.
  Disabled — keep the IQR ratio at --outlier-iqr-multiplier (default 1.5)
             without any adjustment.  Default; recommended when the residual
             distribution is well behaved.)"
};

inline const FlagDef kThreads = {
    "--threads", "INT", "Number of worker threads (default: 1)",
    R"(--threads N requests N worker threads for the chunk-level marker
engine and for the pre-marker work-stealing pools (null-model fits,
Chow-Liu MAF bins, K-means restarts, etc.).
The main thread is *not* counted in N: it dispatches chunks, drains
output, and stays mostly idle, but does briefly occupy one CPU during
load, finalize, and synchronization steps.  Plan for N + 1 logical
cores when sizing a job to physical cores.)"
};

inline const FlagDef kChunkSize = {
    "--chunk-size", "INT", "Markers per chunk (default: 8192, min: 256)",
    R"(Controls the granularity of the chunk-level work-stealing thread pool.
Each chunk is processed end-to-end by a single worker, then handed to a
single writer thread that emits chunks in genomic order.  Smaller chunks
improve load balancing on heterogeneous workloads at the cost of more
synchronization and per-chunk overhead; larger chunks reduce overhead
but may starve workers when the total marker count is small.

The default of 8192 suits whole-genome scans where each chunk amortises
the worker's startup cost over thousands of markers.  On a CLI-supplied
value, the engine honors it verbatim (subject to the min: 256 floor).

Exception — SPAsqr --spasqr-mode wald: per-marker QR refit is far slower
than score-mode batched GEMM, and wald runs are typically restricted to
a curated SNP list via --extract.  When --chunk-size is left at the
8192 default sentinel, wald auto-shrinks the chunk size to
  ceil( nMarkers / (4 * nthreads) )
so the chunk count is at least 4 * nthreads, keeping the worker pool
saturated even for very small marker sets (e.g. nMarkers = 10,
nthreads = 4 → chunk = 1, one marker per chunk).  An explicit
--chunk-size on the command line suppresses the auto-shrink.)"
};

inline const FlagDef kCompression = {
    "--compression", "gz|zst", "Compress output files (default: plain text)",
    nullptr
};

inline const FlagDef kCompressionLevel = {
    "--compression-level", "INT",
    "Compression level (gz: 1-9, default 6; zst: 1-22, default 3)",
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

inline const FlagDef kLeafClusterFile = {
    "--leaf-cluster-file", "FILE",
    "Pre-computed cluster labels for LEAF (skip K-means)",
    R"(Header line required.  Reads two columns by name: IID (or #IID) and cluster.
All other columns in the file are ignored.  Cluster values must be integers
in {1, …, K}; K is inferred as max(cluster) and cross-checked against
--leaf-nclusters when both are provided.)"
};

inline const FlagDef kLeafKmeansNstart = {
    "--leaf-kmeans-nstart", "INT",
    "Number of K-means++ random restarts for LEAF (default: 25). "
    "Ignored when --leaf-cluster-file is supplied.",
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
    "--pred-list", "FILE", "pred.list file for LOCO analysis (phenoName locoPath per line; Regenie or LDAK-KVIK)",
    R"(Space-separated file with two columns: phenotype name and path to
the corresponding .loco file. The .loco file format is auto-detected
per phenotype: Regenie (header begins with FID_IID, chromosome-major
rows) or LDAK-KVIK (header FID IID Chr1 Chr2 ... Chr22, subject-major
rows).)"
};

inline const FlagDef kSpasqrSolver = {
    "--spasqr-solver", "NAME",
    "SPAsqr null-model solver: qmme (default) | conquer",
    R"(Selects the smoothed quantile regression solver for the SPAsqr null model:
  qmme    — Quadratic Majorization-Minimization with Extrapolation
            (Heng & Wang, 2025). Caches Cholesky of the Hessian upper bound
            once per phenotype × bandwidth and reuses it across all τ levels;
            far more robust on ill-conditioned X. Default.
  conquer — Convolution-type smoothed QR (He et al., 2021): Huber init
            followed by BB gradient descent. Refits from scratch per τ.)"
};

inline const FlagDef kSpasqrMode = {
    "--spasqr-mode", "MODE",
    "SPAsqr inference mode: score (default) | wald",
    R"(Selects the inference framework used by --method SPAsqr:
  score — score test on null-model residuals (default).  One null QR fit per
          phenotype × τ; markers are streamed and tested via the score
          statistic S = Σ R_i G_i with M-estimation sandwich variance.
          Output (per marker): CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P
          P_CCT P_tau{val}... Z_tau{val}...
  wald  — full-model Wald test.  For every (marker, τ), the joint smoothed-QR
          model with [X | G] is refit by QMME and β̂_G + SE are computed from
          the (γ,γ) entry of the M-estimation sandwich V = A^{-1} B A^{-1}/n.
          Slower per marker; suited for follow-up effect-size estimation on
          a small SNP list (--extract).  Per-marker QR refit runs on the
          shared marker-engine thread pool (--threads), and --chunk-size
          auto-shrinks at its 8192 default so the pool stays fed even on
          small marker sets — see --chunk-size for details.  Output is
          plink2-style one-marker-per-line wide format:
          CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P
          P_CCT P_tau{val}... Z_tau{val}... BETA_tau{val}... SE_tau{val}...
          (--pred-list gives y_resp = Y − loco_chr; --compression honored.))"
};

inline const FlagDef kPhenoTransform = {
    "--pheno-transform", "MODE",
    "Phenotype pre-transform for SPAsqr: raw | int | standardize "
    "(default: int)",
    R"(Selects how Y is transformed before SQR fitting (and before LOCO offset
subtraction when --pred-list is given):
  raw         — no transform; SQR fits on Y as supplied.
  int         — inverse normal transform (Blom plotting position, average-rank
                ties) applied per phenotype on its non-missing scope. Default
                in both contexts; LOCO PRS should be trained on INT Y so
                the offset is on the same scale.
  standardize — Y is centered and scaled to unit variance per phenotype.
                Default in both contexts; LOCO PRS should be trained on
                a standardized Y so the offset is on the same scale.

When --pred-list is provided, the transformed Y has loco_chr subtracted as
an offset (β=1, α=0). The LOCO PRS scale must match the chosen transform —
mixing scales gives meaningless residuals.)"
};

inline const FlagDef kExtract = {
    "--extract", "FILE", "Restrict analysis to SNPs listed in file (one ID per line)",
    R"(Single-column file of SNP IDs.  Comparison is by marker ID
(.bim column 2, .pvar ID column, BGEN RSID/SNPID, VCF/BCF ID field).
Applied uniformly to all genotype inputs: --bfile, --pfile, --vcf,
--bcf, --bgen, and --admix-bfile.  IDs listed in the file that do not
match any marker are silently ignored.)"
};

inline const FlagDef kExclude = {
    "--exclude", "FILE", "Exclude SNPs listed in file from analysis (one ID per line)",
    R"(Single-column file of SNP IDs.  Comparison is by marker ID
(.bim column 2, .pvar ID column, BGEN RSID/SNPID, VCF/BCF ID field).
Applied uniformly to all genotype inputs: --bfile, --pfile, --vcf,
--bcf, --bgen, and --admix-bfile.  IDs listed in the file that do not
match any marker are silently ignored.)"
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
Used with --vcf or --bcf to produce admixed .abed ancestry tracks.)"
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

inline const FlagDef kSageldX = {
    "--sageld-x", "COL_IDS",
    "Comma-separated environment column name(s) for SAGELD G x E (pheno-input mode)",
    R"(Required for SAGELD's pheno-input mode together with --pheno-name and
--covar-name.  Each name must match a numeric column of --pheno; the env
column also has to be listed in --covar-name so it enters the fixed-effect
design.  For every (phenotype, env) pair the null model
    Y ~ X + (E | IID)            (random intercept + random slope on E)
is fit by EM-ML internally, and the BLUP residuals are aggregated to per-
IID (R_G, R_<E>, R_Gx<E>) before the marker-level G and G x E score tests
run.  Multiple envs trigger a separate model per env.)"
};

inline const FlagDef kSpasqrTol = {
    "--spasqr-tol", "FLOAT",
    "Convergence tolerance for the SPAsqr null-model SQR solver (default: 1e-7). "
    "QMME tightens this internally to min(--spasqr-tol, 1e-9).",
    nullptr
};

inline const FlagDef kSpasqrH = {
    "--spasqr-h", "FLOAT",
    "Explicit bandwidth for the SQR null model (mutually exclusive with --spasqr-h-scale)",
    nullptr
};

inline const FlagDef kSpasqrHScale = {
    "--spasqr-h-scale", "FLOAT",
    "Divisor for IQR-based bandwidth: h = IQR(Y) / SCALE  "
    "(default: 3 in score mode, 10 in --spasqr-mode wald; mutually exclusive with --spasqr-h)",
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
    &kCovar,       &kCovarName,        &kResidName,    &kPhenoName,   &kRegressionModel,  &kSaveResid,
    &kCovarPThresh, &kSpaZThresh,      &kSeed,
    &kThreads,      &kChunkSize,
    &kCompression, &kCompressionLevel,
    &kKeep,         &kRemove,           &kExtract,      &kExclude,
    &kGeno,         &kMaf,         &kMac,        &kHwe,     &kChr,
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
    &kResidName,  &kPhenoName,  &kRegressionModel,           &kSaveResid,
    &kCovar,      &kCovarName,
    &kSpaZThresh, &kOutlierIqr, &kSpagrmControlOutlier,      &kSeed,
    &kThreads, &kChunkSize, &kCompression, &kCompressionLevel,
    &kKeep,       &kRemove,  &kExtract,    &kExclude,
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
    "Generate pairwise IBD file with: grab2 --cal-pairwise-ibd",
};

// ── SAGELD ─────────────────────────────────────────────────────────
inline const FlagDef *const kSAGELDReq[] = {
    &kGeno_input, &kPheno, &kOut, &kSpGrm, &kPairwiseIbd,
    nullptr
};
inline const FlagDef *const kSAGELDOpt[] = {
    &kResidName, &kPhenoName, &kCovarName, &kSageldX,    &kSaveResid,
    &kSpaZThresh, &kThreads, &kChunkSize, &kCompression, &kCompressionLevel,
    &kKeep,       &kRemove,  &kExtract,   &kExclude,
    &kGeno, &kMaf, &kMac, &kHwe, &kChr,
    nullptr
};
inline const MethodDef kSAGELD = {
    "SAGELD",
    "G x E interaction analysis for longitudinal data with GRM correction",
    kSAGELDReq,
    kSAGELDOpt,
    "Residual mode: --pheno FILE --resid-name R_G,R_<E1>,R_Gx<E1>[,...] "
    "(file format: #IID  R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...]).  "
    "Pheno mode: --pheno LONG_FILE --pheno-name Y1,Y2,... --covar-name X1,X2,... --sageld-x E1[,E2,...]",
    R"(Residual mode (single file): PREFIX.SAGELD
Pheno mode (per phenotype):  PREFIX.<phenoName>.SAGELD
  CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
  P_G  P_Gx<E1>  [...]  Z_G  Z_Gx<E1>  [...])",
    R"(Two input modes (mutually exclusive):
  Residual mode — supply lme4::lmer() residuals directly via --resid-name.
                  Column layout: R_G followed by (R_<E>, R_Gx<E>) pairs.
  Pheno mode    — supply long-format Y, X, E and fit  Y ~ X + (E | IID)
                  internally via EM-ML.  --covar-name must include every
                  variable in --sageld-x.

Generate the IBD file once with: grab2 --cal-pairwise-ibd)",
};

// ── SPAmix ─────────────────────────────────────────────────────────
inline const FlagDef *const kSPAmixReq[] = {
    &kGeno_input, &kPheno, &kOut, &kPcCols,
    nullptr
};
inline const FlagDef *const kSPAmixOpt[] = {
    &kCovar,      &kCovarName,  &kResidName,    &kPhenoName,    &kRegressionModel, &kSaveResid,
    &kIndAfCoef, &kSpGrm,       &kOutlierIqr,
    &kSpaZThresh, &kSeed,
    &kThreads, &kChunkSize, &kCompression, &kCompressionLevel,
    &kKeep,       &kRemove,  &kExtract,   &kExclude,
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
    R"(Pre-compute the AF model for speed: grab2 --cal-af-coef.
Optional --sp-grm-plink2 enables GRM-based variance correction (residual mode
absorbs the GRM via the variance ratio).

AF coefficient scope.  With --ind-af-coef, every phenotype in this run
shares the single pre-computed AF model per marker (fit at --cal-af-coef
time on its full subject set).  Without --ind-af-coef, the AF model
is re-fit on the fly for each unique missingness pattern present
among the phenotypes, so each phenotype's AF is computed only from
its own non-missing subjects.  Either way, running one phenotype
alone and running it together with other phenotypes yields the same
per-marker output for that phenotype.)",
};

// ── SPAmixPlus ─────────────────────────────────────────────────────
inline const FlagDef *const kSPAmixPlusReq[] = {
    &kGeno_input, &kPheno, &kOut, &kPcCols, &kSpGrm,
    nullptr
};
inline const FlagDef *const kSPAmixPlusOpt[] = {
    &kCovar,      &kCovarName,  &kResidName,  &kPhenoName,    &kRegressionModel,    &kSaveResid,
    &kIndAfCoef,        &kOutlierIqr, &kSpaZThresh, &kSeed,    &kThreads,
    &kChunkSize, &kCompression, &kCompressionLevel,
    &kKeep,       &kRemove,    &kExtract,    &kExclude,
    &kGeno,       &kMaf,        &kMac,
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
    R"(AF coefficient scope.  With --ind-af-coef, every phenotype in this
run shares the single pre-computed AF model per marker (fit at
--cal-af-coef time on its full subject set).  Without --ind-af-coef,
the AF model is re-fit on the fly for each unique missingness pattern
present among the phenotypes, so each phenotype's AF is computed only
from its own non-missing subjects.  Either way, running one phenotype
alone and running it together with other phenotypes yields the same
per-marker output for that phenotype.)",
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
    &kSpaZThresh,   &kThreads,    &kChunkSize,
    &kCompression,  &kCompressionLevel,
    &kKeep,         &kRemove,     &kExtract,    &kExclude,
    &kGeno, &kMaf,
    &kMac,          &kHwe,        &kChr,        &kPredList,    &kPhenoTransform,
    &kSpasqrSolver, &kSpasqrMode,
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
    &kPhenoName,  &kRegressionModel,    &kSpGrm,  &kBatchPThresh,
    &kSpaZThresh, &kOutlierIqr, &kThreads, &kChunkSize,
    &kCompression, &kCompressionLevel,
    &kKeep,       &kRemove, &kExtract,   &kExclude,
    &kGeno, &kMaf,
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
    &kRegressionModel,
    &kPcCols,    &kNClusters, &kLeafClusterFile, &kLeafKmeansNstart,
    &kSeed,      &kSpGrm,     &kBatchPThresh, &kSpaZThresh,
    &kOutlierIqr,
    &kThreads,   &kChunkSize,
    &kCompression, &kCompressionLevel,
    &kKeep,      &kRemove,    &kExtract,   &kExclude,
    &kGeno,      &kMaf,       &kMac,          &kHwe,        &kChr,
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
    &kKeep,        &kRemove,           &kExtract, &kExclude,
    &kThreads,     &kCompression,      &kCompressionLevel,
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
    &kCovar,            &kCovarName,        &kResidName,  &kPhenoName,  &kRegressionModel, &kSaveResid,
    &kKeep,             &kRemove,           &kExtract,    &kExclude,
    &kOutlierIqr,       &kSpaZThresh,       &kThreads,    &kChunkSize,
    &kCompression,      &kCompressionLevel,
    &kGeno,             &kMaf,              &kMac,        &kHwe,        &kChr,
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
    &kVcf, &kBcf, &kMsp, &kAdmixTextPrefix,
    &kKeep, &kRemove, &kThreads,
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
  --vcf FILE --rfmix-msp FILE  phased VCF + rfmix2 MSP            ->  .abed
  --bcf FILE --rfmix-msp FILE  phased BCF2 + rfmix2 MSP           ->  .abed
  --admix-text-prefix PREFIX   extract_tracts text output         ->  .abed
--out PREFIX writes PREFIX.abed, PREFIX.bim, PREFIX.fam.
Pass PREFIX as --admix-bfile to SPAmixLocalPlus or --cal-phi.)",
};

// ── Utility mode: int-pheno ────────────────────────────────────────
inline const FlagDef *const kIntPhenoReq[] = {
    &kPheno, &kOut,
    nullptr
};

inline const FlagDef *const kIntPhenoOpt[] = {
    &kPhenoName,
    nullptr
};

inline const MethodDef kIntPheno = {
    "int-pheno",
    "Inverse-normal-transform selected trait columns in a phenotype file",
    kIntPhenoReq,
    kIntPhenoOpt,
    nullptr,
    "Output: PREFIX.txt  (ID columns + selected Y columns; missing entries → NA)",
    R"(Reads --pheno FILE under GRAB's standard header conventions
    #IID col...        IID col...
    #FID IID col...    FID IID col...
and writes PREFIX.txt preserving the input ID-column layout.  When
--pheno-name is supplied, only the listed trait columns are
INT-transformed and emitted; when omitted, every trait column in the
input is transformed.  Each retained column is independently
INT-transformed (Blom plotting position, average-rank ties) on its own
non-missing scope; missing entries (NA / "." / blank) stay missing in
the output.)",
};

// ── Utility mode: cal-af-coef ──────────────────────────────────────
inline const FlagDef *const kCalAfReq[] = {
    &kGeno_input, &kPcCols, &kOut,
    nullptr
};

inline const FlagDef *const kCalAfOpt[] = {
    &kPheno,   &kCovar,     &kKeep, &kRemove,
    &kExtract, &kExclude,
    &kCompression, &kCompressionLevel,
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
    R"(One AF model per marker, fit on the subject set defined by
--bfile / --pfile / --bgen / --vcf intersected with --keep / --remove
and the marker QC thresholds (--geno / --maf / --mac / --hwe / --chr).
The file is independent of any later phenotype: every phenotype that
consumes it via --ind-af-coef shares the same per-marker model.

Use this mode when (i) the same genotype + subject set will be reused
across many phenotypes / runs and the AF re-fit cost is significant,
and (ii) it is acceptable to compute each marker's AF on the full
--cal-af-coef subject set rather than on each phenotype's
non-missing subset.  When phenotype-specific missingness should
drive the AF model (e.g. some phenotypes drop a substantial fraction
of subjects), omit --ind-af-coef and let SPAmix / SPAmixPlus re-fit
the AF on the fly per unique missingness pattern; in that mode each
phenotype's AF is computed only from its own non-missing subjects.
In both modes the per-phenotype output is invariant to whether the
phenotype is run alone or jointly with other phenotypes.)",
};

// ── Utility mode: cal-pairwise-ibd ────────────────────────────────
inline const FlagDef *const kCalIbdReq[] = {
    &kGeno_input, &kOut, &kSpGrm,
    nullptr
};

inline const FlagDef *const kCalIbdOpt[] = {
    &kKeep,     &kRemove,
    &kExtract,  &kExclude,  &kChr,
    &kMinMafIbd,
    &kThreads,  &kCompression, &kCompressionLevel,
    nullptr
};

inline const MethodDef kCalPairwiseIbd = {
    "cal-pairwise-ibd",
    "Compute pairwise IBD probabilities from sparse GRM + genotypes",
    kCalIbdReq,
    kCalIbdOpt,
    nullptr,
    R"(Output: PREFIX.ibd[.gz|.zst]
Tab-separated: ID1  ID2  pa  pb  pc
Pass to --pairwise-ibd for SPAGRM.)",
    nullptr,
};

// ════════════════════════════════════════════════════════════════════
//  Lookup tables (null-terminated)
// ════════════════════════════════════════════════════════════════════

// kAllMethods / kAllUtilModes drive method-name canonicalization inside the
// dispatcher and must therefore continue to list every method GRAB knows
// how to run, including the ones that are intentionally hidden from
// --help (SPAmixLocalPlus, --make-abed, --cal-phi).  Help generation uses
// the kVisible* arrays below, which omit the hidden entries; running
// `grab --help SPAmixLocalPlus|make-abed|cal-phi` therefore reports
// "Unknown help topic" while the methods themselves remain functional.
inline const MethodDef *const kAllMethods[] = {
    &kSPACox, &kSPAGRM, &kSAGELD, &kSPAmix, &kSPAmixPlus, &kSPAmixLocalPlus,
    &kSPAsqr, &kWtCoxG, &kLEAF,
    nullptr
};

inline const MethodDef *const kAllUtilModes[] = {
    &kCalAfCoef, &kCalPairwiseIbd, &kCalPhi, &kMakeAbed, &kIntPheno,
    nullptr
};

// Visible methods exclude SPAmixPlus and SPAmixLocalPlus.  Both remain
// callable via --method but are hidden from --help to keep the surface
// area focussed on the seven core GWAS methods.
inline const MethodDef *const kVisibleMethods[] = {
    &kSPACox, &kSPAGRM, &kSAGELD, &kSPAmix,
    &kSPAsqr, &kWtCoxG, &kLEAF,
    nullptr
};

inline const MethodDef *const kVisibleUtilModes[] = {
    &kCalAfCoef, &kCalPairwiseIbd, &kIntPheno,
    nullptr
};

// File-accepting flags (for --help <flag-topic>).  Admix-* topics and
// --sp-grm-grab are omitted because the methods that consume them
// (SPAmixLocalPlus, --cal-phi, --make-abed) and the legacy GRAB sparse
// GRM format are hidden from --help; the flags themselves remain
// accepted by the parser.
inline const FlagDef *const kFileFlags[] = {
    &kPheno,       &kCovar,     &kRefAf,
    &kSpGrmPlink2, &kIndAfCoef, &kPairwiseIbd,
    nullptr
};

// Flags shown by `--help options`.  Admix-mode flags (--admix-bfile,
// --admix-phi, --rfmix-msp, --admix-text-prefix) and --sp-grm-grab are
// excluded because the methods that consume them (SPAmixLocalPlus,
// --cal-phi, --make-abed) and the legacy GRAB sparse-GRM format are
// hidden from --help.
inline const FlagDef *const kInputFlags[] = {
    &kBfile,       &kPfile,       &kVcf,         &kBcf,
    &kBgen,
    &kOut,         &kCompression, &kCompressionLevel,
    &kPheno,       &kCovar,       &kCovarName,
    &kPhenoName,   &kResidName,   &kRegressionModel, &kSaveResid,
    &kPcCols,      &kRefAf,
    &kSpGrmPlink2, &kIndAfCoef,   &kPairwiseIbd,
    &kPredList,    &kPhenoTransform,
    &kLeafClusterFile,
    &kKeep,        &kRemove,
    &kExtract,     &kExclude,     &kChr,
    nullptr
};

inline const FlagDef *const kNumericFlags[] = {
    &kPrevalence, &kBatchPThresh, &kCovarPThresh,     &kSpaZThresh, &kOutlierIqr, &kOutlierAbs,
    &kSpagrmControlOutlier,
    &kThreads,    &kChunkSize,    &kCompressionLevel, &kNClusters,
    &kLeafKmeansNstart,
    &kSeed,       &kGeno,
    &kMaf,        &kMac,          &kHwe,              &kMinMafIbd,
    &kSpasqrTaus, &kSpasqrTol,    &kSpasqrH,          &kSpasqrHScale,
    &kSpasqrSolver, &kSpasqrMode,
    &kSageldX,
    nullptr
};

inline const FlagDef *const kFilterFlags[] = {
    &kExtract, &kExclude, &kChr,
    nullptr
};

} // namespace cli
