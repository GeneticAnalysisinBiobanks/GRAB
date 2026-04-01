// main.cpp — GRAB CLI entry point
// Parses command-line arguments and delegates to the selected method workflow.

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "util/logging.hpp"
#include "spacox/spacox.hpp"
#include "spamix/spamix.hpp"
#include "spamix/spamixplus.hpp"
#include "spamix/indiv_af.hpp"
#include "wtcoxg/wtcoxg.hpp"
#include "wtcoxg/leaf.hpp"
#include "spasqr/spasqr.hpp"
#include "spagrm/ibd.hpp"
#include "spagrm/gt_prob.hpp"

namespace {

struct Args {
  std::string method;
  std::string helpTopic;          // set when --help <topic> is used
  std::string residFile;
  std::string covarFile;
  std::string eigenVecsFile;
  std::string bfilePrefix;
  std::string refAfFile;
  std::string spGrmSaigeFile;     // --sp-grm-saige
  std::string spGrmPlink2Prefix;  // --sp-grm-plink2  (PREFIX → .grm.sp + .grm.id)
  std::string indAfFile;          // --ind-af-coef
  std::string pairwiseIBDFile;
  std::string outputFile;
  bool   calIndAfCoef      = false;  // --cal-ind-af-coef
  bool   calPairwiseIBD    = false;
  double minMafIBD         = 0.01;
  double refPrevalence     = -1.0;
  double cutoff            = 0.05;
  double spaCutoff         = 2.0;
  double pvalCovAdjCut     = 5e-5;
  double missingCutoff     = 0.1;
  double minMafCutoff      = 1e-4;
  double minMacCutoff      = 10.0;
  double outlierRatio      = 1.5;
  double outlierAbsBound   = 0.55;
  int    nthread           = 1;
  int    nSnpPerChunk      = 8192;
};

void printShortHelp() {
  std::cerr <<
    "GRAB — Genome-wide Regression Analysis for Biobanks\n"
    "\n"
    "Usage:\n"
    "  grab --method METHOD  --bfile PREFIX  --null-resid FILE  --out FILE  [OPTIONS]\n"
    "  grab --cal-ind-af-coef  --bfile PREFIX  --eigenvec FILE  --out FILE  [OPTIONS]\n"
    "  grab --cal-pairwise-ibd  (--sp-grm-saige FILE | --sp-grm-plink2 PREFIX)\n"
    "       --bfile PREFIX  --out FILE  [OPTIONS]\n"
    "\n"
    "METHOD is one of: SPACox  SPAGRM  SPAmix  SPAmixPlus  SPAsqr  WtCoxG  LEAF\n"
    "\n"
    "Per-method required inputs (* = one of --sp-grm-saige / --sp-grm-plink2):\n"
    "  SPACox      --null-resid  [--covar]\n"
    "  SPAGRM      --null-resid  --sp-grm-*  --pairwise-ibd\n"
    "  SPAmix      --null-resid  --eigenvec  [--ind-af-coef]\n"
    "  SPAmixPlus  --null-resid  --eigenvec  --sp-grm-*  [--ind-af-coef]\n"
    "  SPAsqr      --null-resid  --sp-grm-*\n"
    "  WtCoxG      --null-resid  --ref-af    --prevalence  [--sp-grm-*]\n"
    "  LEAF        --null-resid  --ref-af    --prevalence  [--sp-grm-*]\n"
    "              (LEAF: comma-separated lists for --null-resid and --ref-af)\n"
    "\n"
    "Shared options (all methods):\n"
    "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
    "  --spa-z-threshold FLOAT\n"
    "\n"
    "Run 'grab --help <topic>' for details.  Topics:\n"
    "  SPACox  SPAGRM  SPAmix  SPAmixPlus  SPAsqr  WtCoxG  LEAF\n"
    "  cal-ind-af-coef  cal-pairwise-ibd\n"
    "  null-resid  sp-grm  options\n";
}

void printTopicHelp(const std::string& topic) {
  if (topic == "null-resid") {
    std::cerr <<
      "--null-resid FILE\n"
      "  Null model residual file.  plink2 --pheno-compatible format:\n"
      "    - Optional ##-prefixed comment lines before the header\n"
      "    - Header line starting with #FID, FID, #IID, or IID\n"
      "    - If header present: FID IID col(s) skipped, remaining cols are residual values\n"
      "    - Legacy headerless format (col 0 = IID, rest = values) also accepted\n"
      "  Multi-column residual files: each value column runs as a separate GWAS.\n"
      "  Single column → output to --out file directly.\n"
      "  Multiple columns → outputs PREFIX.1.gz, PREFIX.2.gz, ... (always gzip).\n"
      "  Column layout by method:\n"
      "    SPACox / SPAGRM / SPAmix / SPAmixPlus : IID  RESID [RESID2 ...]\n"
      "    WtCoxG / LEAF                          : IID  RESID  WEIGHT  INDICATOR\n"
      "    SPAsqr                                 : IID  R_tau1  R_tau2  ...  R_tauK\n";
    return;
  }
  if (topic == "sp-grm") {
    std::cerr <<
      "Sparse GRM input — provide exactly one of:\n"
      "\n"
      "--sp-grm-saige FILE\n"
      "  SAIGE sparse GRM format.  Whitespace-delimited, '#'-lines skipped.\n"
      "  Columns: ID1  ID2  VALUE\n"
      "\n"
      "--sp-grm-plink2 PREFIX\n"
      "  plink2 --make-grm-sparse output.  Reads PREFIX.grm.sp and PREFIX.grm.id.\n"
      "  PREFIX.grm.id : FID  IID per line (defines 0-based subject order)\n"
      "  PREFIX.grm.sp : idx1  idx2  value  (0-based indices, lower-triangle)\n";
    return;
  }
  if (topic == "options") {
    std::cerr <<
      "Shared options (all GWAS methods):\n"
      "  --threads INT               Number of worker threads (default: 1)\n"
      "  --chunk-size INT            Markers per chunk (default: 8192, min: 256)\n"
      "  --geno FLOAT                Per-marker missing rate cutoff (default: 0.1)\n"
      "  --maf  FLOAT                Min minor allele frequency (default: 1e-4)\n"
      "  --mac  FLOAT                Min minor allele count (default: 10)\n"
      "  --spa-z-threshold FLOAT     SPA saddle-point z-score cutoff (default: 2.0)\n"
      "\n"
      "Method-specific options:\n"
      "  --covar FILE                Covariate file, plink2 .cov format (SPACox)\n"
      "                              Tab/space-delimited; header with IID col required\n"
      "                              FID/SID/PAT/MAT/SEX cols skipped automatically\n"
      "                              Remaining numeric cols used as covariates\n"
      "                              An intercept column is added automatically\n"
      "  --covar-p-threshold FLOAT   p-value threshold for covariate adjustment (default: 5e-5, SPACox)\n"
      "  --ref-af FILE               plink2 --freq .afreq file (WtCoxG, LEAF)\n"
      "                              Columns: #CHROM ID REF ALT ALT_FREQS OBS_CT\n"
      "                              For LEAF: comma-separated, one per cluster\n"
      "  --prevalence FLOAT          Disease prevalence in reference population (WtCoxG, LEAF)\n"
      "  --batch-effect-p-threshold FLOAT  Batch-effect test p-value cutoff (default: 0.05, WtCoxG/LEAF)\n"
      "  --eigenvec FILE             plink2 --pca .eigenvec file (SPAmix, SPAmixPlus)\n"
      "                              Header: #FID IID PC1 PC2 ...\n"
      "  --ind-af-coef FILE          Pre-computed individual AF model (SPAmix, SPAmixPlus)\n"
      "                              Produced by --cal-ind-af-coef\n"
      "                              If omitted, AF is computed on-the-fly from --eigenvec\n"
      "                              Text/gz cols: #CHROM ID STATUS BETA0 BETA1 ...\n"
      "                              Binary (.bin): produced by --cal-ind-af-coef\n"
      "  --pairwise-ibd FILE         Pairwise IBD file (SPAGRM); produced by --cal-pairwise-ibd\n"
      "  --outlier-iqr-threshold FLOAT  IQR outlier multiplier (default: 1.5, SPAmix/SPAmixPlus/SPAsqr)\n"
      "  --outlier-abs-bound FLOAT   Absolute outlier cutoff clamp (default: 0.55, SPAsqr)\n";
    return;
  }
  if (topic == "SPACox") {
    std::cerr <<
      "Method: SPACox — Saddlepoint Approximation for Cox proportional hazards\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --null-resid FILE   --out FILE\n"
      "\n"
      "Optional:\n"
      "  --covar FILE                plink2 .cov covariate file\n"
      "  --covar-p-threshold FLOAT   p-value cutoff for covariate correction (default: 5e-5)\n"
      "  --spa-z-threshold FLOAT     SPA z-score cutoff (default: 2.0)\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output columns:\n"
      "  Marker  Chr  Pos  Ref  Alt  AltCount  AltFreq  MissRate  HWEpval  zScore  Pvalue\n";
    return;
  }
  if (topic == "SPAGRM") {
    std::cerr <<
      "Method: SPAGRM — SPA accounting for genetic relatedness via a sparse GRM\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --null-resid FILE   --out FILE\n"
      "  --sp-grm-saige FILE  |  --sp-grm-plink2 PREFIX   (exactly one)\n"
      "  --pairwise-ibd FILE\n"
      "\n"
      "Optional:\n"
      "  --spa-z-threshold FLOAT   (default: 2.0)\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output columns:\n"
      "  Marker  Chr  Pos  Ref  Alt  AltCount  AltFreq  MissRate  HWEpval  zScore  Pvalue\n"
      "\n"
      "Hint: generate the pairwise IBD file with: grab --cal-pairwise-ibd\n";
    return;
  }
  if (topic == "SPAmix") {
    std::cerr <<
      "Method: SPAmix — SPA with individual ancestry-based allele frequencies\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --null-resid FILE   --out FILE\n"
      "  --eigenvec FILE\n"
      "\n"
      "Optional:\n"
      "  --ind-af-coef FILE          Pre-computed individual AF model (recommended for speed)\n"
      "  --outlier-iqr-threshold FLOAT  (default: 1.5)\n"
      "  --spa-z-threshold FLOAT   (default: 2.0)\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output columns:\n"
      "  Marker  Chr  Pos  Ref  Alt  AltCount  AltFreq  MissRate  HWEpval  Pvalue  zScore\n"
      "\n"
      "Hint: pre-compute the AF model with: grab --cal-ind-af-coef\n";
    return;
  }
  if (topic == "SPAmixPlus") {
    std::cerr <<
      "Method: SPAmixPlus — SPAmix with additional sparse GRM relatedness correction\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --null-resid FILE   --out FILE\n"
      "  --eigenvec FILE\n"
      "  --sp-grm-saige FILE  |  --sp-grm-plink2 PREFIX   (exactly one)\n"
      "\n"
      "Optional:\n"
      "  --ind-af-coef FILE          Pre-computed individual AF model\n"
      "  --outlier-iqr-threshold FLOAT  (default: 1.5)\n"
      "  --spa-z-threshold FLOAT   (default: 2.0)\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output columns:\n"
      "  Marker  Chr  Pos  Ref  Alt  AltCount  AltFreq  MissRate  HWEpval  Pvalue  zScore\n";
    return;
  }
  if (topic == "SPAsqr") {
    std::cerr <<
      "Method: SPAsqr — SPA for quantile regression residuals (multiple tau levels)\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --null-resid FILE   --out FILE\n"
      "  --sp-grm-saige FILE  |  --sp-grm-plink2 PREFIX   (exactly one)\n"
      "\n"
      "  --null-resid format: IID  R_tau1  R_tau2  ...  R_tauK\n"
      "    Each column is the residual at one quantile level tau.\n"
      "\n"
      "Optional:\n"
      "  --outlier-iqr-threshold FLOAT  (default: 1.5)\n"
      "  --outlier-abs-bound FLOAT      (default: 0.55)\n"
      "  --spa-z-threshold FLOAT        (default: 2.0)\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output columns:\n"
      "  Marker  Chr  Pos  Ref  Alt  AltCount  AltFreq  MissRate  HWEpval\n"
      "  Z_tau1 ... Z_tauK   P_tau1 ... P_tauK   P_CCT\n";
    return;
  }
  if (topic == "WtCoxG") {
    std::cerr <<
      "Method: WtCoxG — Weighted Cox regression for time-to-event GWAS\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --null-resid FILE   --out FILE\n"
      "  --ref-af FILE    --prevalence FLOAT\n"
      "\n"
      "  --null-resid format: IID  RESID  WEIGHT  INDICATOR\n"
      "\n"
      "Optional:\n"
      "  --sp-grm-saige FILE  |  --sp-grm-plink2 PREFIX\n"
      "  --batch-effect-p-threshold FLOAT   (default: 0.05)\n"
      "  --spa-z-threshold FLOAT            (default: 2.0)\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output columns:\n"
      "  Marker  Chr  Pos  Ref  Alt  AltCount  AltFreq  MissRate  HWEpval\n"
      "  WtCoxG.ext  WtCoxG.noext  zscore.ext  zscore.noext  AF_ref  N_ref  pvalue_bat\n";
    return;
  }
  if (topic == "LEAF") {
    std::cerr <<
      "Method: LEAF — Local Ethnicity-Aware GWAS, WtCoxG per ancestry cluster\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX\n"
      "  --null-resid FILE1,FILE2,...   one residual file per cluster\n"
      "  --ref-af   FILE1,FILE2,...     one .afreq file per cluster (same order)\n"
      "  --prevalence FLOAT\n"
      "  --out FILE\n"
      "\n"
      "  --null-resid format (each file): IID  RESID  WEIGHT  INDICATOR\n"
      "\n"
      "Optional:\n"
      "  --sp-grm-saige FILE  |  --sp-grm-plink2 PREFIX\n"
      "  --batch-effect-p-threshold FLOAT   (default: 0.05)\n"
      "  --spa-z-threshold FLOAT            (default: 2.0)\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output columns:\n"
      "  Marker  Chr  Pos  Ref  Alt  AltCount  AltFreq  MissRate  HWEpval\n"
      "  meta.p_ext  meta.p_noext\n"
      "  cl1.p_ext  cl1.p_noext  cl1.p_batch  [cl2 ...]  ...\n";
    return;
  }
  if (topic == "cal-ind-af-coef") {
    std::cerr <<
      "Mode: --cal-ind-af-coef\n"
      "  Compute per-marker individual-ancestry allele-frequency model coefficients.\n"
      "  Pre-computes the AF model used by SPAmix / SPAmixPlus.\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --eigenvec FILE   --out FILE\n"
      "\n"
      "Optional:\n"
      "  --threads INT  --chunk-size INT  --geno FLOAT  --maf FLOAT  --mac FLOAT\n"
      "\n"
      "Output: binary .bin file.  Pass path to --ind-af-coef when running SPAmix / SPAmixPlus.\n";
    return;
  }
  if (topic == "cal-pairwise-ibd") {
    std::cerr <<
      "Mode: --cal-pairwise-ibd\n"
      "  Compute pairwise IBD probabilities from a sparse GRM and PLINK genotypes.\n"
      "  The output file is required as --pairwise-ibd for SPAGRM.\n"
      "\n"
      "Required:\n"
      "  --bfile PREFIX   --out FILE\n"
      "  --sp-grm-saige FILE  |  --sp-grm-plink2 PREFIX   (exactly one)\n"
      "\n"
      "Optional:\n"
      "  --min-maf-ibd FLOAT   Min MAF for IBD calculation (default: 0.01)\n";
    return;
  }
  // Unknown topic
  std::cerr << "Unknown help topic: '" << topic << "'\n\n";
  printShortHelp();
}

Args parseArgs(int argc, char* argv[]) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    // --help [topic]
    if (arg == "--help" || arg == "-h") {
      if (i + 1 < argc && argv[i + 1][0] != '-') {
        a.helpTopic = argv[++i];
      } else {
        a.helpTopic = "__short__";
      }
      return a;  // skip further parsing
    }

    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Error: missing value for " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };

    if      (arg == "--method")                   a.method              = next();
    else if (arg == "--null-resid")               a.residFile           = next();
    else if (arg == "--covar")                    a.covarFile           = next();
    else if (arg == "--eigenvec")                 a.eigenVecsFile       = next();
    else if (arg == "--bfile")                    a.bfilePrefix         = next();
    else if (arg == "--ref-af")                   a.refAfFile           = next();
    else if (arg == "--sp-grm-saige")             a.spGrmSaigeFile      = next();
    else if (arg == "--sp-grm-plink2")            a.spGrmPlink2Prefix   = next();
    else if (arg == "--ind-af-coef")               a.indAfFile           = next();
    else if (arg == "--pairwise-ibd")             a.pairwiseIBDFile     = next();
    else if (arg == "--out")                      a.outputFile          = next();
    else if (arg == "--prevalence")               a.refPrevalence       = std::stod(next());
    else if (arg == "--batch-effect-p-threshold") a.cutoff              = std::stod(next());
    else if (arg == "--spa-z-threshold")          a.spaCutoff           = std::stod(next());
    else if (arg == "--covar-p-threshold")        a.pvalCovAdjCut       = std::stod(next());
    else if (arg == "--geno")                     a.missingCutoff       = std::stod(next());
    else if (arg == "--maf")                      a.minMafCutoff        = std::stod(next());
    else if (arg == "--mac")                      a.minMacCutoff        = std::stod(next());
    else if (arg == "--outlier-iqr-threshold")    a.outlierRatio        = std::stod(next());
    else if (arg == "--outlier-abs-bound")        a.outlierAbsBound     = std::stod(next());
    else if (arg == "--threads")                  a.nthread             = std::stoi(next());
    else if (arg == "--chunk-size")               a.nSnpPerChunk        = std::stoi(next());
    else if (arg == "--cal-ind-af-coef")          a.calIndAfCoef        = true;
    else if (arg == "--cal-pairwise-ibd")         a.calPairwiseIBD      = true;
    else if (arg == "--min-maf-ibd")              a.minMafIBD           = std::stod(next());
    else {
      std::cerr << "Error: unknown option: " << arg
                << "  (run 'grab --help' for usage)\n";
      std::exit(1);
    }
  }
  return a;
}

} // anon namespace


// ── Validation helpers ────────────────────────────────────────────────────

static void require(const std::string& val, const char* flag, const char* ctx) {
  if (val.empty()) {
    std::cerr << "Error: " << flag << " is required for " << ctx << ".\n";
    std::exit(1);
  }
}

static std::vector<std::string> splitComma(const std::string& s, const char* flag,
                                           size_t minCount) {
  std::vector<std::string> out;
  std::istringstream iss(s);
  std::string tok;
  while (std::getline(iss, tok, ','))
    if (!tok.empty()) out.push_back(tok);
  if (out.size() < minCount) {
    std::cerr << "Error: " << flag << " requires at least " << minCount
              << " comma-separated values.\n";
    std::exit(1);
  }
  return out;
}

static void checkSpGrm(const Args& a, bool required, const char* ctx) {
  if (!a.spGrmSaigeFile.empty() && !a.spGrmPlink2Prefix.empty()) {
    std::cerr << "Error: --sp-grm-saige and --sp-grm-plink2 are mutually exclusive.\n";
    std::exit(1);
  }
  if (required && a.spGrmSaigeFile.empty() && a.spGrmPlink2Prefix.empty()) {
    std::cerr << "Error: --sp-grm-saige or --sp-grm-plink2 is required for "
              << ctx << ".\n";
    std::exit(1);
  }
}


int main(int argc, char* argv[]) {
  if (argc < 2) { printShortHelp(); return 1; }

  Args args = parseArgs(argc, argv);

  // ── Help dispatch ────────────────────────────────────────────────
  if (!args.helpTopic.empty()) {
    if (args.helpTopic == "__short__")
      printShortHelp();
    else
      printTopicHelp(args.helpTopic);
    return 0;
  }

  // ── Mode: --cal-ind-af-coef ──────────────────────────────────────
  if (args.calIndAfCoef) {
    if (!args.method.empty() || args.calPairwiseIBD) {
      std::cerr << "Error: --cal-ind-af-coef cannot be combined with --method"
                   " or --cal-pairwise-ibd.\n";
      return 1;
    }
    require(args.bfilePrefix,   "--bfile",    "--cal-ind-af-coef");
    require(args.eigenVecsFile, "--eigenvec", "--cal-ind-af-coef");
    require(args.outputFile,    "--out",      "--cal-ind-af-coef");
    infoMsg("GRAB starting: --cal-ind-af-coef, nthread=%d", args.nthread);
    try {
      runSPAmixAF(
          args.eigenVecsFile, args.bfilePrefix,
          args.outputFile,
          args.nthread, args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    } catch (const std::exception& e) {
      std::cerr << "[ERROR] " << e.what() << "\n"; return 1;
    }
    return 0;
  }

  // ── Mode: --cal-pairwise-ibd ─────────────────────────────────────
  if (args.calPairwiseIBD) {
    if (!args.method.empty()) {
      std::cerr << "Error: --cal-pairwise-ibd cannot be combined with --method.\n";
      return 1;
    }
    require(args.bfilePrefix, "--bfile", "--cal-pairwise-ibd");
    require(args.outputFile,  "--out",   "--cal-pairwise-ibd");
    checkSpGrm(args, /*required=*/true, "--cal-pairwise-ibd");
    infoMsg("GRAB starting: --cal-pairwise-ibd, min-maf-ibd=%.4f", args.minMafIBD);
    try {
      runPairwiseIBD(
          args.spGrmSaigeFile, args.spGrmPlink2Prefix, args.bfilePrefix,
          args.outputFile, args.minMafIBD);
    } catch (const std::exception& e) {
      std::cerr << "[ERROR] " << e.what() << "\n"; return 1;
    }
    return 0;
  }

  // ── Mode: --method ───────────────────────────────────────────────
  static const char* const METHODS[] = {
    "SPACox","SPAGRM","SPAmix","SPAmixPlus","SPAsqr","WtCoxG","LEAF"
  };
  bool methodOk = false;
  for (const auto* m : METHODS) if (args.method == m) { methodOk = true; break; }
  if (!methodOk) {
    if (args.method.empty())
      std::cerr << "Error: --method is required (or use --cal-ind-af-coef"
                   " / --cal-pairwise-ibd).\n";
    else
      std::cerr << "Error: unknown method '" << args.method
                << "'.  Supported: SPACox SPAGRM SPAmix SPAmixPlus SPAsqr WtCoxG LEAF\n";
    return 1;
  }

  // Common required args for all GWAS methods
  require(args.bfilePrefix, "--bfile",      args.method.c_str());
  require(args.residFile,   "--null-resid", args.method.c_str());
  require(args.outputFile,  "--out",        args.method.c_str());

  if (args.nSnpPerChunk < 256) {
    std::cerr << "Error: --chunk-size must be >= 256 (got " << args.nSnpPerChunk << ")\n";
    return 1;
  }

  infoMsg("GRAB starting: method=%s, nthread=%d", args.method.c_str(), args.nthread);
  try {

    // ── SPACox ──────────────────────────────────────────────────────
    if (args.method == "SPACox") {
      runSPACox(
          args.residFile, args.covarFile, args.bfilePrefix,
          args.outputFile,
          args.pvalCovAdjCut, args.spaCutoff, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    }

    // ── SPAGRM ──────────────────────────────────────────────────────
    else if (args.method == "SPAGRM") {
      checkSpGrm(args, /*required=*/true, "SPAGRM");
      require(args.pairwiseIBDFile, "--pairwise-ibd", "SPAGRM");
      runSPAGRM(
          args.residFile, args.spGrmSaigeFile, args.spGrmPlink2Prefix,
          args.pairwiseIBDFile, args.bfilePrefix, args.outputFile,
          args.spaCutoff, args.nthread, args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    }

    // ── SPAmix ──────────────────────────────────────────────────────
    else if (args.method == "SPAmix") {
      require(args.eigenVecsFile, "--eigenvec", "SPAmix");
      runSPAmix(
          args.residFile, args.eigenVecsFile, args.bfilePrefix,
          args.indAfFile, args.outputFile,
          args.spaCutoff, args.outlierRatio, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    }

    // ── SPAmixPlus ──────────────────────────────────────────────────
    else if (args.method == "SPAmixPlus") {
      require(args.eigenVecsFile, "--eigenvec", "SPAmixPlus");
      checkSpGrm(args, /*required=*/true, "SPAmixPlus");
      runSPAmixPlus(
          args.residFile, args.eigenVecsFile, args.bfilePrefix,
          args.spGrmSaigeFile, args.spGrmPlink2Prefix, args.indAfFile,
          args.outputFile,
          args.spaCutoff, args.outlierRatio, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    }

    // ── SPAsqr ──────────────────────────────────────────────────────
    else if (args.method == "SPAsqr") {
      checkSpGrm(args, /*required=*/true, "SPAsqr");
      runSPAsqr(
          args.residFile, args.spGrmSaigeFile, args.spGrmPlink2Prefix,
          args.bfilePrefix, args.outputFile,
          args.spaCutoff, args.outlierRatio, args.outlierAbsBound,
          args.nthread, args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    }

    // ── WtCoxG ──────────────────────────────────────────────────────
    else if (args.method == "WtCoxG") {
      require(args.refAfFile, "--ref-af", "WtCoxG");
      if (args.refPrevalence <= 0.0) {
        std::cerr << "Error: --prevalence is required for WtCoxG and must be positive.\n";
        return 1;
      }
      checkSpGrm(args, /*required=*/false, "WtCoxG");
      runWtCoxG(
          args.residFile, args.bfilePrefix, args.refAfFile,
          args.spGrmSaigeFile, args.spGrmPlink2Prefix, args.outputFile,
          args.refPrevalence, args.cutoff, args.spaCutoff, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    }

    // ── LEAF ────────────────────────────────────────────────────────
    else {
      require(args.refAfFile, "--ref-af", "LEAF");
      if (args.refPrevalence <= 0.0) {
        std::cerr << "Error: --prevalence is required for LEAF and must be positive.\n";
        return 1;
      }
      checkSpGrm(args, /*required=*/false, "LEAF");

      auto residFiles = splitComma(args.residFile, "--null-resid", 2);
      auto refAfFiles = splitComma(args.refAfFile,  "--ref-af",    2);
      if (refAfFiles.size() != residFiles.size()) {
        std::cerr << "Error: LEAF requires equal numbers of --null-resid and --ref-af files"
                     " (got " << residFiles.size() << " vs " << refAfFiles.size() << ").\n";
        return 1;
      }
      runLEAF(
          residFiles, args.bfilePrefix, refAfFiles,
          args.spGrmSaigeFile, args.spGrmPlink2Prefix, args.outputFile,
          args.refPrevalence, args.cutoff, args.spaCutoff, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    }

  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << "\n";
    return 1;
  }
  return 0;
}

