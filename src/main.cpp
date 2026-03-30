// main.cpp — GRAB CLI entry point
// Parses command-line arguments and delegates to the selected method workflow.

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "util/logging.hpp"
#include "spacox/spacox.hpp"
#include "spamix/spamix.hpp"
#include "spamix/spamixplus.hpp"
#include "spamix/indiv_af.hpp"
#include "wtcoxg/wtcoxg.hpp"
#include "wtcoxg/leaf.hpp"

namespace {

struct Args {
  std::string method;
  std::string residFile;
  std::string designFile;  // --design
  std::string eigenVecsFile;
  std::string bfilePrefix;
  std::string refAfFile;
  std::string sparseGrmFile;
  std::string spamixAfFile;
  std::string outputFile;
  bool   calAfCoef         = false;
  double refPrevalence     = -1.0;  // sentinel: must be set for WtCoxG, LEAF
  double cutoff            = 0.05;
  double spaCutoff         = 2.0;
  double pvalCovAdjCut     = 5e-5;
  double missingCutoff     = 0.1;
  double minMafCutoff      = 1e-4;
  double minMacCutoff      = 10.0;
  double outlierRatio      = 1.5;
  int    nthread           = 1;
  int    nSnpPerChunk      = 8192;
};

void printUsage() {
  std::cerr
    << "Usage: grab --method <SPACox|SPAmix|SPAmixPlus|WtCoxG|LEAF> [options]\n"
    << "       grab --cal-af-coef [options]\n"
    << "\n"
    << "Required:\n"
    << "  --method              NAME   Method: SPACox, SPAmix, SPAmixPlus, WtCoxG, or LEAF\n"
    << "  --bfile               PREFIX PLINK binary prefix (.bed/.bim/.fam)\n"
    << "  --null-resid          FILE   Null model residual file(s); comma-separated for LEAF\n"
    << "                               Whitespace-delimited, '#'-lines skipped, no header\n"
    << "                               Columns: ID  RESIDUAL  [WEIGHT=1]  [INDICATOR=0]\n"
    << "  --out                 FILE   Output file\n"
    << "\n"
    << "Method-specific required:\n"
    << "  --design              FILE   Design (covariate) matrix file (SPACox)\n"
    << "                               Whitespace-delimited, '#'-lines skipped\n"
    << "                               Columns: #IID  COV1  COV2  ...\n"
    << "  --eigenvec            FILE   Eigenvector (PC) file (SPAmix, SPAmixPlus, --cal-af-coef)\n"
    << "                               Whitespace-delimited, '#'-lines skipped\n"
    << "                               Columns: #IID  PC1  PC2  ...\n"
    << "  --ref-af              FILE   Reference allele-frequency file (WtCoxG, LEAF)\n"
    << "                               Whitespace-delimited; first line starting with '#' is the header\n"
    << "                               Columns: #CHROM  POS  A1  A2  A1F_POP1  N_POP1  [A1F_POP2  N_POP2  ...]\n"
    << "                               N_POPk = sample size; multiple pops supported (LEAF only)\n"
    << "  --sparse-grm          FILE   Sparse GRM file (SPAmixPlus)\n"
    << "                               Tab-delimited, '#'-lines skipped, no header\n"
    << "                               Columns: #ID1  ID2  VALUE\n"
    << "  --prevalence          FLOAT  Prevalence in reference (required for WtCoxG, LEAF)\n"
    << "\n"
    << "Optional:\n"
    << "  --af-coef             FILE   Pre-computed AF model file (SPAmix, SPAmixPlus)\n"
    << "                               Text/gz: whitespace-delimited, '#'-lines skipped\n"
    << "                               Columns: #CHROM  ID  STATUS  BETA0  BETA1  ...  (STATUS: 0=uniform 1=OLS 2=logistic)\n"
    << "                               Binary (.bin): produced by --cal-af-coef\n"
    << "  --threads             INT    Number of threads (default: 1)\n"
    << "  --chunk-size          INT    Markers per chunk (default: 8192, min: 256)\n"
    << "  --geno                FLOAT  Per-marker missing rate cutoff (default: 0.1)\n"
    << "  --maf                 FLOAT  Min minor allele frequency (default: 1e-4)\n"
    << "  --mac                 FLOAT  Min minor allele count (default: 10)\n"
    << "  --spa-z-threshold     FLOAT  SPA z-score cutoff (default: 2.0)\n"
    << "  --help                       Print this message\n";
}

Args parseArgs(int argc, char* argv[]) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    if      (arg == "--method")                    a.method             = next();
    else if (arg == "--null-resid")                a.residFile          = next();
    else if (arg == "--design")                    a.designFile         = next();
    else if (arg == "--eigenvec")                  a.eigenVecsFile      = next();
    else if (arg == "--bfile")                     a.bfilePrefix        = next();
    else if (arg == "--ref-af")                    a.refAfFile          = next();
    else if (arg == "--sparse-grm")                a.sparseGrmFile      = next();
    else if (arg == "--af-coef")                   a.spamixAfFile       = next();
    else if (arg == "--out")                       a.outputFile         = next();
    else if (arg == "--prevalence")                a.refPrevalence      = std::stod(next());
    else if (arg == "--batch-effect-p-threshold")  a.cutoff             = std::stod(next());
    else if (arg == "--spa-z-threshold")           a.spaCutoff          = std::stod(next());
    else if (arg == "--covar-p-threshold")         a.pvalCovAdjCut      = std::stod(next());
    else if (arg == "--geno")                      a.missingCutoff      = std::stod(next());
    else if (arg == "--maf")                       a.minMafCutoff       = std::stod(next());
    else if (arg == "--mac")                       a.minMacCutoff       = std::stod(next());
    else if (arg == "--resid-iqr-threshold")       a.outlierRatio       = std::stod(next());
    else if (arg == "--threads")                   a.nthread            = std::stoi(next());
    else if (arg == "--chunk-size")                a.nSnpPerChunk       = std::stoi(next());
    else if (arg == "--cal-af-coef")               a.calAfCoef          = true;
    else if (arg == "--help" || arg == "-h") { printUsage(); std::exit(0); }
    else { std::cerr << "Unknown option: " << arg << "\n"; std::exit(1); }
  }
  return a;
}

} // anon namespace


int main(int argc, char* argv[]) {
  if (argc < 2) { printUsage(); return 1; }
  Args args = parseArgs(argc, argv);

  if (args.calAfCoef && !args.method.empty()) {
    std::cerr << "Error: --cal-af-coef and --method are mutually exclusive.\n";
    return 1;
  }

  if (args.calAfCoef) {
    if (args.bfilePrefix.empty() || args.outputFile.empty()) {
      std::cerr << "Error: --bfile and --out are required.\n";
      return 1;
    }
    if (args.eigenVecsFile.empty()) {
      std::cerr << "Error: --eigenvec is required for --cal-af-coef.\n";
      return 1;
    }
    infoMsg("GRAB starting: --cal-af-coef, nthread=%d", args.nthread);
    try {
      runSPAmixAF(
          args.eigenVecsFile, args.bfilePrefix,
          args.outputFile,
          args.nthread, args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    } catch (const std::exception& e) {
      std::cerr << "[ERROR] " << e.what() << "\n";
      return 1;
    }
    return 0;
  }

  if (args.method != "SPACox" && args.method != "SPAmix" &&
      args.method != "SPAmixPlus" &&
      args.method != "WtCoxG" && args.method != "LEAF") {
    std::cerr << "Unsupported method: " << args.method
              << " (supported: SPACox, SPAmix, SPAmixPlus, WtCoxG, LEAF)\n";
    return 1;
  }
  if (args.bfilePrefix.empty() || args.outputFile.empty()) {
    std::cerr << "Error: --bfile and --out are required.\n";
    return 1;
  }
  if (args.residFile.empty()) {
    std::cerr << "Error: --null-resid is required for " << args.method << ".\n";
    return 1;
  }
  if (args.method == "SPACox" && args.designFile.empty()) {
    std::cerr << "Error: --design is required for SPACox.\n";
    return 1;
  }
  if ((args.method == "SPAmix" || args.method == "SPAmixPlus") &&
      args.eigenVecsFile.empty()) {
    std::cerr << "Error: --eigenvec is required for " << args.method << ".\n";
    return 1;
  }
  if (args.method == "SPAmixPlus" && args.sparseGrmFile.empty()) {
    std::cerr << "Error: --sparse-grm is required for SPAmixPlus.\n";
    return 1;
  }
  if ((args.method == "WtCoxG" || args.method == "LEAF") && args.refAfFile.empty()) {
    std::cerr << "Error: --ref-af is required for "
              << args.method << ".\n";
    return 1;
  }

  if ((args.method == "WtCoxG" || args.method == "LEAF") && args.refPrevalence <= 0.0) {
    std::cerr << "Error: --prevalence is required for "
              << args.method << " and must be positive.\n";
    return 1;
  }

  if (args.nSnpPerChunk < 256) {
    std::cerr << "Error: --chunk-size must be >= 256 (got "
              << args.nSnpPerChunk << ")\n";
    return 1;
  }

  infoMsg("GRAB starting: method=%s, nthread=%d", args.method.c_str(), args.nthread);
  try {
    if (args.method == "SPACox") {
      runSPACox(
          args.residFile, args.designFile, args.bfilePrefix,
          args.outputFile,
          args.pvalCovAdjCut, args.spaCutoff, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    } else if (args.method == "SPAmix") {
      runSPAmix(
          args.residFile, args.eigenVecsFile, args.bfilePrefix,
          args.spamixAfFile, args.outputFile,
          args.spaCutoff, args.outlierRatio, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    } else if (args.method == "SPAmixPlus") {
      runSPAmixPlus(
          args.residFile, args.eigenVecsFile, args.bfilePrefix,
          args.sparseGrmFile, args.spamixAfFile,
          args.outputFile,
          args.spaCutoff, args.outlierRatio, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    } else if (args.method == "WtCoxG") {
      runWtCoxG(
          args.residFile,  args.bfilePrefix, args.refAfFile,
          args.sparseGrmFile, args.outputFile,
          args.refPrevalence, args.cutoff, args.spaCutoff, args.nthread,
          args.nSnpPerChunk,
          args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
    } else {
      // LEAF: split comma-separated resid files
      std::vector<std::string> residFiles;
      {
        std::istringstream iss(args.residFile);
        std::string token;
        while (std::getline(iss, token, ','))
          if (!token.empty()) residFiles.push_back(token);
      }
      if (residFiles.size() < 2) {
        std::cerr << "Error: LEAF requires at least 2 comma-separated resid files.\n";
        return 1;
      }
      runLEAF(
          residFiles, args.bfilePrefix, args.refAfFile,
          args.sparseGrmFile, args.outputFile,
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

