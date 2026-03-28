// main.cpp — GRAB CLI entry point
// Parses command-line arguments and delegates to the selected method workflow.

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "util/logging.hpp"
#include "wtcoxg/wtcoxg.hpp"
#include "wtcoxg/leaf.hpp"

namespace {

struct Args {
  std::string method       = "WtCoxG";
  std::string residFile;
  std::string bfilePrefix;
  std::string refAfFile;
  std::string sparseGrmFile;
  std::string outputFile;
  double refPrevalence = 0.1;
  double cutoff        = 0.05;
  double spaCutoff     = 2.0;
  int    nthread       = 1;
  int    nSnpPerChunk  = 8192;
};

void printUsage() {
  std::cerr
    << "Usage: grab --method <WtCoxG|LEAF> [options]\n"
    << "  --method           NAME   Method: WtCoxG or LEAF\n"
    << "  --resid-file       FILE   Residual file(s), comma-separated for LEAF\n"
    << "  --bfile            PREFIX PLINK binary prefix (.bed/.bim/.fam)\n"
    << "  --ref-af-file      FILE   Reference allele-frequency file\n"
    << "  --sparse-grm-file  FILE   Sparse GRM file (optional)\n"
    << "  --ref-prevalence   FLOAT  Prevalence in reference (default: 0.1)\n"
    << "  --output-file      FILE   Output file\n"
    << "  --nthread          INT    Number of threads (default: 1)\n"
    << "  --nsnp-per-chunk   INT    Markers per chunk (default: 8192, min: 256)\n"
    << "  --cutoff           FLOAT  Batch-effect p-value cutoff (default: 0.05)\n"
    << "  --spa-cutoff       FLOAT  SPA z-score cutoff (default: 2.0)\n"
    << "  --help                    Print this message\n";
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
    if      (arg == "--method")           a.method        = next();
    else if (arg == "--resid-file")       a.residFile     = next();
    else if (arg == "--bfile")            a.bfilePrefix   = next();
    else if (arg == "--ref-af-file")      a.refAfFile     = next();
    else if (arg == "--sparse-grm-file")  a.sparseGrmFile = next();
    else if (arg == "--output-file")      a.outputFile    = next();
    else if (arg == "--ref-prevalence")   a.refPrevalence = std::stod(next());
    else if (arg == "--cutoff")           a.cutoff        = std::stod(next());
    else if (arg == "--spa-cutoff")       a.spaCutoff     = std::stod(next());
    else if (arg == "--nthread")          a.nthread       = std::stoi(next());
    else if (arg == "--nsnp-per-chunk")   a.nSnpPerChunk  = std::stoi(next());
    else if (arg == "--help" || arg == "-h") { printUsage(); std::exit(0); }
    else { std::cerr << "Unknown option: " << arg << "\n"; std::exit(1); }
  }
  return a;
}

} // anon namespace


int main(int argc, char* argv[]) {
  if (argc < 2) { printUsage(); return 1; }
  Args args = parseArgs(argc, argv);

  if (args.method != "WtCoxG" && args.method != "LEAF") {
    std::cerr << "Unsupported method: " << args.method
              << " (supported: WtCoxG, LEAF)\n";
    return 1;
  }
  if (args.residFile.empty() || args.bfilePrefix.empty()
      || args.refAfFile.empty() || args.outputFile.empty()) {
    std::cerr << "Error: --resid-file, --bfile, --ref-af-file, "
                 "and --output-file are required.\n";
    return 1;
  }

  if (args.nSnpPerChunk < 256) {
    std::cerr << "Error: --nsnp-per-chunk must be >= 256 (got "
              << args.nSnpPerChunk << ")\n";
    return 1;
  }

  infoMsg("GRAB starting: method=%s, nthread=%d", args.method.c_str(), args.nthread);
  try {
    if (args.method == "WtCoxG") {
      runWtCoxG(
          args.residFile,  args.bfilePrefix, args.refAfFile,
          args.sparseGrmFile, args.outputFile,
          args.refPrevalence, args.cutoff, args.spaCutoff, args.nthread,
          args.nSnpPerChunk);
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
          args.nSnpPerChunk);
    }
  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << "\n";
    return 1;
  }
  return 0;
}

