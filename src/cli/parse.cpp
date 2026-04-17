// parse.cpp — Command-line argument parsing

#include "cli/cli.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace cli {

// Safe numeric parsers — print a clear error and exit instead of throwing.
static double parseDouble(
    const std::string &val,
    const std::string &flag
) {
    try {
        size_t pos = 0;
        double v = std::stod(val, &pos);
        if (pos != val.size()) {
            std::cerr << "Error: " << flag << " requires a numeric value, got '" << val << "'\n";
            std::exit(1);
        }
        return v;
    } catch (const std::invalid_argument &) {
        std::cerr << "Error: " << flag << " requires a numeric value, got '" << val << "'\n";
        std::exit(1);
    } catch (const std::out_of_range &) {
        std::cerr << "Error: " << flag << " value out of range: '" << val << "'\n";
        std::exit(1);
    }
}

static int parseInt(
    const std::string &val,
    const std::string &flag
) {
    try {
        size_t pos = 0;
        int v = std::stoi(val, &pos);
        if (pos != val.size()) {
            std::cerr << "Error: " << flag << " requires an integer value, got '" << val << "'\n";
            std::exit(1);
        }
        return v;
    } catch (const std::invalid_argument &) {
        std::cerr << "Error: " << flag << " requires an integer value, got '" << val << "'\n";
        std::exit(1);
    } catch (const std::out_of_range &) {
        std::cerr << "Error: " << flag << " value out of range: '" << val << "'\n";
        std::exit(1);
    }
}

static unsigned long long parseULL(
    const std::string &val,
    const std::string &flag
) {
    try {
        size_t pos = 0;
        unsigned long long v = std::stoull(val, &pos);
        if (pos != val.size()) {
            std::cerr << "Error: " << flag << " requires an integer value, got '" << val << "'\n";
            std::exit(1);
        }
        return v;
    } catch (const std::invalid_argument &) {
        std::cerr << "Error: " << flag << " requires an integer value, got '" << val << "'\n";
        std::exit(1);
    } catch (const std::out_of_range &) {
        std::cerr << "Error: " << flag << " value out of range: '" << val << "'\n";
        std::exit(1);
    }
}

Args parseArgs(
    int argc,
    char *argv[]
) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // --help [topic]
        if (arg == "--help" || arg == "-h") {
            if (i + 1 < argc && argv[i + 1][0] != '-')a.helpTopic = argv[++i];
            else a.helpTopic = "__short__";
            return a; // skip further parsing
        }

        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for " << arg << "\n";
                std::exit(1);
            }
            return argv[++i];
        };

        if (arg == "--method")a.method = next();
        else if (arg == "--pheno")a.phenoFile = next();
        else if (arg == "--covar")a.covarFile = next();
        else if (arg == "--covar-name")a.covarName = next();
        else if (arg == "--pheno-name")a.phenoName = next();
        else if (arg == "--resid-name")a.residName = next();
        else if (arg == "--pc-cols")a.pcCols = next();
        else if (arg == "--spasqr-taus")a.spasqrTaus = next();
        else if (arg == "--spasqr-tol")a.spasqrTol = parseDouble(next(), arg);
        else if (arg == "--spasqr-h")a.spasqrH = parseDouble(next(), arg);
        else if (arg == "--spasqr-h-scale")a.spasqrHScale = parseDouble(next(), arg);
        else if (arg == "--bfile")a.bfilePrefix = next();
        else if (arg == "--pfile")a.pfilePrefix = next();
        else if (arg == "--vcf")a.vcfFile = next();
        else if (arg == "--bgen")a.bgenFile = next();
        else if (arg == "--ref-af")a.refAfFile = next();
        else if (arg == "--sp-grm-grab")a.spGrmGrabFile = next();
        else if (arg == "--sp-grm-plink2")a.spGrmPlink2File = next();
        else if (arg == "--ind-af-coef")a.indAfFile = next();
        else if (arg == "--pairwise-ibd")a.pairwiseIBDFile = next();
        else if (arg == "--out")a.outPrefix = next();
        else if (arg == "--prevalence")a.refPrevalence = parseDouble(next(), arg);
        else if (arg == "--batch-effect-p-threshold")a.cutoff = parseDouble(next(), arg);
        else if (arg == "--spa-z-threshold")a.spaCutoff = parseDouble(next(), arg);
        else if (arg == "--covar-p-threshold")a.pvalCovAdjCut = parseDouble(next(), arg);
        else if (arg == "--geno")a.missingCutoff = parseDouble(next(), arg);
        else if (arg == "--maf")a.minMafCutoff = parseDouble(next(), arg);
        else if (arg == "--mac")a.minMacCutoff = parseDouble(next(), arg);
        else if (arg == "--hwe")a.hweCutoff = parseDouble(next(), arg);
        else if (arg == "--outlier-iqr-threshold")a.outlierRatio = parseDouble(next(), arg);
        else if (arg == "--outlier-abs-bound")a.outlierAbsBound = parseDouble(next(), arg);
        else if (arg == "--threads")a.nthread = parseInt(next(), arg);
        else if (arg == "--chunk-size")a.nSnpPerChunk = parseInt(next(), arg);
        else if (arg == "--leaf-nclusters")a.nClusters = parseInt(next(), arg);
        else if (arg == "--seed")a.seed = parseULL(next(), arg);
        else if (arg == "--extract")a.extractFile = next();
        else if (arg == "--exclude")a.excludeFile = next();
        else if (arg == "--chr")a.chrSpec = next();
        else if (arg == "--keep")a.keepFile = next();
        else if (arg == "--remove")a.removeFile = next();
        else if (arg == "--pred-list")a.predListFile = next();
        else if (arg == "--admix-bfile")a.admixBfilePrefix = next();
        else if (arg == "--admix-phi")a.admixPhiFile = next();
        else if (arg == "--rfmix-msp")a.mspFile = next();
        else if (arg == "--admix-text-prefix")a.admixTextPrefix = next();
        else if (arg == "--compression")a.compression = next();
        else if (arg == "--compression-level")a.compressionLevel = parseInt(next(), arg);
        // --phi-maf-cutoff removed: hardcoded to 0.01 inside estimatePhiOneAncestry
        else if (arg == "--cal-af-coef")a.calAfCoef = true;
        else if (arg == "--cal-pairwise-ibd")a.calPairwiseIBD = true;
        else if (arg == "--cal-phi")a.calPhi = true;
        else if (arg == "--make-abed")a.makeAbed = true;
        else if (arg == "--min-maf-ibd")a.minMafIBD = parseDouble(next(), arg);
        else {
            std::cerr << "Error: unknown option: " << arg << "  (run 'grab --help' for usage)\n";
            std::exit(1);
        }
    }
    return a;
}

} // namespace cli
