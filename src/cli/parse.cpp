// parse.cpp — Command-line argument parsing

#include "cli/cli.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace cli {

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
        else if (arg == "--spasqr-tol")a.spasqrTol = std::stod(next());
        else if (arg == "--spasqr-h")a.spasqrH = std::stod(next());
        else if (arg == "--spasqr-h-scale")a.spasqrHScale = std::stod(next());
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
        else if (arg == "--prevalence")a.refPrevalence = std::stod(next());
        else if (arg == "--batch-effect-p-threshold")a.cutoff = std::stod(next());
        else if (arg == "--spa-z-threshold")a.spaCutoff = std::stod(next());
        else if (arg == "--covar-p-threshold")a.pvalCovAdjCut = std::stod(next());
        else if (arg == "--geno")a.missingCutoff = std::stod(next());
        else if (arg == "--maf")a.minMafCutoff = std::stod(next());
        else if (arg == "--mac")a.minMacCutoff = std::stod(next());
        else if (arg == "--hwe")a.hweCutoff = std::stod(next());
        else if (arg == "--outlier-iqr-threshold")a.outlierRatio = std::stod(next());
        else if (arg == "--outlier-abs-bound")a.outlierAbsBound = std::stod(next());
        else if (arg == "--threads")a.nthread = std::stoi(next());
        else if (arg == "--chunk-size")a.nSnpPerChunk = std::stoi(next());
        else if (arg == "--leaf-nclusters")a.nClusters = std::stoi(next());
        else if (arg == "--seed")a.seed = std::stoull(next());
        else if (arg == "--extract")a.extractFile = next();
        else if (arg == "--exclude")a.excludeFile = next();
        else if (arg == "--chr")a.chrSpec = next();
        else if (arg == "--keep")a.keepFile = next();
        else if (arg == "--remove")a.removeFile = next();
        else if (arg == "--admix-bfile")a.admixBfilePrefix = next();
        else if (arg == "--admix-phi")a.admixPhiFile = next();
        else if (arg == "--rfmix-msp")a.mspFile = next();
        else if (arg == "--admix-text-prefix")a.admixTextPrefix = next();
        else if (arg == "--compression")a.compression = next();
        else if (arg == "--compression-level")a.compressionLevel = std::stoi(next());
        // --phi-maf-cutoff removed: hardcoded to 0.01 inside estimatePhiOneAncestry
        else if (arg == "--cal-af-coef")a.calAfCoef = true;
        else if (arg == "--cal-pairwise-ibd")a.calPairwiseIBD = true;
        else if (arg == "--cal-phi")a.calPhi = true;
        else if (arg == "--make-abed")a.makeAbed = true;
        else if (arg == "--min-maf-ibd")a.minMafIBD = std::stod(next());
        else {
            std::cerr << "Error: unknown option: " << arg << "  (run 'grab --help' for usage)\n";
            std::exit(1);
        }
    }
    return a;
}

} // namespace cli
