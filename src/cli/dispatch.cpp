// dispatch.cpp — Argument validation and method dispatch
//
// This is the only cli translation unit that includes the heavy method
// headers (Eigen, Boost, etc.).  help.cpp and parse.cpp stay lightweight.

#include "cli/cli.hpp"
#include "cli/flags.hpp"
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

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cli {

// ── Validation helpers ─────────────────────────────────────────────

static void require(const std::string& val, const char* flag, const char* ctx) {
    if (val.empty()) {
        std::cerr << "Error: " << flag << " is required for " << ctx << ".\n";
        std::exit(1);
    }
}

static std::vector<std::string> splitComma(const std::string& s,
                                           const char* flag,
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
    if (!a.spGrmGrabFile.empty() && !a.spGrmPlink2File.empty()) {
        std::cerr << "Error: --sp-grm-grab and --sp-grm-plink2 are mutually"
                     " exclusive.\n";
        std::exit(1);
    }
    if (required && a.spGrmGrabFile.empty() && a.spGrmPlink2File.empty()) {
        std::cerr << "Error: --sp-grm-grab or --sp-grm-plink2 is required for "
                  << ctx << ".\n";
        std::exit(1);
    }
}

// ── Entry point ────────────────────────────────────────────────────

int run(int argc, char* argv[]) {
    if (argc < 2) { printHelp("__short__"); return 1; }

    Args args = parseArgs(argc, argv);

    // ── Help dispatch ──────────────────────────────────────────────
    if (!args.helpTopic.empty()) {
        printHelp(args.helpTopic);
        return 0;
    }

    // ── Convenience: parse comma-separated column name lists ───────
    auto pcColNames   = args.pcCols.empty()
        ? std::vector<std::string>{}
        : splitComma(args.pcCols,   "--pc-cols",    1);
    auto covarNames   = args.covarName.empty()
        ? std::vector<std::string>{}
        : splitComma(args.covarName, "--covar-name", 1);

    // ── Mode: --cal-ind-af-coef ────────────────────────────────────
    if (args.calIndAfCoef) {
        if (!args.method.empty() || args.calPairwiseIBD) {
            std::cerr << "Error: --cal-ind-af-coef cannot be combined with"
                         " --method or --cal-pairwise-ibd.\n";
            return 1;
        }
        require(args.bfilePrefix, "--bfile",   "--cal-ind-af-coef");
        require(args.outputFile,  "--out",     "--cal-ind-af-coef");
        if (pcColNames.empty()) {
            std::cerr << "Error: --pc-cols is required for --cal-ind-af-coef.\n";
            return 1;
        }
        if (args.phenoFile.empty() && args.covarFile.empty()) {
            std::cerr << "Error: --pheno or --covar required for --cal-ind-af-coef"
                         " (provides PC columns).\n";
            return 1;
        }
        infoMsg("GRAB starting: --cal-ind-af-coef, nthread=%d", args.nthread);
        try {
            runSPAmixAF(
                pcColNames, args.phenoFile, args.covarFile,
                args.bfilePrefix, args.outputFile,
                args.nthread, args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] " << e.what() << "\n"; return 1;
        }
        return 0;
    }

    // ── Mode: --cal-pairwise-ibd ───────────────────────────────────
    if (args.calPairwiseIBD) {
        if (!args.method.empty()) {
            std::cerr << "Error: --cal-pairwise-ibd cannot be combined with"
                         " --method.\n";
            return 1;
        }
        require(args.bfilePrefix, "--bfile", "--cal-pairwise-ibd");
        require(args.outputFile,  "--out",   "--cal-pairwise-ibd");
        checkSpGrm(args, /*required=*/true, "--cal-pairwise-ibd");
        infoMsg("GRAB starting: --cal-pairwise-ibd, min-maf-ibd=%.4f",
                args.minMafIBD);
        try {
            runPairwiseIBD(
                args.spGrmGrabFile, args.spGrmPlink2File, args.bfilePrefix,
                args.outputFile, args.minMafIBD);
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] " << e.what() << "\n"; return 1;
        }
        return 0;
    }

    // ── Mode: --method ─────────────────────────────────────────────
    bool methodOk = false;
    for (const MethodDef* const* p = kAllMethods; *p; ++p)
        if (args.method == (*p)->name) { methodOk = true; break; }

    if (!methodOk) {
        if (args.method.empty())
            std::cerr << "Error: --method is required (or use --cal-ind-af-coef"
                         " / --cal-pairwise-ibd).\n";
        else
            std::cerr << "Error: unknown method '" << args.method
                      << "'.  Run 'grab --help' for supported methods.\n";
        return 1;
    }

    // Common required flags for all GWAS methods
    require(args.bfilePrefix, "--bfile", args.method.c_str());
    require(args.outputFile,  "--out",   args.method.c_str());

    // Methods that always need --null-resid (no --pheno alternative)
    if (args.method != "WtCoxG" && args.method != "LEAF")
        require(args.residFile, "--null-resid", args.method.c_str());

    if (args.nSnpPerChunk < 256) {
        std::cerr << "Error: --chunk-size must be >= 256 (got "
                  << args.nSnpPerChunk << ")\n";
        return 1;
    }

    infoMsg("GRAB starting: method=%s, nthread=%d",
            args.method.c_str(), args.nthread);

    try {
        // ── SPACox ─────────────────────────────────────────────────
        if (args.method == "SPACox") {
            runSPACox(
                args.residFile, covarNames,
                args.phenoFile, args.covarFile,
                args.bfilePrefix, args.outputFile,
                args.pvalCovAdjCut, args.spaCutoff, args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── SPAGRM ────────────────────────────────────────────────
        else if (args.method == "SPAGRM") {
            checkSpGrm(args, /*required=*/true, "SPAGRM");
            require(args.pairwiseIBDFile, "--pairwise-ibd", "SPAGRM");
            runSPAGRM(
                args.residFile, args.spGrmGrabFile, args.spGrmPlink2File,
                args.pairwiseIBDFile, args.bfilePrefix, args.outputFile,
                args.spaCutoff, args.nthread, args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── SPAmix ─────────────────────────────────────────────────
        else if (args.method == "SPAmix") {
            if (pcColNames.empty()) {
                std::cerr << "Error: --pc-cols is required for SPAmix.\n";
                return 1;
            }
            if (args.phenoFile.empty() && args.covarFile.empty()) {
                std::cerr << "Error: --pheno or --covar required for SPAmix"
                             " (provides PC columns).\n";
                return 1;
            }
            runSPAmix(
                args.residFile, pcColNames,
                args.phenoFile, args.covarFile,
                args.bfilePrefix, args.indAfFile, args.outputFile,
                args.spaCutoff, args.outlierRatio, args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── SPAmixPlus ─────────────────────────────────────────────
        else if (args.method == "SPAmixPlus") {
            if (pcColNames.empty()) {
                std::cerr << "Error: --pc-cols is required for SPAmixPlus.\n";
                return 1;
            }
            if (args.phenoFile.empty() && args.covarFile.empty()) {
                std::cerr << "Error: --pheno or --covar required for SPAmixPlus"
                             " (provides PC columns).\n";
                return 1;
            }
            checkSpGrm(args, /*required=*/true, "SPAmixPlus");
            runSPAmixPlus(
                args.residFile, pcColNames,
                args.phenoFile, args.covarFile,
                args.bfilePrefix,
                args.spGrmGrabFile, args.spGrmPlink2File, args.indAfFile,
                args.outputFile,
                args.spaCutoff, args.outlierRatio, args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── SPAsqr ─────────────────────────────────────────────────
        else if (args.method == "SPAsqr") {
            checkSpGrm(args, /*required=*/true, "SPAsqr");
            runSPAsqr(
                args.residFile, args.spGrmGrabFile, args.spGrmPlink2File,
                args.bfilePrefix, args.outputFile,
                args.spaCutoff, args.outlierRatio, args.outlierAbsBound,
                args.nthread, args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── WtCoxG ─────────────────────────────────────────────────
        else if (args.method == "WtCoxG") {
            require(args.refAfFile, "--ref-af", "WtCoxG");
            if (args.refPrevalence <= 0.0) {
                std::cerr << "Error: --prevalence is required for WtCoxG"
                             " and must be positive.\n";
                return 1;
            }
            if (args.residFile.empty() && args.phenoFile.empty()) {
                std::cerr << "Error: --null-resid or --pheno is required"
                             " for WtCoxG.\n";
                return 1;
            }
            if (!args.residFile.empty() && !args.phenoFile.empty()) {
                std::cerr << "Error: --null-resid and --pheno are mutually"
                             " exclusive for WtCoxG.\n";
                return 1;
            }
            checkSpGrm(args, /*required=*/false, "WtCoxG");

            if (!args.residFile.empty()) {
                // Existing --null-resid path
                runWtCoxG(
                    args.residFile, args.bfilePrefix, args.refAfFile,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            } else {
                // New --pheno path: compute regression internally
                if (args.binaryPheno.empty() && args.survPheno.empty()) {
                    std::cerr << "Error: --pheno-binary or --pheno-survival is"
                                 " required when using --pheno for WtCoxG.\n";
                    return 1;
                }
                runWtCoxGPheno(
                    args.phenoFile, args.covarFile, covarNames,
                    args.binaryPheno, args.survPheno,
                    args.bfilePrefix, args.refAfFile,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            }
        }

        // ── LEAF ───────────────────────────────────────────────────
        else {
            require(args.refAfFile, "--ref-af", "LEAF");
            if (args.refPrevalence <= 0.0) {
                std::cerr << "Error: --prevalence is required for LEAF"
                             " and must be positive.\n";
                return 1;
            }
            if (args.residFile.empty() && args.phenoFile.empty()) {
                std::cerr << "Error: --null-resid or --pheno is required"
                             " for LEAF.\n";
                return 1;
            }
            if (!args.residFile.empty() && !args.phenoFile.empty()) {
                std::cerr << "Error: --null-resid and --pheno are mutually"
                             " exclusive for LEAF.\n";
                return 1;
            }
            checkSpGrm(args, /*required=*/false, "LEAF");

            if (!args.residFile.empty()) {
                // Existing --null-resid path
                auto residFiles = splitComma(args.residFile, "--null-resid", 2);
                auto refAfFiles = splitComma(args.refAfFile,  "--ref-af",    1);
                runLEAF(
                    residFiles, args.bfilePrefix, refAfFiles,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            } else {
                // New --pheno path: K-means + per-cluster regression
                if (args.binaryPheno.empty() && args.survPheno.empty()) {
                    std::cerr << "Error: --pheno-binary or --pheno-survival is"
                                 " required when using --pheno for LEAF.\n";
                    return 1;
                }
                if (pcColNames.empty()) {
                    std::cerr << "Error: --pc-cols is required when using"
                                 " --pheno for LEAF.\n";
                    return 1;
                }
                auto refAfFiles = splitComma(args.refAfFile, "--ref-af", 1);
                int nClusters = args.nClusters > 0
                    ? args.nClusters
                    : static_cast<int>(refAfFiles.size());
                if (nClusters < 2) {
                    std::cerr << "Error: LEAF --pheno path needs at least"
                                 " 2 clusters (--leaf-nclusters or 2+ --ref-af).\n";
                    return 1;
                }
                runLEAFPheno(
                    args.phenoFile, args.covarFile, covarNames,
                    args.binaryPheno, args.survPheno, pcColNames,
                    nClusters, args.seed,
                    args.bfilePrefix, refAfFiles,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
    return 0;
}

} // namespace cli
