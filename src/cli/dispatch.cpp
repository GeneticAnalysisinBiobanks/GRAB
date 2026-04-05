// dispatch.cpp — Argument validation and method dispatch
//
// This is the only cli translation unit that includes the heavy method
// headers (Eigen, Boost, etc.).  help.cpp and parse.cpp stay lightweight.

#include "cli/cli.hpp"
#include "cli/flags.hpp"
#include "util/logging.hpp"
#include "io/geno_data.hpp"

#include "spacox/spacox.hpp"
#include "spamix/spamixplus.hpp"
#include "spamix/spamixlocalp.hpp"
#include "spamix/indiv_af.hpp"
#include "wtcoxg/wtcoxg.hpp"
#include "wtcoxg/leaf.hpp"
#include "spasqr/spasqr.hpp"
#include "polmm/polmm.hpp"
#include "spagrm/ibd.hpp"
#include "spagrm/spagrm.hpp"
#include "spagrm/sageld.hpp"
#include "io/admix_convert.hpp"
#include "io/admix_msp.hpp"

#include <algorithm>
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

// Resolve which genotype format was specified.
// Exactly one of --bfile / --pfile / --vcf / --bgen must be given.
static GenoSpec resolveGenoSpec(const Args& a, const char* ctx) {
    int count = !a.bfilePrefix.empty() + !a.pfilePrefix.empty()
              + !a.vcfFile.empty()     + !a.bgenFile.empty();
    if (count == 0) {
        std::cerr << "Error: a genotype input (--bfile, --pfile, --vcf, or --bgen)"
                     " is required for " << ctx << ".\n";
        std::exit(1);
    }
    if (count > 1) {
        std::cerr << "Error: --bfile, --pfile, --vcf, and --bgen are mutually"
                     " exclusive.\n";
        std::exit(1);
    }
    GenoSpec spec;
    spec.extractFile = a.extractFile;
    spec.excludeFile = a.excludeFile;
    if (!a.bfilePrefix.empty()) { spec.format = GenoFormat::Plink; spec.path = a.bfilePrefix; }
    else if (!a.pfilePrefix.empty()) { spec.format = GenoFormat::Pgen;  spec.path = a.pfilePrefix; }
    else if (!a.vcfFile.empty())     { spec.format = GenoFormat::Vcf;   spec.path = a.vcfFile; }
    else { spec.format = GenoFormat::Bgen; spec.path = a.bgenFile; }
    return spec;
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

// ── Parameter logging (plink2-style "Options in effect:") ──────────

static void logArgsInEffect(const Args& args) {
    infoMsg("Options in effect:");
    // modes
    if (args.calIndAfCoef)               std::fprintf(stderr, "  --cal-ind-af-coef\n");
    if (args.calPairwiseIBD)             std::fprintf(stderr, "  --cal-pairwise-ibd\n");
    if (args.calAdmixPhi)               std::fprintf(stderr, "  --cal-admix-phi\n");
    if (args.makeAbed)                   std::fprintf(stderr, "  --make-abed\n");
    if (!args.method.empty())            std::fprintf(stderr, "  --method %s\n",       args.method.c_str());
    // input
    if (!args.bfilePrefix.empty())       std::fprintf(stderr, "  --bfile %s\n",        args.bfilePrefix.c_str());
    if (!args.pfilePrefix.empty())       std::fprintf(stderr, "  --pfile %s\n",        args.pfilePrefix.c_str());
    if (!args.vcfFile.empty())           std::fprintf(stderr, "  --vcf %s\n",          args.vcfFile.c_str());
    if (!args.bgenFile.empty())          std::fprintf(stderr, "  --bgen %s\n",         args.bgenFile.c_str());
    if (!args.residFile.empty())         std::fprintf(stderr, "  --null-resid %s\n",   args.residFile.c_str());
    if (!args.phenoFile.empty())         std::fprintf(stderr, "  --pheno %s\n",        args.phenoFile.c_str());
    if (!args.covarFile.empty())         std::fprintf(stderr, "  --covar %s\n",        args.covarFile.c_str());
    if (!args.covarName.empty())         std::fprintf(stderr, "  --covar-name %s\n",   args.covarName.c_str());
    if (!args.binaryPheno.empty())       std::fprintf(stderr, "  --pheno-binary %s\n", args.binaryPheno.c_str());
    if (!args.survPheno.empty())         std::fprintf(stderr, "  --pheno-surv %s\n",   args.survPheno.c_str());
    if (!args.quantPheno.empty())        std::fprintf(stderr, "  --pheno-quant %s\n",  args.quantPheno.c_str());
    if (!args.ordinalPheno.empty())      std::fprintf(stderr, "  --pheno-ordinal %s\n", args.ordinalPheno.c_str());
    // pc-cols: relevant for SPAmix/SPAmixPlus/LEAF and cal-ind-af-coef
    {
        bool usesPcCols = args.calIndAfCoef ||
            args.method == "SPAmix" || args.method == "SPAmixPlus" || args.method == "LEAF";
        if (usesPcCols && !args.pcCols.empty())
            std::fprintf(stderr, "  --pc-cols %s\n", args.pcCols.c_str());
    }
    // spasqr-taus: relevant for SPAsqr pheno path
    if (args.method == "SPAsqr" && !args.quantPheno.empty() && !args.spasqrTaus.empty())
        std::fprintf(stderr, "  --spasqr-taus %s\n", args.spasqrTaus.c_str());
    // other files
    if (!args.refAfFile.empty())         std::fprintf(stderr, "  --ref-af %s\n",        args.refAfFile.c_str());
    if (!args.spGrmGrabFile.empty())     std::fprintf(stderr, "  --sp-grm-grab %s\n",   args.spGrmGrabFile.c_str());
    if (!args.spGrmPlink2File.empty())   std::fprintf(stderr, "  --sp-grm-plink2 %s\n", args.spGrmPlink2File.c_str());
    if (!args.indAfFile.empty())         std::fprintf(stderr, "  --ind-af-coef %s\n",   args.indAfFile.c_str());
    if (!args.pairwiseIBDFile.empty())   std::fprintf(stderr, "  --pairwise-ibd %s\n",  args.pairwiseIBDFile.c_str());
    if (!args.extractFile.empty())       std::fprintf(stderr, "  --extract %s\n",       args.extractFile.c_str());
    if (!args.excludeFile.empty())       std::fprintf(stderr, "  --exclude %s\n",       args.excludeFile.c_str());
    if (!args.admixBfilePrefix.empty())  std::fprintf(stderr, "  --admix-bfile %s\n",   args.admixBfilePrefix.c_str());
    if (!args.admixPhiFile.empty())      std::fprintf(stderr, "  --admix-phi %s\n",     args.admixPhiFile.c_str());
    if (!args.mspFile.empty())           std::fprintf(stderr, "  --rfmix-msp %s\n",     args.mspFile.c_str());
    if (!args.admixTextPrefix.empty())   std::fprintf(stderr, "  --admix-text-prefix %s\n", args.admixTextPrefix.c_str());
    if (!args.outputFile.empty())        std::fprintf(stderr, "  --out %s\n",           args.outputFile.c_str());
    if (!args.outPrefix.empty())         std::fprintf(stderr, "  --out-prefix %s\n",    args.outPrefix.c_str());
    // numeric: log only when non-default
    if (args.refPrevalence > 0.0)        std::fprintf(stderr, "  --prevalence %g\n",               args.refPrevalence);
    if (args.nthread != 1)               std::fprintf(stderr, "  --threads %d\n",                  args.nthread);
    if (args.nSnpPerChunk != 8192)       std::fprintf(stderr, "  --chunk-size %d\n",               args.nSnpPerChunk);
    if (args.spaCutoff != 2.0)           std::fprintf(stderr, "  --spa-z-threshold %g\n",          args.spaCutoff);
    if (args.outlierRatio != 1.5)        std::fprintf(stderr, "  --outlier-iqr-threshold %g\n",    args.outlierRatio);
    if (args.outlierAbsBound != 0.55)    std::fprintf(stderr, "  --outlier-abs-bound %g\n",        args.outlierAbsBound);
    if (args.pvalCovAdjCut != 5e-5)      std::fprintf(stderr, "  --covar-p-threshold %g\n",        args.pvalCovAdjCut);
    if (args.cutoff != 0.05)             std::fprintf(stderr, "  --batch-effect-p-threshold %g\n", args.cutoff);
    if (args.missingCutoff != 0.1)       std::fprintf(stderr, "  --geno %g\n",                     args.missingCutoff);
    if (args.minMafCutoff != 1e-4)       std::fprintf(stderr, "  --maf %g\n",                      args.minMafCutoff);
    if (args.minMacCutoff != 10.0)       std::fprintf(stderr, "  --mac %g\n",                      args.minMacCutoff);
    if (args.seed != 0)                  std::fprintf(stderr, "  --seed %llu\n",                   (unsigned long long)args.seed);
    if (args.nClusters != 0)             std::fprintf(stderr, "  --leaf-nclusters %d\n",            args.nClusters);
    if (args.minMafIBD != 0.01)          std::fprintf(stderr, "  --min-maf-ibd %g\n",              args.minMafIBD);
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
        require(args.outputFile,  "--out",     "--cal-ind-af-coef");
        auto geno = resolveGenoSpec(args, "--cal-ind-af-coef");
        if (pcColNames.empty()) {
            std::cerr << "Error: --pc-cols is required for --cal-ind-af-coef.\n";
            return 1;
        }
        if (args.phenoFile.empty() && args.covarFile.empty()) {
            std::cerr << "Error: --pheno or --covar required for --cal-ind-af-coef"
                         " (provides PC columns).\n";
            return 1;
        }
        logArgsInEffect(args);
        try {
            runSPAmixAF(
                pcColNames, args.phenoFile, args.covarFile,
                geno, args.outputFile,
                args.nthread, args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] " << e.what() << "\n"; return 1;
        }
        return 0;
    }

    // ── Mode: --cal-admix-phi ──────────────────────────────────────
    if (args.calAdmixPhi) {
        if (!args.method.empty() || args.calIndAfCoef || args.calPairwiseIBD) {
            std::cerr << "Error: --cal-admix-phi cannot be combined with"
                         " --method, --cal-ind-af-coef, or --cal-pairwise-ibd.\n";
            return 1;
        }
        require(args.admixBfilePrefix, "--admix-bfile", "--cal-admix-phi");
        require(args.outputFile,       "--out",         "--cal-admix-phi");
        checkSpGrm(args, /*required=*/true, "--cal-admix-phi");
        logArgsInEffect(args);
        try {
            runPhiEstimation(
                args.admixBfilePrefix,
                args.spGrmGrabFile, args.spGrmPlink2File,
                args.outputFile,
                args.extractFile, args.excludeFile);
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] " << e.what() << "\n"; return 1;
        }
        return 0;
    }

    // ── Mode: --make-abed ──────────────────────────────────────────
    if (args.makeAbed) {
        if (!args.method.empty() || args.calIndAfCoef ||
            args.calPairwiseIBD  || args.calAdmixPhi) {
            std::cerr << "Error: --make-abed cannot be combined with"
                         " --method, --cal-ind-af-coef, --cal-pairwise-ibd,"
                         " or --cal-admix-phi.\n";
            return 1;
        }
        bool hasVcfMsp  = !args.vcfFile.empty() && !args.mspFile.empty();
        bool hasTextPre = !args.admixTextPrefix.empty();
        if (hasVcfMsp && hasTextPre) {
            std::cerr << "Error: --make-abed requires either (--vcf + --rfmix-msp) or"
                         " --admix-text-prefix, not both.\n";
            return 1;
        }
        if (!hasVcfMsp && !hasTextPre) {
            std::cerr << "Error: --make-abed requires either (--vcf FILE --rfmix-msp FILE)"
                         " or (--admix-text-prefix PREFIX).\n";
            return 1;
        }
        require(args.outPrefix, "--out-prefix", "--make-abed");
        logArgsInEffect(args);
        try {
            if (hasVcfMsp)
                convertVcfMspToAbed(args.vcfFile, args.mspFile, args.outPrefix);
            else
                convertTextToAbed(args.admixTextPrefix, args.outPrefix);
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
        require(args.outputFile,  "--out",   "--cal-pairwise-ibd");
        auto geno = resolveGenoSpec(args, "--cal-pairwise-ibd");
        checkSpGrm(args, /*required=*/true, "--cal-pairwise-ibd");
        logArgsInEffect(args);
        try {
            runPairwiseIBD(
                args.spGrmGrabFile, args.spGrmPlink2File, geno,
                args.outputFile, args.minMafIBD);
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] " << e.what() << "\n"; return 1;
        }
        return 0;
    }

    // Guard: admix GWAS flags present but --method omitted
    if (args.method.empty() &&
        (!args.admixBfilePrefix.empty() || !args.admixPhiFile.empty())) {
        std::cerr << "Error: --admix-bfile and --admix-phi require"
                     " --method SPAmixLocalPlus.\n";
        return 1;
    }

    // ── Mode: --method ─────────────────────────────────────────────
    // Normalize to canonical case (accept spacox, SPACOX, SPACox, etc.)
    bool methodOk = false;
    {
        std::string lower = args.method;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        for (const MethodDef* const* p = kAllMethods; *p; ++p) {
            std::string canonLower = (*p)->name;
            std::transform(canonLower.begin(), canonLower.end(), canonLower.begin(), ::tolower);
            if (lower == canonLower) { args.method = (*p)->name; methodOk = true; break; }
        }
    }

    if (!methodOk) {
        if (args.method.empty())
            std::cerr << "Error: --method is required (or use"
                         " --cal-ind-af-coef / --cal-pairwise-ibd / --cal-admix-phi"
                         " / --make-abed).\n";
        else
            std::cerr << "Error: unknown method '" << args.method
                      << "'.  Run 'grab --help' for supported methods.\n";
        return 1;
    }

    // Common required flags for all GWAS methods
    // Exactly one of --out / --out-prefix must be provided
    if (!args.outputFile.empty() && !args.outPrefix.empty()) {
        std::cerr << "Error: --out and --out-prefix are mutually exclusive.\n";
        return 1;
    }
    if (args.outputFile.empty() && args.outPrefix.empty()) {
        std::cerr << "Error: --out or --out-prefix is required for "
                  << args.method << ".\n";
        return 1;
    }

    // --out-prefix is only supported for SPAmixLocalPlus
    bool supportsOutPrefix = (args.method == "SPAmixLocalPlus");
    if (!args.outPrefix.empty() && !supportsOutPrefix) {
        std::cerr << "Error: --out-prefix is not supported for " << args.method
                  << ". Use --out instead.\n";
        return 1;
    }

    // SPAmixLocalPlus uses --admix-bfile instead of standard geno input
    GenoSpec geno{};
    if (args.method != "SPAmixLocalPlus")
        geno = resolveGenoSpec(args, args.method.c_str());

    // Methods that always need --null-resid (no --pheno alternative)
    if (args.method != "WtCoxG" && args.method != "LEAF" &&
        args.method != "SPAsqr" && args.method != "POLMM" &&
        args.method != "SPAmixLocalPlus")
        require(args.residFile, "--null-resid", args.method.c_str());

    // Method-specific phenotype flag validation
    {
        const bool hasQuantPheno  = !args.quantPheno.empty();
        const bool hasBinaryPheno = !args.binaryPheno.empty();
        const bool hasSurvPheno   = !args.survPheno.empty();
        const bool hasOrdinalPheno = !args.ordinalPheno.empty();
        // --pheno-quant/--pheno-binary/--pheno-surv/--pheno-ordinal require --pheno (not --covar)
        if (hasQuantPheno || hasBinaryPheno || hasSurvPheno || hasOrdinalPheno) {
            if (args.phenoFile.empty()) {
                std::cerr << "Error: --pheno-quant, --pheno-binary, --pheno-surv,"
                             " and --pheno-ordinal require --pheno.\n";
                return 1;
            }
            if (!args.covarFile.empty()) {
                std::cerr << "Error: --pheno-quant, --pheno-binary, --pheno-surv,"
                             " and --pheno-ordinal cannot be combined with --covar.\n"
                             "  Use --covar-name to select covariate columns from --pheno.\n";
                return 1;
            }
        }
        if (args.method == "SPACox" || args.method == "SPAGRM" ||
            args.method == "SAGELD" ||
            args.method == "SPAmix" || args.method == "SPAmixPlus" ||
            args.method == "SPAmixLocalPlus") {
            if (hasQuantPheno || hasBinaryPheno || hasSurvPheno || hasOrdinalPheno) {
                std::cerr << "Error: " << args.method << " only accepts --null-resid."
                             " --pheno-quant, --pheno-binary, --pheno-surv, and"
                             " --pheno-ordinal are not supported for this method.\n";
                return 1;
            }
        } else if (args.method == "POLMM") {
            if (hasBinaryPheno || hasSurvPheno || hasQuantPheno) {
                std::cerr << "Error: POLMM does not support --pheno-binary,"
                             " --pheno-surv, or --pheno-quant. Use --pheno-ordinal.\n";
                return 1;
            }
            if (!hasOrdinalPheno) {
                std::cerr << "Error: POLMM requires --pheno-ordinal.\n";
                return 1;
            }
        } else if (args.method == "SPAsqr") {
            if (hasBinaryPheno || hasSurvPheno) {
                std::cerr << "Error: SPAsqr does not support --pheno-binary or"
                             " --pheno-surv. Use --pheno-quant for phenotype-based"
                             " workflow.\n";
                return 1;
            }
        } else { // WtCoxG, LEAF
            if (hasQuantPheno || hasOrdinalPheno) {
                std::cerr << "Error: " << args.method << " does not support"
                             " --pheno-quant or --pheno-ordinal."
                             " Use --pheno-binary or --pheno-surv.\n";
                return 1;
            }
        }
    }

    if (args.nSnpPerChunk < 256) {
        std::cerr << "Error: --chunk-size must be >= 256 (got "
                  << args.nSnpPerChunk << ")\n";
        return 1;
    }

    logArgsInEffect(args);

    try {
        // ── SPACox ─────────────────────────────────────────────────
        if (args.method == "SPACox") {
            runSPACox(
                args.residFile, covarNames,
                args.phenoFile, args.covarFile,
                geno, args.outputFile,
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
                args.pairwiseIBDFile, geno, args.outputFile,
                args.spaCutoff, args.nthread, args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── SAGELD ───────────────────────────────────────────────
        else if (args.method == "SAGELD") {
            checkSpGrm(args, /*required=*/true, "SAGELD");
            require(args.pairwiseIBDFile, "--pairwise-ibd", "SAGELD");
            runSAGELD(
                args.residFile, args.spGrmGrabFile, args.spGrmPlink2File,
                args.pairwiseIBDFile, geno, args.outputFile,
                args.spaCutoff, args.nthread, args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── SPAmix / SPAmixPlus (unified) ──────────────────────
        else if (args.method == "SPAmix" || args.method == "SPAmixPlus") {
            if (pcColNames.empty()) {
                std::cerr << "Error: --pc-cols is required for "
                          << args.method << ".\n";
                return 1;
            }
            if (args.phenoFile.empty() && args.covarFile.empty()) {
                std::cerr << "Error: --pheno or --covar required for "
                          << args.method
                          << " (provides PC columns).\n";
                return 1;
            }
            if (args.method == "SPAmixPlus")
                checkSpGrm(args, /*required=*/true, "SPAmixPlus");
            else
                checkSpGrm(args, /*required=*/false, "SPAmix");
            runSPAmixPlus(
                args.residFile, pcColNames,
                args.phenoFile, args.covarFile,
                geno,
                args.spGrmGrabFile, args.spGrmPlink2File, args.indAfFile,
                args.outputFile,
                args.spaCutoff, args.outlierRatio, args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── POLMM ───────────────────────────────────────────────────
        else if (args.method == "POLMM") {
            checkSpGrm(args, /*required=*/true, "POLMM");
            runPOLMM(
                args.phenoFile, args.ordinalPheno, covarNames,
                geno,
                args.spGrmGrabFile, args.spGrmPlink2File,
                args.outputFile,
                args.spaCutoff, args.nthread, args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
        }

        // ── SPAsqr ─────────────────────────────────────────────────
        else if (args.method == "SPAsqr") {
            checkSpGrm(args, /*required=*/true, "SPAsqr");
            bool hasPhenoPath = !args.quantPheno.empty();
            bool hasResidPath = !args.residFile.empty();
            if (!hasPhenoPath && !hasResidPath) {
                std::cerr << "Error: SPAsqr requires either --null-resid"
                             " or --pheno + --pheno-quant.\n";
                return 1;
            }
            if (hasPhenoPath && hasResidPath) {
                std::cerr << "Error: --null-resid and --pheno-quant"
                             " are mutually exclusive for SPAsqr.\n";
                return 1;
            }
            if (hasPhenoPath) {
                if (args.phenoFile.empty() && args.covarFile.empty()) {
                    std::cerr << "Error: --pheno or --covar required when"
                                 " using --pheno-quant.\n";
                    return 1;
                }
                if (args.spasqrTaus.empty()) {
                    std::cerr << "Error: --spasqr-taus is required when using"
                                 " --pheno-quant.\n";
                    return 1;
                }
                auto tauStrs = splitComma(args.spasqrTaus, "--spasqr-taus", 1);
                std::vector<double> taus;
                taus.reserve(tauStrs.size());
                for (const auto& s : tauStrs) {
                    double t = std::stod(s);
                    if (t <= 0.0 || t >= 1.0) {
                        std::cerr << "Error: tau values must be in (0,1), got "
                                  << s << "\n";
                        return 1;
                    }
                    taus.push_back(t);
                }
                runSPAsqrPheno(
                    args.phenoFile, args.covarFile,
                    args.quantPheno, covarNames, taus,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    geno, args.outputFile,
                    args.spaCutoff, args.outlierRatio, args.outlierAbsBound,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            } else {
                runSPAsqr(
                    args.residFile, args.spGrmGrabFile, args.spGrmPlink2File,
                    geno, args.outputFile,
                    args.spaCutoff, args.outlierRatio, args.outlierAbsBound,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            }
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
                    args.residFile, geno, args.refAfFile,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            } else {
                // New --pheno path: compute regression internally
                if (args.binaryPheno.empty() && args.survPheno.empty()) {
                    std::cerr << "Error: --pheno-binary or --pheno-surv is"
                                 " required when using --pheno for WtCoxG.\n";
                    return 1;
                }
                runWtCoxGPheno(
                    args.phenoFile, args.covarFile, covarNames,
                    args.binaryPheno, args.survPheno,
                    geno, args.refAfFile,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            }
        }

        // ── LEAF ───────────────────────────────────────────────────
        else if (args.method == "LEAF") {
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
                    residFiles, geno, refAfFiles,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            } else {
                // New --pheno path: K-means + per-cluster regression
                if (args.binaryPheno.empty() && args.survPheno.empty()) {
                    std::cerr << "Error: --pheno-binary or --pheno-surv is"
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
                    geno, refAfFiles,
                    args.spGrmGrabFile, args.spGrmPlink2File,
                    args.outputFile,
                    args.refPrevalence, args.cutoff, args.spaCutoff,
                    args.nthread, args.nSnpPerChunk,
                    args.missingCutoff, args.minMafCutoff, args.minMacCutoff);
            }
        }

        // ── SPAmixLocalPlus ────────────────────────────────────────
        else if (args.method == "SPAmixLocalPlus") {
            require(args.admixBfilePrefix, "--admix-bfile", "SPAmixLocalPlus");
            require(args.admixPhiFile,     "--admix-phi",   "SPAmixLocalPlus");
            require(args.residFile,        "--null-resid",  "SPAmixLocalPlus");
            runSPAmixLocalPlus(
                args.residFile, args.admixBfilePrefix,
                args.admixPhiFile, args.outputFile, args.outPrefix,
                args.spaCutoff, args.outlierRatio, args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff, args.minMafCutoff, args.minMacCutoff,
                args.extractFile, args.excludeFile);
        }

    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
    return 0;
}

} // namespace cli
