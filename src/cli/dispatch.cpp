// dispatch.cpp — Argument validation and method dispatch
//
// This is the only cli translation unit that includes the heavy method
// headers (Eigen, Boost, etc.).  help.cpp and parse.cpp stay lightweight.

#include "cli/cli.hpp"
#include "cli/flags.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"

#include "localplus/abed_convert_msp.hpp"
#include "localplus/abed_convert_txt.hpp"
#include "localplus/spamixlocalp.hpp"
#include "polmm/polmm.hpp"
#include "spacox/spacox.hpp"
#include "spagrm/ibd.hpp"
#include "spagrm/sageld.hpp"
#include "spagrm/spagrm.hpp"
#include "spamix/indiv_af.hpp"
#include "spamix/spamixplus.hpp"
#include "spasqr/spasqr.hpp"
#include "wtcoxg/leaf.hpp"
#include "wtcoxg/wtcoxg.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace cli {

// ── Validation helpers ─────────────────────────────────────────────

// Parse --chr spec: "5", "2,3", "1-4,6-8,22" → set of chromosome strings.
static std::unordered_set<std::string> parseChrSpec(const std::string &spec) {
    std::unordered_set<std::string> result;
    std::istringstream iss(spec);
    std::string token;
    while (std::getline(iss, token, ',')) {
        if (token.empty()) continue;
        auto dash = token.find('-');
        if (dash != std::string::npos && dash > 0 && dash < token.size() - 1) {
            // Range: e.g. "1-4"
            int lo, hi;
            try {
                lo = std::stoi(token.substr(0, dash));
                hi = std::stoi(token.substr(dash + 1));
            } catch (const std::exception &) {
                std::cerr << "Error: invalid --chr range '" << token << "' (non-numeric).\n";
                std::exit(1);
            }
            if (lo > hi) {
                std::cerr << "Error: invalid --chr range '" << token << "' (start > end).\n";
                std::exit(1);
            }
            for (int c = lo; c <= hi; ++c)
                result.insert(std::to_string(c));
        } else {
            result.insert(token);
        }
    }
    if (result.empty()) {
        std::cerr << "Error: --chr produced an empty chromosome set.\n";
        std::exit(1);
    }
    return result;
}

static void require(
    const std::string &val,
    const char *flag,
    const char *ctx
) {
    if (val.empty()) {
        std::cerr << "Error: " << flag << " is required for " << ctx << ".\n";
        std::exit(1);
    }
}

// Resolve which genotype format was specified.
// Exactly one of --bfile / --pfile / --vcf / --bgen must be given.
static GenoSpec resolveGenoSpec(
    const Args &a,
    const char *ctx
) {
    int count = !a.bfilePrefix.empty() + !a.pfilePrefix.empty() + !a.vcfFile.empty() + !a.bgenFile.empty();
    if (count == 0) {
        std::cerr << "Error: a genotype input (--bfile, --pfile, --vcf, or --bgen)"
            " is required for "
                  << ctx << ".\n";
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
    if (!a.chrSpec.empty()) spec.chrFilter = parseChrSpec(a.chrSpec);
    if (!a.bfilePrefix.empty()) {
        spec.format = GenoFormat::Plink;
        spec.path = a.bfilePrefix;
    } else if (!a.pfilePrefix.empty()) {
        spec.format = GenoFormat::Pgen;
        spec.path = a.pfilePrefix;
    } else if (!a.vcfFile.empty()) {
        spec.format = GenoFormat::Vcf;
        spec.path = a.vcfFile;
    } else {
        spec.format = GenoFormat::Bgen;
        spec.path = a.bgenFile;
    }
    return spec;
}

static std::vector<std::string> splitComma(
    const std::string &s,
    const char *flag,
    size_t minCount
) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string tok;
    while (std::getline(iss, tok, ','))
        if (!tok.empty()) out.push_back(tok);
    if (out.size() < minCount) {
        std::cerr << "Error: " << flag << " requires at least " << minCount << " comma-separated values.\n";
        std::exit(1);
    }
    return out;
}

static void checkSpGrm(
    const Args &a,
    bool required,
    const char *ctx
) {
    if (!a.spGrmGrabFile.empty() && !a.spGrmPlink2File.empty()) {
        std::cerr << "Error: --sp-grm-grab and --sp-grm-plink2 are mutually"
            " exclusive.\n";
        std::exit(1);
    }
    if (required && a.spGrmGrabFile.empty() && a.spGrmPlink2File.empty()) {
        std::cerr << "Error: --sp-grm-grab or --sp-grm-plink2 is required for " << ctx << ".\n";
        std::exit(1);
    }
}

// ── Parameter logging (plink2-style "Options in effect:") ──────────

static void logArgsInEffect(const Args &args) {
    infoMsg("Options in effect:");
    // modes
    if (args.calAfCoef) std::fprintf(stderr, "  --cal-af-coef\n");
    if (args.calPairwiseIBD) std::fprintf(stderr, "  --cal-pairwise-ibd\n");
    if (args.calPhi) std::fprintf(stderr, "  --cal-phi\n");
    if (args.makeAbed) std::fprintf(stderr, "  --make-abed\n");
    if (!args.method.empty()) std::fprintf(stderr, "  --method %s\n", args.method.c_str());
    // input
    if (!args.bfilePrefix.empty()) std::fprintf(stderr, "  --bfile %s\n", args.bfilePrefix.c_str());
    if (!args.pfilePrefix.empty()) std::fprintf(stderr, "  --pfile %s\n", args.pfilePrefix.c_str());
    if (!args.vcfFile.empty()) std::fprintf(stderr, "  --vcf %s\n", args.vcfFile.c_str());
    if (!args.bgenFile.empty()) std::fprintf(stderr, "  --bgen %s\n", args.bgenFile.c_str());
    if (!args.phenoFile.empty()) std::fprintf(stderr, "  --pheno %s\n", args.phenoFile.c_str());
    if (!args.covarFile.empty()) std::fprintf(stderr, "  --covar %s\n", args.covarFile.c_str());
    if (!args.covarName.empty()) std::fprintf(stderr, "  --covar-name %s\n", args.covarName.c_str());
    if (!args.phenoName.empty()) std::fprintf(stderr, "  --pheno-name %s\n", args.phenoName.c_str());
    if (!args.residName.empty()) std::fprintf(stderr, "  --resid-name %s\n", args.residName.c_str());
    // pc-cols: relevant for SPAmix/SPAmixPlus/LEAF and cal-af-coef
    {
        bool usesPcCols =
            args.calAfCoef || args.method == "SPAmix" || args.method == "SPAmixPlus" || args.method == "LEAF";
        if (usesPcCols && !args.pcCols.empty()) std::fprintf(stderr, "  --pc-cols %s\n", args.pcCols.c_str());
    }
    // spasqr-taus: relevant for SPAsqr pheno path
    if (args.method == "SPAsqr" && !args.phenoName.empty() && !args.spasqrTaus.empty())std::fprintf(
            stderr,
            "  --spasqr-taus %s\n",
            args.spasqrTaus.c_str()
    );
    // other files
    if (!args.refAfFile.empty()) std::fprintf(stderr, "  --ref-af %s\n", args.refAfFile.c_str());
    if (!args.spGrmGrabFile.empty()) std::fprintf(stderr, "  --sp-grm-grab %s\n", args.spGrmGrabFile.c_str());
    if (!args.spGrmPlink2File.empty()) std::fprintf(stderr, "  --sp-grm-plink2 %s\n", args.spGrmPlink2File.c_str());
    if (!args.indAfFile.empty()) std::fprintf(stderr, "  --ind-af-coef %s\n", args.indAfFile.c_str());
    if (!args.pairwiseIBDFile.empty()) std::fprintf(stderr, "  --pairwise-ibd %s\n", args.pairwiseIBDFile.c_str());
    if (!args.extractFile.empty()) std::fprintf(stderr, "  --extract %s\n", args.extractFile.c_str());
    if (!args.excludeFile.empty()) std::fprintf(stderr, "  --exclude %s\n", args.excludeFile.c_str());
    if (!args.chrSpec.empty()) std::fprintf(stderr, "  --chr %s\n", args.chrSpec.c_str());
    if (!args.admixBfilePrefix.empty()) std::fprintf(stderr, "  --admix-bfile %s\n", args.admixBfilePrefix.c_str());
    if (!args.admixPhiFile.empty()) std::fprintf(stderr, "  --admix-phi %s\n", args.admixPhiFile.c_str());
    if (!args.mspFile.empty()) std::fprintf(stderr, "  --rfmix-msp %s\n", args.mspFile.c_str());
    if (!args.admixTextPrefix.empty()) std::fprintf(stderr, "  --admix-text-prefix %s\n", args.admixTextPrefix.c_str());
    if (!args.outPrefix.empty()) std::fprintf(stderr, "  --out %s\n", args.outPrefix.c_str());
    if (!args.compression.empty()) std::fprintf(stderr, "  --compression %s\n", args.compression.c_str());
    // numeric: log only when non-default
    if (args.refPrevalence > 0.0) std::fprintf(stderr, "  --prevalence %g\n", args.refPrevalence);
    if (args.nthread != 1) std::fprintf(stderr, "  --threads %d\n", args.nthread);
    if (args.nSnpPerChunk != 8192) std::fprintf(stderr, "  --chunk-size %d\n", args.nSnpPerChunk);
    if (args.spaCutoff != 2.0) std::fprintf(stderr, "  --spa-z-threshold %g\n", args.spaCutoff);
    if (args.outlierRatio != 1.5) std::fprintf(stderr, "  --outlier-iqr-threshold %g\n", args.outlierRatio);
    if (args.outlierAbsBound != 0.55) std::fprintf(stderr, "  --outlier-abs-bound %g\n", args.outlierAbsBound);
    if (args.pvalCovAdjCut != 5e-5) std::fprintf(stderr, "  --covar-p-threshold %g\n", args.pvalCovAdjCut);
    if (args.cutoff != 0.05) std::fprintf(stderr, "  --batch-effect-p-threshold %g\n", args.cutoff);
    if (args.missingCutoff != 0.1) std::fprintf(stderr, "  --geno %g\n", args.missingCutoff);
    if (args.minMafCutoff != 1e-5) std::fprintf(stderr, "  --maf %g\n", args.minMafCutoff);
    if (args.minMacCutoff != 10.0) std::fprintf(stderr, "  --mac %g\n", args.minMacCutoff);
    if (args.hweCutoff != 0.0) std::fprintf(stderr, "  --hwe %g\n", args.hweCutoff);
    if (args.seed != 0) std::fprintf(stderr, "  --seed %llu\n", (unsigned long long)args.seed);
    if (args.nClusters != 0) std::fprintf(stderr, "  --leaf-nclusters %d\n", args.nClusters);
    if (args.compressionLevel != 0) std::fprintf(stderr, "  --compression-level %d\n", args.compressionLevel);
    if (args.minMafIBD != 0.01) std::fprintf(stderr, "  --min-maf-ibd %g\n", args.minMafIBD);
}

// ── Entry point ────────────────────────────────────────────────────

int run(
    int argc,
    char *argv[]
) {
    if (argc < 2) {
        printHelp("__short__");
        return 1;
    }

    const auto wallStart = std::chrono::steady_clock::now();
    const std::clock_t cpuStart = std::clock();

    Args args = parseArgs(argc, argv);

    auto printTimer = [&]() {
        const double wallSec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() -
            wallStart
        ).count();
        const double cpuSec = static_cast<double>(std::clock() - cpuStart) / CLOCKS_PER_SEC;
        infoMsg("Wall time: %.1f seconds, CPU time: %.1f seconds", wallSec, cpuSec);
    };

    // ── Help dispatch ──────────────────────────────────────────────
    if (!args.helpTopic.empty()) {
        printHelp(args.helpTopic);
        return 0;
    }

    // ── Convenience: parse comma-separated column name lists ───────
    auto pcColNames = args.pcCols.empty() ? std::vector<std::string>{} : splitComma(args.pcCols, "--pc-cols", 1);
    auto covarNames =
        args.covarName.empty() ? std::vector<std::string>{} : splitComma(args.covarName, "--covar-name", 1);
    auto phenoNames =
        args.phenoName.empty() ? std::vector<std::string>{} : splitComma(args.phenoName, "--pheno-name", 1);
    auto residNames =
        args.residName.empty() ? std::vector<std::string>{} : splitComma(args.residName, "--resid-name", 1);

    // Covar-from-pheno fallback: when --covar is absent, --covar-name
    // selects columns from --pheno instead.
    std::string effectiveCovarFile;
    if (!args.covarFile.empty())effectiveCovarFile = args.covarFile;
    else if (!covarNames.empty())effectiveCovarFile = args.phenoFile;

    // ── Mode: --cal-af-coef ────────────────────────────────────────────
    if (args.calAfCoef) {
        if (!args.method.empty() || args.calPairwiseIBD) {
            std::cerr << "Error: --cal-af-coef cannot be combined with"
                " --method or --cal-pairwise-ibd.\n";
            return 1;
        }
        require(args.outPrefix, "--out", "--cal-af-coef");
        auto geno = resolveGenoSpec(args, "--cal-af-coef");
        if (pcColNames.empty()) {
            std::cerr << "Error: --pc-cols is required for --cal-af-coef.\n";
            return 1;
        }
        if (args.phenoFile.empty() && args.covarFile.empty()) {
            std::cerr << "Error: --pheno or --covar required for --cal-af-coef"
                " (provides PC columns).\n";
            return 1;
        }
        logArgsInEffect(args);
        std::string afcSuffix = ".afc";
        if (args.compression == "gz")afcSuffix += ".gz";
        else if (args.compression == "zst")afcSuffix += ".zst";
        std::string afcOutput = args.outPrefix + afcSuffix;
        try {
            runSPAmixAF(
                pcColNames,
                args.phenoFile,
                args.covarFile,
                geno,
                afcOutput,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] " << e.what() << "\n";
            return 1;
        }
        printTimer();
        return 0;
    }

    // ── Mode: --cal-phi ─────────────────────────────────────────────
    if (args.calPhi) {
        if (!args.method.empty() || args.calAfCoef || args.calPairwiseIBD) {
            std::cerr << "Error: --cal-phi cannot be combined with"
                " --method, --cal-af-coef, or --cal-pairwise-ibd.\n";
            return 1;
        }
        require(args.admixBfilePrefix, "--admix-bfile", "--cal-phi");
        require(args.outPrefix, "--out", "--cal-phi");
        checkSpGrm(args, /*required=*/ true, "--cal-phi");
        logArgsInEffect(args);
        std::string phiSuffix = ".phi";
        if (args.compression == "gz")phiSuffix += ".gz";
        else if (args.compression == "zst")phiSuffix += ".zst";
        std::string phiOutput = args.outPrefix + phiSuffix;
        try {
            runPhiEstimation(
                args.admixBfilePrefix,
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                phiOutput,
                args.keepFile,
                args.removeFile,
                args.extractFile,
                args.excludeFile,
                args.nthread
            );
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] " << e.what() << "\n";
            return 1;
        }
        printTimer();
        return 0;
    }

    // ── Mode: --make-abed ──────────────────────────────────────────
    if (args.makeAbed) {
        if (!args.method.empty() || args.calAfCoef || args.calPairwiseIBD || args.calPhi) {
            std::cerr << "Error: --make-abed cannot be combined with"
                " --method, --cal-af-coef, --cal-pairwise-ibd,"
                " or --cal-phi.\n";
            return 1;
        }
        bool hasVcfMsp = !args.vcfFile.empty() && !args.mspFile.empty();
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
        require(args.outPrefix, "--out", "--make-abed");
        logArgsInEffect(args);
        if (hasVcfMsp)infoMsg(
                "Parsing MSP file into memory; this may take several minutes and memory scales with sample count "
                "and window count."
        );
        try {
            if (hasVcfMsp)convertVcfMspToAbed(
                    args.vcfFile,
                    args.mspFile,
                    args.outPrefix,
                    args.keepFile,
                    args.removeFile,
                    args.nthread
            );
            else convertTextToAbed(args.admixTextPrefix, args.outPrefix, args.keepFile, args.removeFile, args.nthread);
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] " << e.what() << "\n";
            return 1;
        }
        printTimer();
        return 0;
    }

    // ── Mode: --cal-pairwise-ibd ───────────────────────────────────
    if (args.calPairwiseIBD) {
        if (!args.method.empty()) {
            std::cerr << "Error: --cal-pairwise-ibd cannot be combined with"
                " --method.\n";
            return 1;
        }
        require(args.outPrefix, "--out", "--cal-pairwise-ibd");
        auto geno = resolveGenoSpec(args, "--cal-pairwise-ibd");
        checkSpGrm(args, /*required=*/ true, "--cal-pairwise-ibd");
        logArgsInEffect(args);
        std::string ibdSuffix = ".ibd";
        if (args.compression == "gz")ibdSuffix += ".gz";
        else if (args.compression == "zst")ibdSuffix += ".zst";
        std::string ibdOutput = args.outPrefix + ibdSuffix;
        try {
            runPairwiseIBD(
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                geno,
                ibdOutput,
                args.keepFile,
                args.removeFile,
                args.minMafIBD,
                args.nthread
            );
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] " << e.what() << "\n";
            return 1;
        }
        printTimer();
        return 0;
    }

    // Guard: admix GWAS flags present but --method omitted
    if (args.method.empty() && (!args.admixBfilePrefix.empty() || !args.admixPhiFile.empty())) {
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
        for (const MethodDef *const *p = kAllMethods; *p; ++p) {
            std::string canonLower = (*p)->name;
            std::transform(canonLower.begin(), canonLower.end(), canonLower.begin(), ::tolower);
            if (lower == canonLower) {
                args.method = (*p)->name;
                methodOk = true;
                break;
            }
        }
    }

    if (!methodOk) {
        if (args.method.empty())
            std::cerr << "Error: --method is required (or use"
                " --cal-af-coef / --cal-pairwise-ibd / --cal-phi"
                " / --make-abed).\n";
        else std::cerr << "Error: unknown method '" << args.method << "'.  Run 'grab --help' for supported methods.\n";
        return 1;
    }

    // Common required flags for all GWAS methods
    require(args.outPrefix, "--out", args.method.c_str());

    // Validate --compression
    if (!args.compression.empty() && args.compression != "gz" && args.compression != "zst") {
        std::cerr << "Error: --compression must be 'gz' or 'zst', got '" << args.compression << "'\n";
        return 1;
    }

    // SPAmixLocalPlus uses --admix-bfile instead of standard geno input
    GenoSpec geno{};
    if (args.method != "SPAmixLocalPlus") geno = resolveGenoSpec(args, args.method.c_str());

    // All GWAS methods require --pheno (except SPAmixLocalPlus which uses --admix-bfile)
    if (args.method != "SPAmixLocalPlus")require(args.phenoFile, "--pheno", args.method.c_str());

    // Method-specific phenotype flag validation
    {
        const bool hasPhenoName = !args.phenoName.empty();
        const bool hasResidName = !args.residName.empty();
        if (args.method == "SPAsqr" || args.method == "POLMM" ||
            args.method == "WtCoxG" || args.method == "LEAF") {
            // Auto-populate from pheno file header when --pheno-name is absent
            if (!hasPhenoName && !args.phenoFile.empty()) {
                phenoNames = SubjectData::readColumnNames(args.phenoFile);
                if (args.method == "WtCoxG" || args.method == "LEAF") {
                    infoMsg("--pheno-name not specified; using all %zu phenotype"
                            " columns as binary traits", phenoNames.size());
                } else {
                    infoMsg("--pheno-name not specified; using all %zu phenotype columns",
                            phenoNames.size());
                }
            }
            if (phenoNames.empty()) {
                std::cerr << "Error: " << args.method << " requires --pheno-name"
                    " (or a --pheno file with named columns).\n";
                return 1;
            }
        }
        // Residual-based methods require --resid-name
        if (args.method == "SPACox" || args.method == "SPAGRM" ||
            args.method == "SAGELD" || args.method == "SPAmix" ||
            args.method == "SPAmixPlus" || args.method == "SPAmixLocalPlus") {
            if (!hasResidName) {
                std::cerr << "Error: " << args.method << " requires --resid-name.\n";
                return 1;
            }
        }
    }

    if (args.nSnpPerChunk < 256) {
        std::cerr << "Error: --chunk-size must be >= 256 (got " << args.nSnpPerChunk << ")\n";
        return 1;
    }

    logArgsInEffect(args);

    // Build suffix-based output path for single-file methods (SAGELD, POLMM, SPAsqr, WtCoxG, LEAF)
    auto buildOutputPath = [&](const std::string &suffix) -> std::string {
        std::string path = args.outPrefix + suffix;
        if (args.compression == "gz") path += ".gz";
        else if (args.compression == "zst") path += ".zst";
        return path;
    };

    try {
        // ── SPACox ─────────────────────────────────────────────────
        if (args.method == "SPACox") {
            runSPACox(
                residNames,
                covarNames,
                args.phenoFile,
                effectiveCovarFile,
                geno,
                args.outPrefix,
                args.compression,
                args.compressionLevel,
                args.pvalCovAdjCut,
                args.spaCutoff,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        }

        // ── SPAGRM ────────────────────────────────────────────────
        else if (args.method == "SPAGRM") {
            checkSpGrm(args, /*required=*/ true, "SPAGRM");
            require(args.pairwiseIBDFile, "--pairwise-ibd", "SPAGRM");
            runSPAGRM(
                args.phenoFile,
                residNames,
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                args.pairwiseIBDFile,
                geno,
                args.outPrefix,
                args.compression,
                args.compressionLevel,
                args.spaCutoff,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        }

        // ── SAGELD ───────────────────────────────────────────────
        else if (args.method == "SAGELD") {
            checkSpGrm(args, /*required=*/ true, "SAGELD");
            require(args.pairwiseIBDFile, "--pairwise-ibd", "SAGELD");
            runSAGELD(
                args.phenoFile,
                residNames,
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                args.pairwiseIBDFile,
                geno,
                buildOutputPath(".SAGELD"),
                args.spaCutoff,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        }

        // ── SPAmix / SPAmixPlus (unified) ──────────────────────
        else if (args.method == "SPAmix" || args.method == "SPAmixPlus") {
            if (pcColNames.empty()) {
                std::cerr << "Error: --pc-cols is required for " << args.method << ".\n";
                return 1;
            }
            if (args.phenoFile.empty() && args.covarFile.empty()) {
                std::cerr << "Error: --pheno or --covar required for " << args.method << " (provides PC columns).\n";
                return 1;
            }
            if (args.method == "SPAmixPlus")checkSpGrm(args, /*required=*/ true, "SPAmixPlus");
            else checkSpGrm(args, /*required=*/ false, "SPAmix");
            runSPAmixPlus(
                residNames,
                pcColNames,
                args.phenoFile,
                args.covarFile,
                geno,
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                args.indAfFile,
                args.outPrefix,
                args.compression,
                args.compressionLevel,
                args.spaCutoff,
                args.outlierRatio,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        }

        // ── POLMM ───────────────────────────────────────────────────
        else if (args.method == "POLMM") {
            checkSpGrm(args, /*required=*/ true, "POLMM");
            runPOLMM(
                args.phenoFile,
                effectiveCovarFile,
                phenoNames,
                covarNames,
                geno,
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                args.outPrefix,
                args.compression,
                args.compressionLevel,
                args.spaCutoff,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        }

        // ── SPAsqr ─────────────────────────────────────────────────
        else if (args.method == "SPAsqr") {
            checkSpGrm(args, /*required=*/ true, "SPAsqr");
            if (args.phenoFile.empty() && args.covarFile.empty()) {
                std::cerr << "Error: --pheno or --covar required for SPAsqr.\n";
                return 1;
            }
            if (args.spasqrH >= 0.0 && args.spasqrHScale >= 0.0) {
                std::cerr << "Error: --spasqr-h and --spasqr-h-scale"
                    " are mutually exclusive.\n";
                return 1;
            }
            if (args.spasqrTaus.empty()) {
                std::cerr << "Error: --spasqr-taus is required for SPAsqr.\n";
                return 1;
            }
            auto tauStrs = splitComma(args.spasqrTaus, "--spasqr-taus", 1);
            if (tauStrs.size() > 20) {
                std::cerr << "Error: --spasqr-taus accepts at most 20 tau levels, got "
                          << tauStrs.size() << "\n";
                return 1;
            }
            std::vector<double> taus;
            taus.reserve(tauStrs.size());
            for (const auto &s : tauStrs) {
                double t = std::stod(s);
                if (t <= 0.0 || t >= 1.0) {
                    std::cerr << "Error: tau values must be in (0,1), got " << s << "\n";
                    return 1;
                }
                taus.push_back(t);
            }
            if (!args.predListFile.empty()) {
                // Early validation: check all phenotypes have LOCO entries
                {
                    std::ifstream pf(args.predListFile);
                    if (!pf.is_open()) {
                        std::cerr << "Error: cannot open --pred-list file: "
                                  << args.predListFile << "\n";
                        return 1;
                    }
                    std::unordered_set<std::string> predPhenos;
                    std::string ln;
                    while (std::getline(pf, ln)) {
                        if (ln.empty()) continue;
                        std::istringstream lss(ln);
                        std::string ph;
                        lss >> ph;
                        if (!ph.empty()) predPhenos.insert(ph);
                    }
                    for (const auto &pn : phenoNames) {
                        if (predPhenos.find(pn) == predPhenos.end()) {
                            std::cerr << "Error: phenotype '" << pn
                                      << "' not found in --pred-list file: "
                                      << args.predListFile << "\n";
                            return 1;
                        }
                    }
                }
                runSPAsqrLoco(
                    args.phenoFile,
                    effectiveCovarFile,
                    phenoNames,
                    covarNames,
                    taus,
                    args.spGrmGrabFile,
                    args.spGrmPlink2File,
                    geno,
                    args.predListFile,
                    args.outPrefix,
                    args.compression,
                    args.compressionLevel,
                    args.spaCutoff,
                    args.outlierRatio,
                    args.outlierAbsBound,
                    args.nthread,
                    args.nSnpPerChunk,
                    args.missingCutoff,
                    args.minMafCutoff,
                    args.minMacCutoff,
                    args.hweCutoff,
                    args.spasqrTol,
                    args.spasqrH,
                    args.spasqrHScale,
                    args.keepFile,
                    args.removeFile
                );
            } else {
                runSPAsqr(
                    args.phenoFile,
                    effectiveCovarFile,
                    phenoNames,
                    covarNames,
                    taus,
                    args.spGrmGrabFile,
                    args.spGrmPlink2File,
                    geno,
                    args.outPrefix,
                    args.compression,
                    args.compressionLevel,
                    args.spaCutoff,
                    args.outlierRatio,
                    args.outlierAbsBound,
                    args.nthread,
                    args.nSnpPerChunk,
                    args.missingCutoff,
                    args.minMafCutoff,
                    args.minMacCutoff,
                    args.hweCutoff,
                    args.spasqrTol,
                    args.spasqrH,
                    args.spasqrHScale,
                    args.keepFile,
                    args.removeFile
                );
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
            checkSpGrm(args, /*required=*/ false, "WtCoxG");
            // Multi-phenotype entry point: shared loading + parallel null models +
            // one matched-marker scan + parallel batch-effect + single engine.
            runWtCoxG(
                args.phenoFile,
                effectiveCovarFile,
                covarNames,
                phenoNames,
                geno,
                args.refAfFile,
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                args.outPrefix,
                args.compression,
                args.compressionLevel,
                args.refPrevalence,
                args.cutoff,
                args.spaCutoff,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        }

        // ── LEAF ───────────────────────────────────────────────────
        else if (args.method == "LEAF") {
            require(args.refAfFile, "--ref-af", "LEAF");
            if (args.refPrevalence <= 0.0) {
                std::cerr << "Error: --prevalence is required for LEAF"
                    " and must be positive.\n";
                return 1;
            }
            if (pcColNames.empty()) {
                std::cerr << "Error: --pc-cols is required for LEAF.\n";
                return 1;
            }
            checkSpGrm(args, /*required=*/ false, "LEAF");
            auto refAfFiles = splitComma(args.refAfFile, "--ref-af", 1);
            int nClusters = args.nClusters > 0 ? args.nClusters : static_cast<int>(refAfFiles.size());
            if (nClusters < 2) {
                std::cerr << "Error: LEAF needs at least"
                    " 2 clusters (--leaf-nclusters or 2+ --ref-af).\n";
                return 1;
            }
            // Multi-phenotype entry point: shared K-means + geno scan + summix + GRM,
            // parallel regression and batch-effect, single engine.
            runLEAF(
                args.phenoFile,
                effectiveCovarFile,
                covarNames,
                phenoNames,
                pcColNames,
                nClusters,
                args.seed,
                geno,
                refAfFiles,
                args.spGrmGrabFile,
                args.spGrmPlink2File,
                args.outPrefix,
                args.compression,
                args.compressionLevel,
                args.refPrevalence,
                args.cutoff,
                args.spaCutoff,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile
            );
        }

        // ── SPAmixLocalPlus ────────────────────────────────────────
        else if (args.method == "SPAmixLocalPlus") {
            require(args.admixBfilePrefix, "--admix-bfile", "SPAmixLocalPlus");
            require(args.admixPhiFile, "--admix-phi", "SPAmixLocalPlus");
            require(args.phenoFile, "--pheno", "SPAmixLocalPlus");
            runSPAmixLocalPlus(
                args.phenoFile,
                residNames,
                args.admixBfilePrefix,
                args.admixPhiFile,
                args.outPrefix,
                args.compression,
                args.compressionLevel,
                args.spaCutoff,
                args.outlierRatio,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile,
                args.extractFile,
                args.excludeFile
            );
        }

    } catch (const std::exception &e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
    printTimer();
    return 0;
}

} // namespace cli
