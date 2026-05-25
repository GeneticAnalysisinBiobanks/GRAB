// dispatch.cpp — Argument validation and method dispatch
//
// This is the only cli translation unit that includes the heavy method
// headers (Eigen, Boost, etc.).  help.cpp and parse.cpp stay lightweight.

#include "cli/cli.hpp"
#include "cli/flags.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "util/int_pheno.hpp"
#include "util/logging.hpp"
#include "util/null_model.hpp"
#include "util/text_stream.hpp"

#include "localplus/abed_convert_msp.hpp"
#include "localplus/abed_convert_txt.hpp"
#include "localplus/spamixlocalp.hpp"
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
// Exactly one of --bfile / --pfile / --vcf / --bcf / --bgen must be given.
static GenoSpec resolveGenoSpec(
    const Args &a,
    const char *ctx
) {
    int count = !a.bfilePrefix.empty() + !a.pfilePrefix.empty() + !a.vcfFile.empty() +
                !a.bcfFile.empty() + !a.bgenFile.empty();
    if (count == 0) {
        std::cerr << "Error: a genotype input (--bfile, --pfile, --vcf, --bcf, or --bgen)"
            " is required for "
                  << ctx << ".\n";
        std::exit(1);
    }
    if (count > 1) {
        std::cerr << "Error: --bfile, --pfile, --vcf, --bcf, and --bgen are mutually"
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
    } else if (!a.bcfFile.empty()) {
        spec.format = GenoFormat::Bcf;
        spec.path = a.bcfFile;
    } else {
        spec.format = GenoFormat::Bgen;
        spec.path = a.bgenFile;
        // The CLI parser guarantees bgenRefMode ∈ {ref-first, ref-last, ref-unknown}
        // whenever --bgen is given.  Translate to the internal binary flag:
        // ref-last and ref-unknown both place ALT at alleles[0]; ref-unknown
        // additionally triggers a "REF is provisional" warning because GRAB
        // has no provisional-REF output column to propagate.
        spec.bgenAltFirst =
            (a.bgenRefMode == "ref-last" || a.bgenRefMode == "ref-unknown");
        if (a.bgenRefMode == "ref-unknown")
            warnMsg("--bgen ref-unknown: REF allele is provisional;"
                    " GRAB treats this as ref-last (alleles[0]=ALT) for the"
                    " score test and does not emit a provisional-REF marker"
                    " in the output.");
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

// Build the per-phenotype PhenoSpec list for WtCoxG / LEAF from the
// dispatch-level phenoNames vector.  Each entry is parsed via Auto mode
// (binary "COL" or survival "TIME:EVENT").  If the user supplied an
// optional --regression-model override (logistic | cox), every spec is
// checked for consistency; linear / ordinal are rejected — neither WtCoxG
// nor LEAF supports them.
static std::vector<nullmodel::PhenoSpec> buildPhenoSpecsForWtCoxGLEAF(
    const std::vector<std::string> &phenoNames,
    const std::string &regressionModelStr,
    const char *methodLabel
) {
    std::vector<nullmodel::PhenoSpec> out;
    out.reserve(phenoNames.size());
    for (const auto &tok : phenoNames)
        out.push_back(nullmodel::parsePhenoSpecAuto(tok));

    if (regressionModelStr.empty()) return out;

    auto rm = nullmodel::parseRegressionModel(regressionModelStr);
    if (rm == nullmodel::RegressionModel::Linear || rm == nullmodel::RegressionModel::Ordinal) {
        std::cerr << "Error: " << methodLabel
                  << " supports --regression-model auto | logistic | cox only"
                     " (got '" << regressionModelStr << "').\n";
        std::exit(1);
    }
    if (rm == nullmodel::RegressionModel::Auto) return out;
    const bool requireSurv = (rm == nullmodel::RegressionModel::Cox);
    for (const auto &spec : out) {
        const bool isCox = nullmodel::isCoxSpec(spec);
        if (requireSurv && !isCox) {
            std::cerr << "Error: " << methodLabel
                      << " --regression-model cox requires TIME:EVENT syntax;"
                         " got '" << spec.name << "'.\n";
            std::exit(1);
        }
        if (!requireSurv && isCox) {
            std::cerr << "Error: " << methodLabel
                      << " --regression-model logistic does not allow"
                         " TIME:EVENT syntax; got '" << spec.name << "'.\n";
            std::exit(1);
        }
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
    if (args.intPheno) std::fprintf(stderr, "  --int-pheno\n");
    if (!args.method.empty()) std::fprintf(stderr, "  --method %s\n", args.method.c_str());
    // input
    if (!args.bfilePrefix.empty()) std::fprintf(stderr, "  --bfile %s\n", args.bfilePrefix.c_str());
    if (!args.pfilePrefix.empty()) std::fprintf(stderr, "  --pfile %s\n", args.pfilePrefix.c_str());
    if (!args.vcfFile.empty()) std::fprintf(stderr, "  --vcf %s\n", args.vcfFile.c_str());
    if (!args.bcfFile.empty()) std::fprintf(stderr, "  --bcf %s\n", args.bcfFile.c_str());
    if (!args.bgenFile.empty())
        std::fprintf(stderr, "  --bgen %s %s\n",
                     args.bgenFile.c_str(), args.bgenRefMode.c_str());
    if (!args.phenoFile.empty()) std::fprintf(stderr, "  --pheno %s\n", args.phenoFile.c_str());
    if (!args.covarFile.empty()) std::fprintf(stderr, "  --covar %s\n", args.covarFile.c_str());
    if (!args.covarName.empty()) std::fprintf(stderr, "  --covar-name %s\n", args.covarName.c_str());
    if (!args.phenoName.empty()) std::fprintf(stderr, "  --pheno-name %s\n", args.phenoName.c_str());
    if (!args.residName.empty()) std::fprintf(stderr, "  --resid-name %s\n", args.residName.c_str());
    if (!args.regressionModel.empty()) std::fprintf(stderr, "  --regression-model %s\n", args.regressionModel.c_str());
    if (args.saveResid) std::fprintf(stderr, "  --save-resid\n");
    // pc-cols: relevant for SPAmix/SPAmixPlus/LEAF and cal-af-coef
    {
        bool usesPcCols =
            args.calAfCoef || args.method == "SPAmix" || args.method == "SPAmixPlus" || args.method == "LEAF";
        if (usesPcCols && !args.pcCols.empty()) std::fprintf(stderr, "  --pc-cols %s\n", args.pcCols.c_str());
    }
    if (args.method == "SAGELD" && !args.sageldX.empty())
        std::fprintf(stderr, "  --sageld-x %s\n", args.sageldX.c_str());
    // SPAsqr-specific knobs (relevant for SPAsqr pheno path)
    if (args.method == "SPAsqr" && !args.phenoName.empty()) {
        // --spasqr-mode and --spasqr-solver are always logged because the
        // two modes (score vs wald) and two solvers (qmme vs conquer)
        // have substantially different code paths, output formats, and
        // runtime behavior; surfacing them in the run header makes it
        // unambiguous which path was exercised.
        std::fprintf(stderr, "  --spasqr-mode %s\n", args.spasqrMode.c_str());
        std::fprintf(stderr, "  --spasqr-solver %s\n", args.spasqrSolver.c_str());
        if (!args.spasqrTaus.empty())
            std::fprintf(stderr, "  --spasqr-taus %s\n", args.spasqrTaus.c_str());
        if (args.spasqrTol != 1e-7)
            std::fprintf(stderr, "  --spasqr-tol %g\n", args.spasqrTol);
        if (args.spasqrH >= 0.0)
            std::fprintf(stderr, "  --spasqr-h %g\n", args.spasqrH);
        if (args.spasqrHScale >= 0.0)
            std::fprintf(stderr, "  --spasqr-h-scale %g\n", args.spasqrHScale);
    }
    // other files
    if (!args.refAfFile.empty()) std::fprintf(stderr, "  --ref-af %s\n", args.refAfFile.c_str());
    if (!args.spGrmGrabFile.empty()) std::fprintf(stderr, "  --sp-grm-grab %s\n", args.spGrmGrabFile.c_str());
    if (!args.spGrmPlink2File.empty()) std::fprintf(stderr, "  --sp-grm-plink2 %s\n", args.spGrmPlink2File.c_str());
    if (!args.indAfFile.empty()) std::fprintf(stderr, "  --ind-af-coef %s\n", args.indAfFile.c_str());
    if (!args.pairwiseIBDFile.empty()) std::fprintf(stderr, "  --pairwise-ibd %s\n", args.pairwiseIBDFile.c_str());
    if (!args.extractFile.empty()) std::fprintf(stderr, "  --extract %s\n", args.extractFile.c_str());
    if (!args.excludeFile.empty()) std::fprintf(stderr, "  --exclude %s\n", args.excludeFile.c_str());
    if (!args.chrSpec.empty()) std::fprintf(stderr, "  --chr %s\n", args.chrSpec.c_str());
    // LOCO (SPAsqr): only meaningful when --pred-list is given
    if (!args.predListFile.empty())
        std::fprintf(stderr, "  --pred-list %s\n", args.predListFile.c_str());
    if (!args.phenoTransform.empty())
        std::fprintf(stderr, "  --pheno-transform %s\n", args.phenoTransform.c_str());
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
    if (args.outlierRatio != 1.5) std::fprintf(stderr, "  --outlier-iqr-multiplier %g\n", args.outlierRatio);
    if (args.outlierAbsBound != 0.55) std::fprintf(stderr, "  --spasqr-outlier-abs-bound %g\n", args.outlierAbsBound);
    if (args.spagrmControlOutlier) std::fprintf(stderr, "  --spagrm-control-outlier\n");
    if (args.pvalCovAdjCut != 5e-5) std::fprintf(stderr, "  --covar-p-threshold %g\n", args.pvalCovAdjCut);
    if (args.cutoff != 0.1) std::fprintf(stderr, "  --batch-effect-p-threshold %g\n", args.cutoff);
    if (args.missingCutoff != 0.1) std::fprintf(stderr, "  --geno %g\n", args.missingCutoff);
    if (args.minMafCutoff != 1e-5) std::fprintf(stderr, "  --maf %g\n", args.minMafCutoff);
    if (args.minMacCutoff != 10.0) std::fprintf(stderr, "  --mac %g\n", args.minMacCutoff);
    if (args.hweCutoff != 0.0) std::fprintf(stderr, "  --hwe %g\n", args.hweCutoff);
    if (args.seed != 0) std::fprintf(stderr, "  --seed %llu\n", (unsigned long long)args.seed);
    if (args.nClusters != 0) std::fprintf(stderr, "  --leaf-nclusters %d\n", args.nClusters);
    if (!args.leafClusterFile.empty())
        std::fprintf(stderr, "  --leaf-cluster-file %s\n", args.leafClusterFile.c_str());
    if (args.leafKmeansNstart != 25)
        std::fprintf(stderr, "  --leaf-kmeans-nstart %d\n", args.leafKmeansNstart);
    if (args.compressionLevelExplicit)
        std::fprintf(stderr, "  --compression-level %d\n", args.compressionLevel);
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

    // ── Version short-circuit ──────────────────────────────────────
    // Honor --version / -V before any other validation so the call is
    // always cheap and free of side effects, even when accompanied by
    // ill-formed flags.  Stdout output (single-line "GRAB <version>")
    // is intended for grep / awk consumption in CI scripts.
    if (args.showVersion) {
        printVersion();
        return 0;
    }

    // Resolve --compression-level default per codec.  The runtime library
    // applies its built-in default (zlib 6, zstd 3) when the value is left
    // at 0; we materialize the same value into args so that "Options in
    // effect" logging, downstream messages, and any caller inspecting
    // args.compressionLevel see the explicit level that will be used.
    if (!args.compressionLevelExplicit) {
        if (args.compression == "zst") args.compressionLevel = 3;
        else if (args.compression == "gz") args.compressionLevel = 6;
    } else {
        if (args.compression == "gz" && (args.compressionLevel < 1 || args.compressionLevel > 9)) {
            std::cerr << "Error: --compression-level for gz must be in 1..9, got "
                      << args.compressionLevel << "\n";
            return 1;
        }
        if (args.compression == "zst" && (args.compressionLevel < 1 || args.compressionLevel > 22)) {
            std::cerr << "Error: --compression-level for zst must be in 1..22, got "
                      << args.compressionLevel << "\n";
            return 1;
        }
        if (args.compression.empty()) {
            std::cerr << "Error: --compression-level requires --compression gz|zst.\n";
            return 1;
        }
    }

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
        if (!args.method.empty() || args.calPairwiseIBD || args.intPheno) {
            std::cerr << "Error: --cal-af-coef cannot be combined with"
                " --method, --cal-pairwise-ibd, or --int-pheno.\n";
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
            TextWriter::assertWritable(afcOutput);
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
        if (!args.method.empty() || args.calAfCoef || args.calPairwiseIBD ||
            args.intPheno) {
            std::cerr << "Error: --cal-phi cannot be combined with"
                " --method, --cal-af-coef, --cal-pairwise-ibd, or --int-pheno.\n";
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
            TextWriter::assertWritable(phiOutput);
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
        if (!args.method.empty() || args.calAfCoef || args.calPairwiseIBD ||
            args.calPhi || args.intPheno) {
            std::cerr << "Error: --make-abed cannot be combined with"
                " --method, --cal-af-coef, --cal-pairwise-ibd,"
                " --cal-phi, or --int-pheno.\n";
            return 1;
        }
        // --vcf and --bcf both feed convertVcfMspToAbed; mutual exclusion is
        // enforced here (only one of them may be combined with --rfmix-msp at
        // a time), and content-validation happens inside the converter.
        if (!args.vcfFile.empty() && !args.bcfFile.empty()) {
            std::cerr << "Error: --vcf and --bcf are mutually exclusive.\n";
            return 1;
        }
        const bool hasVcfMsp = !args.vcfFile.empty() && !args.mspFile.empty();
        const bool hasBcfMsp = !args.bcfFile.empty() && !args.mspFile.empty();
        const bool hasGenoMsp = hasVcfMsp || hasBcfMsp;
        const bool hasTextPre = !args.admixTextPrefix.empty();
        if (hasGenoMsp && hasTextPre) {
            std::cerr << "Error: --make-abed requires either (--vcf/--bcf + --rfmix-msp) or"
                " --admix-text-prefix, not both.\n";
            return 1;
        }
        if (!hasGenoMsp && !hasTextPre) {
            std::cerr << "Error: --make-abed requires either"
                " (--vcf FILE --rfmix-msp FILE),"
                " (--bcf FILE --rfmix-msp FILE),"
                " or (--admix-text-prefix PREFIX).\n";
            return 1;
        }
        require(args.outPrefix, "--out", "--make-abed");
        logArgsInEffect(args);
        if (hasGenoMsp)infoMsg(
                "Parsing MSP file into memory; this may take several minutes and memory scales with sample count "
                "and window count."
        );
        try {
            TextWriter::assertWritable(args.outPrefix + ".abed");
            if (hasGenoMsp) {
                const std::string &genoPath = hasBcfMsp ? args.bcfFile : args.vcfFile;
                convertVcfMspToAbed(
                    genoPath,
                    /*expectBcf=*/ hasBcfMsp,
                    args.mspFile,
                    args.outPrefix,
                    args.keepFile,
                    args.removeFile,
                    args.nthread
                );
            }
            else convertTextToAbed(args.admixTextPrefix, args.outPrefix, args.keepFile, args.removeFile, args.nthread);
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] " << e.what() << "\n";
            return 1;
        }
        printTimer();
        return 0;
    }

    // ── Mode: --int-pheno ──────────────────────────────────────────
    if (args.intPheno) {
        if (!args.method.empty() || args.calAfCoef || args.calPairwiseIBD ||
            args.calPhi || args.makeAbed) {
            std::cerr << "Error: --int-pheno cannot be combined with"
                " --method or any other utility mode.\n";
            return 1;
        }
        require(args.phenoFile, "--pheno",  "--int-pheno");
        require(args.outPrefix, "--out",    "--int-pheno");
        logArgsInEffect(args);
        const std::string intOutput = args.outPrefix + ".txt";
        auto traitNames =
            args.phenoName.empty() ? std::vector<std::string>{}
                                   : splitComma(args.phenoName, "--pheno-name", 1);
        try {
            TextWriter::assertWritable(intOutput);
            runIntPheno(args.phenoFile, intOutput, traitNames);
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] " << e.what() << "\n";
            return 1;
        }
        printTimer();
        return 0;
    }

    // ── Mode: --cal-pairwise-ibd ───────────────────────────────────
    if (args.calPairwiseIBD) {
        if (!args.method.empty() || args.intPheno) {
            std::cerr << "Error: --cal-pairwise-ibd cannot be combined with"
                " --method or --int-pheno.\n";
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
            TextWriter::assertWritable(ibdOutput);
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
        else std::cerr << "Error: unknown method '" << args.method << "'.  Run 'grab2 --help' for supported methods.\n";
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
        if (args.method == "SPAsqr" || args.method == "WtCoxG" || args.method == "LEAF") {
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
        // Residual-based methods: SPACox / SPAGRM / SPAmix / SPAmixPlus /
        // SPAmixLocalPlus additionally accept --pheno-name (mutually
        // exclusive with --resid-name).  --regression-model is optional and
        // defaults to 'auto' (per-spec inference inside the engine).
        // SAGELD remains residual-only and rejects any non-auto --regression-model.
        const bool isFitCapableMethod =
            args.method == "SPACox" || args.method == "SPAGRM" ||
            args.method == "SPAmix" || args.method == "SPAmixPlus" ||
            args.method == "SPAmixLocalPlus";
        const bool hasRegressionModel = !args.regressionModel.empty();
        // Validate the regression-model value itself (rejects legacy names
        // and unknown strings) before any method-specific whitelist check.
        nullmodel::RegressionModel parsedRegressionModel = nullmodel::RegressionModel::Auto;
        if (hasRegressionModel) {
            try {
                parsedRegressionModel = nullmodel::parseRegressionModel(args.regressionModel);
            } catch (const std::exception &ex) {
                std::cerr << "Error: " << ex.what() << "\n";
                return 1;
            }
        }
        const bool regModelNonAuto = hasRegressionModel && parsedRegressionModel != nullmodel::RegressionModel::Auto;
        if (isFitCapableMethod) {
            const bool isFitPath = hasPhenoName;
            if (isFitPath && hasResidName) {
                std::cerr << "Error: " << args.method
                          << ": --pheno-name and --resid-name are"
                             " mutually exclusive.\n";
                return 1;
            }
            if (!isFitPath && !hasResidName) {
                std::cerr << "Error: " << args.method
                          << " requires either --resid-name or --pheno-name.\n";
                return 1;
            }
            if (!isFitPath && regModelNonAuto) {
                std::cerr << "Error: " << args.method
                          << ": --regression-model other than 'auto' requires --pheno-name.\n";
                return 1;
            }
            if (args.saveResid && !isFitPath) {
                std::cerr << "Error: " << args.method
                          << ": --save-resid requires --pheno-name.\n";
                return 1;
            }
            if (args.saveResid && args.outPrefix.empty()) {
                std::cerr << "Error: --save-resid requires --out.\n";
                return 1;
            }
        } else if (args.method == "SAGELD") {
            const bool hasSageldX = !args.sageldX.empty();
            if (hasPhenoName && hasResidName) {
                std::cerr << "Error: SAGELD: --pheno-name and --resid-name are mutually exclusive.\n";
                return 1;
            }
            if (!hasPhenoName && !hasResidName) {
                std::cerr << "Error: SAGELD requires either --resid-name (residual mode)"
                             " or --pheno-name + --covar-name + --sageld-x (pheno mode).\n";
                return 1;
            }
            if (hasResidName && hasSageldX) {
                std::cerr << "Error: SAGELD: --sageld-x is incompatible with --resid-name"
                             " (residual mode reads env layout from the file header).\n";
                return 1;
            }
            if (hasPhenoName && !hasSageldX) {
                std::cerr << "Error: SAGELD pheno mode requires --sageld-x (env column name).\n";
                return 1;
            }
            if (hasPhenoName && covarNames.empty()) {
                std::cerr << "Error: SAGELD pheno mode requires --covar-name"
                             " (must list every --sageld-x variable as a fixed-effect covariate).\n";
                return 1;
            }
            if (regModelNonAuto) {
                std::cerr << "Error: SAGELD does not support --regression-model (other than 'auto').\n";
                return 1;
            }
            if (args.saveResid && !hasPhenoName) {
                std::cerr << "Error: SAGELD --save-resid requires --pheno-name"
                             " (residual-input mode has no null model to save).\n";
                return 1;
            }
            if (args.saveResid && args.outPrefix.empty()) {
                std::cerr << "Error: --save-resid requires --out.\n";
                return 1;
            }
        } else if (args.method == "SPAsqr") {
            // SPAsqr only accepts auto/linear.  It does not invoke the
            // null-model engine; the regression-model value is purely
            // declarative.
            if (regModelNonAuto && parsedRegressionModel != nullmodel::RegressionModel::Linear) {
                std::cerr << "Error: SPAsqr supports --regression-model auto | linear only"
                             " (got '" << args.regressionModel << "').\n";
                return 1;
            }
        } else if (args.method == "WtCoxG" || args.method == "LEAF") {
            // WtCoxG / LEAF only accept auto / logistic / cox.
            // Detailed per-spec consistency is checked downstream in
            // buildPhenoSpecsForWtCoxGLEAF.
            if (regModelNonAuto &&
                parsedRegressionModel != nullmodel::RegressionModel::Logistic &&
                parsedRegressionModel != nullmodel::RegressionModel::Cox) {
                std::cerr << "Error: " << args.method
                          << " supports --regression-model auto | logistic | cox only"
                             " (got '" << args.regressionModel << "').\n";
                return 1;
            }
        }
    }

    if (args.nSnpPerChunk < 256) {
        std::cerr << "Error: --chunk-size must be >= 256 (got " << args.nSnpPerChunk << ")\n";
        return 1;
    }

    logArgsInEffect(args);

    // Fail fast if the output prefix's parent directory is not writable.
    // Per-phenotype output names are not known until the method runner
    // parses --pheno-name, so we probe a sentinel path that shares the
    // same parent directory.  This catches the common case of a missing
    // output directory before any computation is launched.
    try {
        TextWriter::assertWritable(args.outPrefix + ".write_test");
    } catch (const std::exception &e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }

    // Fit-path safety net: warn when the null model will be fit with an
    // intercept-only design.  The detection rule is method-specific:
    //   - SPAmix / SPAmixPlus:  null model uses ONLY --covar-name; --pc-cols
    //     feeds the per-individual AF model and is *not* added to the null-
    //     model design.  Warn whenever --covar-name is empty.
    //   - All other fit-capable methods (SPACox, SPAGRM, SAGELD, ...): the
    //     null model uses --covar-name when given, otherwise falls back to
    //     all columns of --covar.  Warn only when both are empty.
    // Intercept-only fits are statistically valid but rarely the intended
    // GWAS configuration; emit [WARN] rather than silently proceeding.
    if (!args.phenoName.empty()) {
        bool interceptOnly;
        if (args.method == "SPAmix" || args.method == "SPAmixPlus") {
            interceptOnly = covarNames.empty();
        } else {
            interceptOnly = args.covarFile.empty() && covarNames.empty();
        }
        if (interceptOnly) {
            const char *extra =
                (args.method == "SPAmix" || args.method == "SPAmixPlus")
                    ? " For SPAmix / SPAmixPlus, --pc-cols is *not* added to"
                      " the null model; list any PCs you want adjusted in"
                      " --covar-name explicitly."
                    : "";
            warnMsg(
                "Null model for '%s' will be fit with intercept only "
                "(no --covar-name covariates).%s "
                "If your data has population structure, ancestry PCs and "
                "other adjustments are strongly recommended.",
                args.phenoName.c_str(),
                extra
            );
        }
    }

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
                args.removeFile,
                args.regressionModel,
                args.phenoName,
                args.saveResid,
                args.seed
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
                args.outlierRatio,
                args.spagrmControlOutlier,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile,
                args.regressionModel,
                args.phenoName,
                effectiveCovarFile,
                covarNames,
                args.saveResid,
                args.seed
            );
        }

        // ── SAGELD ───────────────────────────────────────────────
        else if (args.method == "SAGELD") {
            checkSpGrm(args, /*required=*/ true, "SAGELD");
            require(args.pairwiseIBDFile, "--pairwise-ibd", "SAGELD");
            auto sageldXNames =
                args.sageldX.empty() ? std::vector<std::string>{} : splitComma(args.sageldX, "--sageld-x", 1);
            runSAGELD(
                args.phenoFile,
                residNames,
                phenoNames,
                covarNames,
                sageldXNames,
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
                args.saveResid,
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
            if (args.method == "SPAmixPlus") {
                checkSpGrm(args, /*required=*/ true, "SPAmixPlus");
            } else {
                // SPAmix does not consume a sparse GRM.  Reject any attempt
                // to attach one rather than silently dispatching elsewhere.
                if (!args.spGrmGrabFile.empty() || !args.spGrmPlink2File.empty()) {
                    std::cerr << "Error: SPAmix does not accept"
                                 " --sp-grm-grab or --sp-grm-plink2.\n";
                    return 1;
                }
            }
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
                args.removeFile,
                args.regressionModel,
                args.phenoName,
                covarNames,
                args.saveResid,
                args.seed
            );
        }

        // ── SPAsqr ─────────────────────────────────────────────────
        else if (args.method == "SPAsqr") {
            checkSpGrm(args, /*required=*/ false, "SPAsqr");
            if (args.spGrmGrabFile.empty() && args.spGrmPlink2File.empty()) {
                warnMsg("SPAsqr: no --sp-grm-grab/--sp-grm-plink2 provided;"
                        " falling back to identity GRM (assumes unrelated subjects)");
            }
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
            // Validate --spasqr-solver
            if (args.spasqrSolver != "qmme" && args.spasqrSolver != "conquer") {
                std::cerr << "Error: --spasqr-solver must be 'qmme' or 'conquer', got '"
                          << args.spasqrSolver << "'\n";
                return 1;
            }
            // Resolve --pheno-transform default. Default is 'int'
            // in both contexts (with and without --pred-list); LOCO PRS
            // should be trained on a Y on the same scale.
            if (args.phenoTransform.empty()) {
                args.phenoTransform = "int";
            }
            if (args.phenoTransform != "raw" &&
                args.phenoTransform != "int" &&
                args.phenoTransform != "standardize") {
                std::cerr << "Error: --pheno-transform must be 'raw', 'int', or 'standardize', got '"
                          << args.phenoTransform << "'\n";
                return 1;
            }
            if (!args.predListFile.empty() &&
                (args.phenoTransform == "raw" || args.phenoTransform == "standardize")) {
                std::cerr << "Warning: --pheno-transform " << args.phenoTransform
                          << " with --pred-list — ensure your LOCO PRS was trained"
                          << " on the same transform; otherwise scales mismatch.\n";
            }
            // Validate --spasqr-mode
            if (args.spasqrMode != "score" && args.spasqrMode != "wald") {
                std::cerr << "Error: --spasqr-mode must be 'score' or 'wald', got '"
                          << args.spasqrMode << "'\n";
                return 1;
            }
            if (args.spasqrMode == "wald") {
                runSPAsqrWald(
                    args.phenoFile,
                    effectiveCovarFile,
                    phenoNames,
                    covarNames,
                    taus,
                    geno,
                    args.predListFile,
                    args.outPrefix,
                    args.spasqrTol,
                    args.spasqrH,
                    args.spasqrHScale,
                    args.missingCutoff,
                    args.minMafCutoff,
                    args.minMacCutoff,
                    args.hweCutoff,
                    args.keepFile,
                    args.removeFile,
                    args.phenoTransform,
                    args.nthread,
                    args.nSnpPerChunk,
                    args.compression,
                    args.compressionLevel
                );
                printTimer();
                return 0;
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
                    args.removeFile,
                    args.phenoTransform,
                    args.spasqrSolver
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
                    args.removeFile,
                    args.phenoTransform,
                    args.spasqrSolver
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
            // Parse phenoNames as PhenoSpec (Auto mode); honour optional
            // --regression-model override.
            auto wtcoxgSpecs = buildPhenoSpecsForWtCoxGLEAF(
                phenoNames, args.regressionModel, "WtCoxG");
            // Multi-phenotype entry point: shared loading + parallel null models +
            // one matched-marker scan + parallel batch-effect + single engine.
            runWtCoxG(
                args.phenoFile,
                effectiveCovarFile,
                covarNames,
                wtcoxgSpecs,
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

        // ── LEAF ───────────────────────────────────────────────────
        else if (args.method == "LEAF") {
            require(args.refAfFile, "--ref-af", "LEAF");
            if (args.refPrevalence <= 0.0) {
                std::cerr << "Error: --prevalence is required for LEAF"
                    " and must be positive.\n";
                return 1;
            }
            if (pcColNames.empty() && args.leafClusterFile.empty()) {
                std::cerr << "Error: --pc-cols is required for LEAF"
                    " (unless --leaf-cluster-file is given).\n";
                return 1;
            }
            checkSpGrm(args, /*required=*/ false, "LEAF");
            auto refAfFiles = splitComma(args.refAfFile, "--ref-af", 1);
            // Cluster count: file > --leaf-nclusters > #refAfFiles.
            // When --leaf-cluster-file is given, K is inferred there; the
            // dispatcher passes a hint (0 = infer) and cross-checks inside
            // runLEAF.
            int nClusters = args.nClusters;
            if (nClusters == 0 && args.leafClusterFile.empty())
                nClusters = static_cast<int>(refAfFiles.size());
            if (args.leafClusterFile.empty() && nClusters < 2) {
                std::cerr << "Error: LEAF needs at least"
                    " 2 clusters (--leaf-nclusters or 2+ --ref-af).\n";
                return 1;
            }
            // Parse phenoNames as PhenoSpec (Auto mode); honour optional
            // --regression-model override.
            auto leafSpecs = buildPhenoSpecsForWtCoxGLEAF(
                phenoNames, args.regressionModel, "LEAF");
            // Multi-phenotype entry point: shared K-means + geno scan + summix + GRM,
            // parallel regression and batch-effect, single engine.
            runLEAF(
                args.phenoFile,
                effectiveCovarFile,
                covarNames,
                leafSpecs,
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
                args.outlierRatio,
                args.nthread,
                args.nSnpPerChunk,
                args.missingCutoff,
                args.minMafCutoff,
                args.minMacCutoff,
                args.hweCutoff,
                args.keepFile,
                args.removeFile,
                args.leafClusterFile,
                args.leafKmeansNstart
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
                args.excludeFile,
                args.covarFile,
                covarNames,
                args.regressionModel,
                args.phenoName,
                args.saveResid
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
