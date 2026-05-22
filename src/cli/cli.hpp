// cli.hpp — Public interface for the GRAB command-line layer
#pragma once

#include <cstdint>
#include <string>

namespace cli {

struct Args {
    std::string method;
    std::string helpTopic; // set when --help <topic> is used
    std::string phenoFile;
    std::string covarFile;
    std::string covarName;                          // comma-separated covariate column names
    std::string phenoName;                          // comma-separated phenotype column names
    std::string residName;                          // comma-separated residual column names
    std::string regressionModel;                    // --regression-model: auto|linear|logistic|cox|ordinal
    bool saveResid = false;                         // --save-resid: write fitted residuals to PREFIX.null.resid

    std::string pcCols = "PC1,PC2,PC3,PC4";         // comma-separated PC column names (default: PC1,PC2,PC3,PC4)
    std::string spasqrTaus = "0.1,0.3,0.5,0.7,0.9"; // default tau levels (SPAsqr)
    std::string sageldX;                            // --sageld-x: comma-separated env names for SAGELD pheno mode
    double spasqrTol = 1e-7;                          // --spasqr-tol (conquer convergence tolerance)
    double spasqrH = -1.0;                            // --spasqr-h (explicit bandwidth; -1 = auto)
    double spasqrHScale = -1.0;                       // --spasqr-h-scale (IQR divisor; -1 = auto → 3)
    std::string bfilePrefix;
    std::string pfilePrefix; // --pfile (pgen/pvar/psam)
    std::string vcfFile;     // --vcf (VCF text, or BGZF-compressed .vcf.gz)
    std::string bcfFile;     // --bcf (BCF2 binary).  Mutually exclusive with --vcf;
                             // both flags route to the same htslib backend, but
                             // content is verified against the chosen flag to
                             // match plink2's --vcf / --bcf separation.
    std::string bgenFile;    // --bgen <filename>
    // --bgen <REF/ALT mode>: ref-first | ref-last | ref-unknown
    // (matches plink2 --bgen syntax; mandatory whenever --bgen is given).
    // ref-first   → alleles[0] is REF (IMPUTE / UK Biobank convention)
    // ref-last    → alleles[0] is ALT (plink2 default --export bgen-1.x)
    // ref-unknown → REF status unknown; alleles are kept in BGEN order with a
    //               "REF is provisional" warning (treated identically to
    //               ref-last for downstream computation).
    std::string bgenRefMode;
    std::string refAfFile;
    std::string spGrmGrabFile;   // --sp-grm-grab
    std::string spGrmPlink2File; // --sp-grm-plink2 (.grm.sp file)
    std::string indAfFile;       // --ind-af-coef
    std::string pairwiseIBDFile;
    std::string outPrefix;        // --out PREFIX (always a prefix)
    std::string compression;      // --compression (gz|zst, default: empty = plain text)
    std::string extractFile;      // --extract (SNP include list)
    std::string excludeFile;      // --exclude (SNP exclude list)
    std::string chrSpec;          // --chr (chromosome filter, e.g. "1-4,6,22")
    std::string admixBfilePrefix; // --admix-bfile
    std::string admixPhiFile;     // --admix-phi
    std::string mspFile;          // --rfmix-msp
    std::string admixTextPrefix;  // --admix-text-prefix
    std::string keepFile;         // --keep (subject include list)
    std::string removeFile;       // --remove (subject exclude list)
    std::string predListFile;     // --pred-list (Regenie / LDAK-KVIK pred.list for LOCO)
    // --pheno-transform {raw,int,standardize}; empty sentinel = "use context default"
    // (int when --pred-list given, raw otherwise; resolved in dispatch).
    std::string phenoTransform;
    // --spasqr-solver {conquer, qmme}: null-model SQR solver. Default qmme.
    std::string spasqrSolver = "qmme";
    // --spasqr-mode {score, wald}: score test (default, multi-tau CCT) vs
    // per-marker full-model Wald (β̂_G + SE per τ via M-estimation sandwich).
    std::string spasqrMode = "score";
    bool calAfCoef = false;
    bool calPairwiseIBD = false;
    bool calPhi = false;
    bool makeAbed = false;
    bool intPheno = false;
    double minMafIBD = 0.01;
    double refPrevalence = -1.0;
    double cutoff = 0.1;
    double spaCutoff = 2.0;
    double pvalCovAdjCut = 5e-5;
    double missingCutoff = 0.1;
    double minMafCutoff = 1e-5;
    double minMacCutoff = 10.0;
    double hweCutoff = 0.0;
    double outlierRatio = 1.5;
    double outlierAbsBound = 0.55;
    bool spagrmControlOutlier = false; // --spagrm-control-outlier (flag, no argument): enable iterative SPAGRM outlier-ratio adjustment (default off)
    int nthread = 1;
    int nSnpPerChunk = 8192;
    int compressionLevel = 0; // --compression-level (0 = library default)
    int nClusters = 0;        // 0 = auto (from --ref-af count)
    int leafKmeansNstart = 25; // --leaf-kmeans-nstart: K-means++ restarts
    uint64_t seed = 0;        // 0 = use std::random_device
    std::string leafClusterFile; // --leaf-cluster-file: pre-computed cluster
                                 // labels (two columns: IID, cluster).  When
                                 // set, LEAF skips K-means and uses these
                                 // labels verbatim.
};

// Entry point: parse argv, print help or dispatch the selected method.
int run(
    int argc,
    char *argv[]
);

// ── Shared helpers (used across cli translation units) ─────────────

Args parseArgs(
    int argc,
    char *argv[]
);

void printHelp(const std::string &topic);

} // namespace cli
