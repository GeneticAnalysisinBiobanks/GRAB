// help.cpp — Help text generation from flag/method definitions

#include "cli/cli.hpp"
#include "cli/flags.hpp"

#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>

namespace cli {

// Lowercase a string for display.  Method names are stored mixed-case in
// MethodDef::name (e.g. "SPACox", "WtCoxG") for downstream string-equality
// dispatch and output-filename suffixes, but the user-facing --help text
// uniformly lowercases them per the project convention.
static std::string toLower(const char *s) {
    std::string out(s ? s : "");
    for (auto &c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

// Format "  --flag METAVAR    brief\n" with aligned columns.
// The label buffer must accommodate the combined display strings used for
// the genotype-input and sparse-GRM "exactly one" entries, which exceed
// 80 characters.
static void printFlag(const FlagDef *f) {
    char label[160];
    if (f->metavar)
        std::snprintf(label, sizeof(label), "%s %s", f->flag, f->metavar);
    else
        std::snprintf(label, sizeof(label), "%s", f->flag);
    std::fprintf(stderr, "  %-38s %s\n", label, f->brief);
}

// ── Short help (no topic) ──────────────────────────────────────────

static void printShortHelp() {
    std::fputs(
        R"(GRAB 2.0-alpha -- Genome-Wide Robust Analysis for Biobank Data

Usage:
  grab2 --method spacox|spagrm|spamix|sageld|spasqr|wtcoxg|leaf \
        --bfile PREFIX --pheno FILE \
        --pheno-name COL_IDS \
        --covar-name COL_IDS \
        [--sp-grm-plink2 FILE] \
        [options] \
        --out PREFIX

For more help, run 'grab2 --help TOPIC':
  Methods:       spacox  spagrm  sageld  spamix  spasqr  wtcoxg  leaf
  Utilities:     cal-af-coef  cal-pairwise-ibd  int-pheno
  Main inputs:   geno  pheno  covar
  Other inputs:  sp-grm  pairwise-ibd  ind-af-coef  ref-af
  List options:  options
)",
        stderr);
}

// ── Method help (generated from MethodDef) ─────────────────────────

static void printMethodHelp(const MethodDef *m) {
    bool isUtil = (std::strcmp(m->name, "cal-af-coef") == 0 || std::strcmp(m->name, "cal-pairwise-ibd") == 0 ||
                   std::strcmp(m->name, "cal-phi") == 0 || std::strcmp(m->name, "make-abed") == 0 ||
                   std::strcmp(m->name, "int-pheno") == 0);
    const std::string nameLc = toLower(m->name);
    if (isUtil)
        std::fprintf(stderr, "Mode: --%s -- %s\n", nameLc.c_str(), m->desc);
    else
        std::fprintf(stderr, "Method: %s -- %s\n", nameLc.c_str(), m->desc);

    std::fprintf(stderr, "\nRequired:\n");
    for (const FlagDef *const *p = m->required; *p; ++p)
        printFlag(*p);

    if (m->residNote) std::fprintf(stderr, "\n  Residual columns (via --pheno): %s\n", m->residNote);

    if (m->optional && m->optional[0]) {
        std::fprintf(stderr, "\nOptional:\n");
        for (const FlagDef *const *p = m->optional; *p; ++p)
            printFlag(*p);
    }

    if (m->extra) std::fprintf(stderr, "\n%s\n", m->extra);

    std::fprintf(stderr, "\nOutput:\n  %s\n", m->outputCols);
}

// ── Flag topic help (--help pheno, --help sp-grm, etc.) ───────

static void printFlagHelp(const char *topic) {
    // Special case: "sp-grm" describes the sparse-GRM input.  Two
    // mutually exclusive formats are accepted; the plink2 .grm.sp form
    // is the documented default and is the only one mentioned in the
    // short usage and in `--help options`.  The legacy GRAB three-column
    // form remains supported and is shown here for users who already have
    // a file in that format.
    if (std::strcmp(topic, "sp-grm") == 0) {
        std::fprintf(stderr, "Sparse GRM input -- provide exactly one of:\n\n");
        std::fprintf(stderr, "%s %s\n  %s\n", kSpGrmPlink2.flag, kSpGrmPlink2.metavar, kSpGrmPlink2.brief);
        if (kSpGrmPlink2.fileInfo) std::fprintf(stderr, "\n%s\n", kSpGrmPlink2.fileInfo);
        std::fprintf(stderr, "\n%s %s\n  %s\n", kSpGrmGrab.flag, kSpGrmGrab.metavar, kSpGrmGrab.brief);
        if (kSpGrmGrab.fileInfo) std::fprintf(stderr, "\n%s\n", kSpGrmGrab.fileInfo);
        return;
    }

    // Special case: "geno" describes the genotype-input alternatives.
    if (std::strcmp(topic, "geno") == 0) {
        std::fputs(
            R"(Genotype input -- provide exactly one of:

  --bfile PREFIX               PLINK binary genotype prefix (.bed/.bim/.fam)
  --pfile PREFIX               PLINK2 genotype prefix (.pgen/.pvar/.psam)
  --vcf FILE                   VCF text file (.vcf or BGZF-compressed .vcf.gz)
  --bcf FILE                   BCF2 binary file
  --bgen FILE <REF/ALT mode>   BGEN file; mode is ref-first | ref-last | ref-unknown

Per-format details follow.
)",
            stderr);
        std::fprintf(stderr, "\n%s\n  %s\n", kVcf.flag, kVcf.brief);
        if (kVcf.fileInfo) std::fprintf(stderr, "\n%s\n", kVcf.fileInfo);
        std::fprintf(stderr, "\n%s\n  %s\n", kBcf.flag, kBcf.brief);
        if (kBcf.fileInfo) std::fprintf(stderr, "\n%s\n", kBcf.fileInfo);
        std::fprintf(stderr, "\n%s\n  %s\n", kBgen.flag, kBgen.brief);
        if (kBgen.fileInfo) std::fprintf(stderr, "\n%s\n", kBgen.fileInfo);
        return;
    }

    // Match "--X" flag names to topic "X" across all flag groups.  The
    // extended-help field (`fileInfo`) was originally introduced for file
    // flags but is now also used by numeric flags (e.g. --outlier-iqr-
    // threshold) that need a paragraph of context beyond the one-line
    // brief.
    const FlagDef *const *flagGroups[] = { kFileFlags, kInputFlags,
                                           kNumericFlags, nullptr };
    for (const FlagDef *const *const *g = flagGroups; *g; ++g) {
        for (const FlagDef *const *p = *g; *p; ++p) {
            const char *fName = (*p)->flag + 2; // skip "--"
            if (std::strcmp(topic, fName) == 0) {
                std::fprintf(stderr, "%s %s\n  %s\n", (*p)->flag,
                             (*p)->metavar ? (*p)->metavar : "", (*p)->brief);
                if ((*p)->fileInfo) std::fprintf(stderr, "\n%s\n", (*p)->fileInfo);
                return;
            }
        }
    }

    std::fprintf(stderr, "Unknown help topic: '%s'\n\n", topic);
    printShortHelp();
}

// ── Options help (all flags grouped) ───────────────────────────────

static void printOptionsHelp() {
    std::fprintf(stderr, "All options:\n\nFile inputs:\n");
    for (const FlagDef *const *p = kInputFlags; *p; ++p)
        printFlag(*p);
    std::fprintf(stderr, "\nNumeric options:\n");
    for (const FlagDef *const *p = kNumericFlags; *p; ++p)
        printFlag(*p);
}

// ── Public entry point ─────────────────────────────────────────────

void printHelp(const std::string &topic) {
    if (topic == "__short__") {
        printShortHelp();
        return;
    }
    if (topic == "options") {
        printOptionsHelp();
        return;
    }

    // Check method names.  Iterate kVisibleMethods / kVisibleUtilModes
    // (which exclude SPAmixPlus, SPAmixLocalPlus, --cal-phi, --make-abed)
    // so those hidden modes also fall through to "Unknown help topic"
    // when looked up by name.  The dispatcher continues to recognize them
    // via kAllMethods / kAllUtilModes.  Matching is case-insensitive: the
    // user may type any case (e.g. spacox, SPACox, SPACOX).
    const std::string topicLc = toLower(topic.c_str());
    for (const MethodDef *const *p = kVisibleMethods; *p; ++p)
        if (topicLc == toLower((*p)->name)) {
            printMethodHelp(*p);
            return;
        }
    for (const MethodDef *const *p = kVisibleUtilModes; *p; ++p)
        if (topicLc == toLower((*p)->name)) {
            printMethodHelp(*p);
            return;
        }

    // Check flag topics
    printFlagHelp(topic.c_str());
}

} // namespace cli
