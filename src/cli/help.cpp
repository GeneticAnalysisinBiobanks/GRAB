// help.cpp — Help text generation from flag/method definitions

#include "cli/cli.hpp"
#include "cli/flags.hpp"

#include <cstdio>
#include <cstring>

namespace cli {

// Format "  --flag METAVAR    brief\n" with aligned columns
static void printFlag(const FlagDef* f) {
    char label[72];
    if (f->metavar)
        std::snprintf(label, sizeof(label), "%s %s", f->flag, f->metavar);
    else
        std::snprintf(label, sizeof(label), "%s", f->flag);
    std::fprintf(stderr, "  %-38s %s\n", label, f->brief);
}

// ── Short help (no topic) ──────────────────────────────────────────

static void printShortHelp() {
    std::fprintf(stderr,
        "GRAB -- Genome-wide Regression Analysis for Biobanks\n"
        "\n"
        "Usage:\n"
        "  grab --method METHOD  --bfile PREFIX  --null-resid FILE  --out FILE  [OPTIONS]\n"
        "  grab --cal-ind-af-coef   --bfile PREFIX  --eigenvec FILE  --out FILE\n"
        "  grab --cal-pairwise-ibd  --sp-grm-* FILE  --bfile PREFIX  --out FILE\n"
        "\n"
        "Methods:\n");
    for (const MethodDef* const* p = kAllMethods; *p; ++p)
        std::fprintf(stderr, "  %-14s %s\n", (*p)->name, (*p)->desc);
    std::fprintf(stderr, "\nUtility modes:\n");
    for (const MethodDef* const* p = kAllUtilModes; *p; ++p)
        std::fprintf(stderr, "  --%-20s %s\n", (*p)->name, (*p)->desc);
    std::fprintf(stderr,
        "\nRun 'grab --help <topic>' for details.  Topics:\n"
        "  SPACox  SPAGRM  SPAmix  SPAmixPlus  POLMM  SPAsqr  WtCoxG  LEAF\n"
        "  cal-ind-af-coef  cal-pairwise-ibd  options\n"
        "  null-resid  covar  eigenvec  ref-af  sp-grm  pairwise-ibd  ind-af-coef\n");
}

// ── Method help (generated from MethodDef) ─────────────────────────

static void printMethodHelp(const MethodDef* m) {
    bool isUtil = (std::strcmp(m->name, "cal-ind-af-coef") == 0 ||
                   std::strcmp(m->name, "cal-pairwise-ibd") == 0);
    if (isUtil)
        std::fprintf(stderr, "Mode: --%s -- %s\n", m->name, m->desc);
    else
        std::fprintf(stderr, "Method: %s -- %s\n", m->name, m->desc);

    std::fprintf(stderr, "\nRequired:\n");
    for (const FlagDef* const* p = m->required; *p; ++p)
        printFlag(*p);

    if (m->residNote)
        std::fprintf(stderr, "\n  --null-resid columns: %s\n", m->residNote);

    if (m->optional && m->optional[0]) {
        std::fprintf(stderr, "\nOptional:\n");
        for (const FlagDef* const* p = m->optional; *p; ++p)
            printFlag(*p);
    }

    if (m->extra)
        std::fprintf(stderr, "\n%s\n", m->extra);

    std::fprintf(stderr, "\nOutput:\n  %s\n", m->outputCols);
}

// ── Flag topic help (--help null-resid, --help sp-grm, etc.) ───────

static void printFlagHelp(const char* topic) {
    // Special case: "sp-grm" covers both flags
    if (std::strcmp(topic, "sp-grm") == 0) {
        std::fprintf(stderr, "Sparse GRM input -- provide exactly one of:\n\n");
        std::fprintf(stderr, "%s %s\n  %s\n", kSpGrmGrab.flag,
                     kSpGrmGrab.metavar, kSpGrmGrab.brief);
        if (kSpGrmGrab.fileInfo)
            std::fprintf(stderr, "\n%s\n", kSpGrmGrab.fileInfo);
        std::fprintf(stderr, "\n%s %s\n  %s\n", kSpGrmPlink2.flag,
                     kSpGrmPlink2.metavar, kSpGrmPlink2.brief);
        if (kSpGrmPlink2.fileInfo)
            std::fprintf(stderr, "\n%s\n", kSpGrmPlink2.fileInfo);
        return;
    }

    // Match "--X" flag names to topic "X"
    for (const FlagDef* const* p = kFileFlags; *p; ++p) {
        const char* fName = (*p)->flag + 2;   // skip "--"
        if (std::strcmp(topic, fName) == 0) {
            std::fprintf(stderr, "%s %s\n  %s\n", (*p)->flag,
                         (*p)->metavar ? (*p)->metavar : "", (*p)->brief);
            if ((*p)->fileInfo)
                std::fprintf(stderr, "\n%s\n", (*p)->fileInfo);
            return;
        }
    }

    std::fprintf(stderr, "Unknown help topic: '%s'\n\n", topic);
    printShortHelp();
}

// ── Options help (all flags grouped) ───────────────────────────────

static void printOptionsHelp() {
    std::fprintf(stderr, "All options:\n\nFile inputs:\n");
    for (const FlagDef* const* p = kInputFlags; *p; ++p)
        printFlag(*p);
    std::fprintf(stderr, "\nNumeric options:\n");
    for (const FlagDef* const* p = kNumericFlags; *p; ++p)
        printFlag(*p);
}

// ── Public entry point ─────────────────────────────────────────────

void printHelp(const std::string& topic) {
    if (topic == "__short__") { printShortHelp(); return; }
    if (topic == "options")   { printOptionsHelp(); return; }

    // Check method names
    for (const MethodDef* const* p = kAllMethods; *p; ++p)
        if (topic == (*p)->name) { printMethodHelp(*p); return; }
    for (const MethodDef* const* p = kAllUtilModes; *p; ++p)
        if (topic == (*p)->name) { printMethodHelp(*p); return; }

    // Check flag topics
    printFlagHelp(topic.c_str());
}

} // namespace cli
