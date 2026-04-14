// loco.hpp — Generic LOCO (Leave-One-Chromosome-Out) engine
//
// Iterates chromosomes serially, calling a user-supplied callback to
// rebuild PhenoTasks per chromosome (because LOCO covariates differ),
// then processes that chromosome's marker chunks via multi-phenotype
// work-stealing parallelism with persistent per-phenotype writers.
//
// Applicable to any method that fits a null model from covariates:
// SPAsqr, SPACox, POLMM, WtCoxG, etc.
#pragma once

#include "engine/marker.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class GenoMeta;

// ======================================================================
// LocoData — LOCO polygenic scores from Regenie step 1
// ======================================================================

struct LocoData {
    // scores[phenoName][chrStr] = VectorXd of LOCO values (in usedIIDs order)
    std::unordered_map<std::string,
                       std::unordered_map<std::string, Eigen::VectorXd> > scores;

    // Load from a pred.list file (space-separated: phenoName locoFilePath).
    // Validates that all phenoNames are present in the pred.list.
    // Remaps subject columns from .loco FID_IID format to usedIIDs order.
    static LocoData load(
        const std::string &predListFile,
        const std::vector<std::string> &phenoNames,
        const std::vector<std::string> &usedIIDs,
        const std::vector<std::string> &famIIDs
    );

    // Return chromosome strings present for ALL phenotypes.
    std::unordered_set<std::string> availableChromosomes() const;

};

// ======================================================================
// LocoTaskBuilder — per-chromosome callback to rebuild PhenoTasks
// ======================================================================

// Called once per chromosome.  The callback should fill `tasks` with
// K PhenoTasks, each containing a freshly built MethodBase for that
// chromosome's LOCO-adjusted covariates.
using LocoTaskBuilder = std::function<void (
                                          const std::string &chr,
                                          std::vector<PhenoTask> &tasks
                                      )>;

// ======================================================================
// locoEngine — generic LOCO marker association engine
// ======================================================================

// Iterates chromosomes serially.  For each chromosome:
//   1. Calls buildTasks(chr, tasks) to get K PhenoTasks
//   2. Processes that chromosome's chunks via work-stealing workers
//   3. Writes results to persistent per-phenotype output files
//
// locoChroms: chromosomes available in LOCO data (used for O4 filtering)
// phenoNames: used to build output file paths and write headers
void locoEngine(
    const GenoMeta &genoData,
    const std::unordered_set<std::string> &locoChroms,
    const std::vector<std::string> &phenoNames,
    LocoTaskBuilder buildTasks,
    const std::string &outPrefix,
    const std::string &methodName,
    const std::string &compression,
    int compressionLevel,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff
);
