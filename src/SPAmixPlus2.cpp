
#include <RcppArmadillo.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>

using namespace Rcpp;
using namespace arma;

#include "SPAmixPlus.h"
#include "Main.h"

// Extern declaration for the global pointer defined in Main.cpp
extern SPAmixPlus::SPAmixPlusClass* ptr_gSPAmixPlusobj;

// [[Rcpp::export]]
void exportAFModelInCPP(
    std::string t_method,
    std::string t_genoType,
    std::vector<uint64_t> t_genoIndex,
    std::string t_outputFile,
    std::string t_impute_method
) {
    if (t_method != "SPAmixPlus" || ptr_gSPAmixPlusobj == nullptr) {
        Rcpp::stop("exportAFModelInCPP only supports SPAmixPlus method and initialized object.");
    }
    
    // Create file if not exists, to check/ensure we can open in random access mode
    {
        std::ifstream testOpen(t_outputFile);
        if(!testOpen.good()){
            std::ofstream create(t_outputFile, std::ios::binary);
            create.close();
        }
    }

    // Open binary file in update mode
    std::fstream outFile(t_outputFile, std::ios::binary | std::ios::in | std::ios::out);
    if (!outFile) Rcpp::stop("Failed to open output file " + t_outputFile);
    
    int nPCs = ptr_gSPAmixPlusobj->getNPCs();
    long long recordSize = sizeof(int) + (long long)(nPCs + 1) * sizeof(double);
    // Write header? No, just raw records to allow seek.
    // Assuming t_genoIndex matches the file structure or we append.
    // If we append, we can't seek easily unless we know order.
    // We assume the user manages the file creation properly (e.g. per chunk output).
    
    int q = t_genoIndex.size();
    
    // Progress
    int nPercent = q / 100;
    if(nPercent == 0) nPercent = 1;
    
    for(int i = 0; i < q; i++){
        if(i % nPercent == 0){
             Rcpp::checkUserInterrupt();
             Rprintf("Processed %d / %d markers (%.1f%%) \r", i, q, (100.0 * i) / q); 
        }

        uint64_t genoIndex = t_genoIndex[i];
        
        std::string marker;
        std::string chr;
        uint32_t pd;
        std::string ref;
        std::string alt;
        double altFreq;
        double altCounts; // Added
        double missingRate;
        std::vector<uint32_t> indexForMissing;
        std::vector<uint32_t> indexForNonZero;
        arma::vec GVec;
        
        double imputeInfo; // QC metrics

        // Use Unified getter
        GVec = Unified_getOneMarker(
            t_genoType,           // Genotype file format ("PLINK", "BGEN")
            genoIndex,            // Marker index in genotype file
            ref,                  // Reference allele (output)
            alt,                  // Alternate allele (output)
            marker,               // Marker ID (output)
            pd,                   // Physical position (output)
            chr,                  // Chromosome (output)
            altFreq,              // Alternate allele frequency (output)
            altCounts,            // Alternate allele count (output)
            missingRate,          // Missing genotype rate (output)
            imputeInfo,           // Imputation quality score (output)
            true,                 // Whether to output missing data indices
            indexForMissing,      // Indices of missing genotypes (output)
            false,                // Whether to output only non-zero genotypes
            indexForNonZero       // Indices of non-zero genotypes (output)
        );
        // Impute
        bool flip = imputeGenoAndFlip(
            GVec,                // Genotype vector
            altFreq,             // Alternate allele frequency
            indexForMissing,     // Indices of missing genotypes
            missingRate,         // Missing genotype rate
            t_impute_method,     // Imputation method ("mean", "minor", "drop")
            t_method             // Statistical method
        );
            
        // Compute Model
        SPAmixPlus::SPAmixPlusClass::AFModelInfo model = ptr_gSPAmixPlusobj->computeAFModel(GVec, altFreq);
        
        // Write to file with random access
        long long pos = (long long)genoIndex * recordSize;
        outFile.seekp(pos, std::ios::beg);
        
        // Format: [Status: int] [Betas: double array of size K+1]
        outFile.write(reinterpret_cast<const char*>(&model.status), sizeof(int));
        outFile.write(reinterpret_cast<const char*>(model.betas.memptr()), (nPCs + 1) * sizeof(double));
    }
    outFile.close();
}



