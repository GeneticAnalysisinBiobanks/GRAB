// Main2.cpp - Async I/O + Multithreading for GRAB.Marker2
// Simple implementation: reader thread -> worker pool -> writer thread

#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>
#include <string>
#include <fstream>

#include "BGEN.h"
#include "PLINK.h"
#include "UTIL.h"
#include "SPAmix.h"
#include "SPACox.h"
#include "Main.h"

//==============================================================================
// SECTION 1: OBJECT DECLARATIONS
//==============================================================================

// commont variables for marker level analysis
static std::string g_impute_method;          // Imputation method: "mean", "minor", or "drop"
static double g_missingRate_cutoff;          // Maximum allowed missing rate for markers
static double g_marker_minMAF_cutoff;         // Minimum Minor Allele Frequency for single markers
static double g_marker_minMAC_cutoff;         // Minimum Minor Allele Count for single markers
static unsigned int g_omp_num_threads;       // Number of OpenMP threads for parallel processing

PLINK::PlinkClass* ptr_gPLINKobj = nullptr;
static SPAmix::SPAmixClass* ptr_gSPAmixobj = nullptr;


//==============================================================================
// SECTION 2: THREAD-SAFE IO BUFFER CLASSES
//==============================================================================

// Thread-safe input buffer for genotype data
class GenoInputBuffer {
private:
    std::queue<int> markerQueue;
    std::mutex mtx;
    std::condition_variable cv;
    bool finished;
    int maxSize;

public:
    GenoInputBuffer(int size) : finished(false), maxSize(size) {}

    void push(int markerIdx) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return markerQueue.size() < maxSize || finished; });
        if (!finished) {
            markerQueue.push(markerIdx);
        }
        cv.notify_all();
    }

    bool pop(int& markerIdx) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return !markerQueue.empty() || finished; });
        if (markerQueue.empty()) {
            return false;
        }
        markerIdx = markerQueue.front();
        markerQueue.pop();
        cv.notify_all();
        return true;
    }

    void setFinished() {
        std::unique_lock<std::mutex> lock(mtx);
        finished = true;
        cv.notify_all();
    }

    bool isFinished() {
        std::unique_lock<std::mutex> lock(mtx);
        return finished && markerQueue.empty();
    }
};


// Thread-safe output buffer for results
class ResultOutputBuffer {
private:
    struct Result {
        int markerIdx;
        std::string line;
        
        Result(int idx, const std::string& l) : markerIdx(idx), line(l) {}
        
        bool operator<(const Result& other) const {
            return markerIdx > other.markerIdx; // Min-heap by markerIdx
        }
    };

    std::priority_queue<Result> resultHeap;
    std::mutex mtx;
    std::condition_variable cv;
    int nextExpectedIdx;
    bool finished;
    int maxSize;

public:
    ResultOutputBuffer(int size) : nextExpectedIdx(0), finished(false), maxSize(size) {}

    void push(int markerIdx, const std::string& line) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return resultHeap.size() < maxSize || finished; });
        if (!finished) {
            resultHeap.push(Result(markerIdx, line));
        }
        cv.notify_all();
    }

    bool popInOrder(std::string& line) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { 
            return (!resultHeap.empty() && resultHeap.top().markerIdx == nextExpectedIdx) || finished; 
        });
        
        if (resultHeap.empty()) {
            return false;
        }
        
        if (resultHeap.top().markerIdx == nextExpectedIdx) {
            line = resultHeap.top().line;
            resultHeap.pop();
            nextExpectedIdx++;
            cv.notify_all();
            return true;
        }
        
        return false;
    }

    void setFinished() {
        std::unique_lock<std::mutex> lock(mtx);
        finished = true;
        cv.notify_all();
    }

    bool isFinished() {
        std::unique_lock<std::mutex> lock(mtx);
        return finished && resultHeap.empty();
    }
};


// ==============================================================================
// SECTION 3: THREAD FUNCTIONS
// ==============================================================================

// Reader thread: read genotype data and push to input buffer
static void readerThread(
    GenoInputBuffer* inputBuffer,
    int nMarkers
) {
    for (int i = 0; i < nMarkers; i++) {
        inputBuffer->push(i);
    }
    inputBuffer->setFinished();
}


// Writer thread: pop results in order and write to file
static void writerThread(
    ResultOutputBuffer* outputBuffer,
    const std::string& outputFile,
    const std::string& header
) {
    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        Rcpp::stop("Cannot open output file: " + outputFile);
    }
    
    // Write header
    outFile << header << "\n";
    
    // Write results in order
    std::string line;
    while (outputBuffer->popInOrder(line)) {
        outFile << line << "\n";
    }
    
    outFile.close();
}


// Worker thread: pop from input buffer, test marker, push to output buffer
static void workerThread(
    GenoInputBuffer* inputBuffer,
    ResultOutputBuffer* outputBuffer
) {
    int markerIdx;
    
    while (inputBuffer->pop(markerIdx)) {
        
        // Get genotype data for one marker
        std::string ref, alt, marker, chr;
        uint32_t pd;
        double altFreq, altCounts, missingRate, imputeInfo;
        std::vector<uint32_t> indexForMissing, indexForNonZero;
        
        arma::vec GVec = Unified_getOneMarker(
            g_genoType,
            markerIdx,
            ref,
            alt,
            marker,
            pd,
            chr,
            altFreq,
            altCounts,
            missingRate,
            imputeInfo,
            true,               // output missing indices
            indexForMissing,
            false,              // don't need non-zero indices
            indexForNonZero
        );
        
        // Impute missing genotypes
        bool flip = imputeGenoAndFlip(
            GVec,
            altFreq,
            indexForMissing,
            missingRate,
            g_impute_method,
            g_method
        );
        
        std::string resultLine;
        std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;
        
        if (g_method == "SPAmix") {
            
            double pval = g_SPAmixObj->getMarkerPval(GVec, altFreq);
            arma::vec pvalVec = g_SPAmixObj->getpvalVec();
            arma::vec zScoreVec = g_SPAmixObj->getzScoreVec();
            
            for (int p = 0; p < g_nPheno; p++) {
                if (p > 0) resultLine += "\n";
                resultLine += "pheno_" + std::to_string(p + 1) + "\t";
                resultLine += marker + "\t";
                resultLine += info + "\t";
                resultLine += std::to_string(altFreq) + "\t";
                resultLine += std::to_string(altCounts) + "\t";
                resultLine += std::to_string(missingRate) + "\t";
                resultLine += std::to_string(pvalVec[p]) + "\t";
                resultLine += std::to_string(zScoreVec[p]);
            }
            
        } else if (g_method == "SPACox") {
            
            double pval, zScore, Beta, seBeta;
            bool isSPAConverge;
            pval = g_SPACoxObj->getMarkerPval(GVec, altFreq, zScore);
            
            resultLine = marker + "\t";
            resultLine += info + "\t";
            resultLine += std::to_string(altFreq) + "\t";
            resultLine += std::to_string(altCounts) + "\t";
            resultLine += std::to_string(missingRate) + "\t";
            resultLine += std::to_string(pval) + "\t";
            resultLine += std::to_string(Beta * (1 - 2 * flip)) + "\t";
            resultLine += std::to_string(seBeta) + "\t";
            resultLine += (isSPAConverge ? "TRUE" : "FALSE");
        }
        
        outputBuffer->push(markerIdx, resultLine);
    }
}




// ==============================================================================
// SECTION 4: MAIN FUNCTION
// ==============================================================================

void setPLINKobjInCPP(
  std::string t_bimFile,                    // Path to PLINK .bim file (marker information)
  std::string t_famFile,                    // Path to PLINK .fam file (sample information)
  std::string t_bedFile,                    // Path to PLINK .bed file (binary genotype data)
  std::vector<std::string> t_SampleInModel, // Vector of sample IDs to include in analysis
  std::string t_AlleleOrder                 // Allele ordering convention ("alt-first" or "ref-first")
) {

}

// [[Rcpp::export]]
void mainMarkerInCPP2(
  // IO and threading parameters
  const std::string outputFile,
  const int inputBufferSize,
  const int outputBufferSize,
  const int nWorkers,

  // PLINK inputs
  const std::string t_bimFile,
  const std::string t_famFile,
  const std::string t_bedFile,
  const std::vector<std::string> t_SampleInModel,
  const std::string t_AlleleOrder,

  // setMarker_GlobalVarsInCPP
  const std::string t_impute_method,
  const double t_missing_cutoff,
  const double t_min_maf_marker,
  const double t_min_mac_marker,

  // setMarker.SPAmix
  arma::mat t_resid,
  arma::mat t_PCs,
  int t_N,
  double t_SPA_Cutoff,
  Rcpp::List t_outlierList
) {


  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;

  const std::string t_genoType = "PLINK"; // For now, only PLINK is supported
  const std::string t_method = "SPAmix"; // For now, only SPAmix is supported

  if (ptr_gPLINKobj)
    delete ptr_gPLINKobj;

  ptr_gPLINKobj = new PLINK::PlinkClass(
    t_bimFile,                              // Path to PLINK .bim file (marker information)
    t_famFile,                              // Path to PLINK .fam file (sample information)
    t_bedFile,                              // Path to PLINK .bed file (binary genotype data)
    t_SampleInModel,                        // Vector of sample IDs to include in analysis
    t_AlleleOrder                           // Allele ordering convention ("alt-first" or "ref-first")
  );

  int n = ptr_gPLINKobj->getN();
  Rcpp::Rcout << "    Number of subjects with genotype: " << n << std::endl;

  const int q = wc -l t_bimFile; // Number of markers in the .bim file
  Rcpp::Rcout << "    Number of markers to analyze: " << q << std::endl;


  if (ptr_gSPAmixobj)
    delete ptr_gSPAmixobj;

  ptr_gSPAmixobj = new SPAmix::SPAmixClass(
    t_resid,
    t_PCs, 
    t_N,
    t_SPA_Cutoff,
    t_outlierList
  );
  const int const Npheno = ptr_gSPAmixobj->getNpheno();

  // Create buffers
  GenoInputBuffer inputBuffer(inputBufferSize);
  ResultOutputBuffer outputBuffer(outputBufferSize);
  
  // Launch reader thread
  std::thread reader(readerThread, &inputBuffer, nMarkers);
  // Launch writer thread
  std::thread writer(writerThread, &outputBuffer, outputFile, header);
  // Launch worker threads
  std::vector<std::thread> workers;
  for (int i = 0; i < nWorkers; i++) {
      workers.emplace_back(workerThread, &inputBuffer, &outputBuffer);
  }

  // Print header
  print("Pheno\tMarker\tInfo\tAltFreq\tAltCounts\tMissingRate\tPvalue1\tzScore1\t...\tPvalue$Npheno\tzScore$Npheno\n");

  for (int i = 0; i < q; i++) {

    if (i % 1000 == 0)
      Rcpp::Rcout << "    Completed " << i << "/" << q << " markers in the chunk." << std::endl;

    // Variables to store marker-specific information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;  // Physical position
    bool flip = false;  // Whether a locus is flipped

    uint64_t gIndex;

    // Extract genotype vector and marker information
    arma::vec GVec = Unified_getOneMarker(
      t_genoType,                      // Genotype file format ("PLINK", "BGEN")
      gIndex,                          // Marker index in genotype file
      ref,                             // Reference allele (output)
      alt,                             // Alternate allele (output)
      marker,                          // Marker ID (output)
      pd,                              // Physical position (output)
      chr,                             // Chromosome (output)
      altFreq,                         // Alternate allele frequency (output)
      altCounts,                       // Alternate allele count (output)
      missingRate,                     // Missing genotype rate (output)
      imputeInfo,                      // Imputation quality score (output)
      true,                            // Whether to output missing data indices
      indexForMissing,                 // Indices of missing genotypes (output)
      false,                           // Whether to output only non-zero genotypes
      indexForNonZero                  // Indices of non-zero genotypes (output)
    );
    int n = GVec.size();  // Sample size


    // Format marker information string (CHR:POS:REF:ALT)
    std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;

    // Store basic marker information in output vectors
    print(marker, info, altFreq, altCounts, missingRate);

    // Calculate MAF and MAC for QC
    double MAF = std::min(
      altFreq,      // Alternate allele frequency
      1 - altFreq   // Reference allele frequency
    );  // MAF is always ≤ 0.5
    double MAC = 2 * MAF * n * (1 - missingRate); // Account for diploid genotypes and missing data

    // Quality Control: check if marker passes MAF, MAC, and missing rate thresholds
    if ((missingRate > g_missingRate_cutoff) ||
        (MAF < g_marker_minMAF_cutoff) ||
        (MAC < g_marker_minMAC_cutoff))
      continue;

    // Check UTIL.cpp
    flip = imputeGenoAndFlip(
      GVec,                // Genotype vector
      altFreq,             // Alternate allele frequency
      indexForMissing,     // Indices of missing genotypes
      missingRate,         // Missing genotype rate
      g_impute_method,     // Imputation method ("mean", "minor", "drop")
      t_method             // Statistical method
    );

    double Beta, seBeta, pval, zScore, hwepval;
    pval = ptr_gSPAmixobj->getMarkerPval(GVec, altFreq);
    arma::vec pvalVecTemp = ptr_gSPAmixobj->getpvalVec();
    arma::vec zScoreVecTemp = ptr_gSPAmixobj->getzScoreVec();
    
    for (int j = 0; j < Npheno; j++) {
      print(pvalVecTemp.at(j), zScoreVecTemp.at(j));
    }
    
    print("\n");
    
  } // End of one marker loop

  // Wait for reader to finish
  reader.join();
  
  // Wait for all workers to finish
  for (auto& w : workers) {
      w.join();
  }
  
  // Signal output buffer finished
  outputBuffer.setFinished();
    
  // Wait for writer to finish
  writer.join();
    
  // Cleanup
  delete ptr_gPLINKobj;
  delete ptr_gSPAmixobj;

  Rcpp::Rcout << "Analysis complete. Results written to: " << outputFile << std::endl;
}

