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
#include <memory>

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
static std::string g_method;
static std::string g_genoType;
static std::string g_impute_method;          // Imputation method: "mean", "minor", or "drop"

static double g_missingRate_cutoff;          // Maximum allowed missing rate for markers
static double g_marker_minMAF_cutoff;         // Minimum Minor Allele Frequency for single markers
static double g_marker_minMAC_cutoff;         // Minimum Minor Allele Count for single markers
static int g_nPheno;

static PLINK::PlinkClass* ptr_gPLINKobj = nullptr;
// static BGEN::BgenClass* ptr_gBGENobj = nullptr;
// static std::list<SPAmix::SPAmixClass*> ptr_gSPAmixobjList;

// static SPAmix::SPAmixClass* ptr_gSPAmixobj = nullptr;
// static SPACox::SPACoxClass* ptr_gSPACoxobj = nullptr;

// Container for one marker's genotypes and metadata
struct MarkerData {
    int markerIdx;
    std::string ref, alt, marker, chr;
    uint32_t pd = 0;
    double altFreq = 0.0, altCounts = 0.0, missingRate = 0.0, imputeInfo = 0.0;
    std::vector<uint32_t> indexForMissing;
    std::vector<uint32_t> indexForNonZero;
    arma::vec GVec;
};


//==============================================================================
// SECTION 2: THREAD-SAFE IO BUFFER CLASSES
//==============================================================================

// Thread-safe input buffer for genotype data
class GenoInputBuffer {
private:
    std::queue<std::shared_ptr<MarkerData>> markerQueue;
    std::mutex mtx;
    std::condition_variable cv;
    bool finished;
    int maxSize;

public:
    GenoInputBuffer(int size) : finished(false), maxSize(size) {}

    void push(std::shared_ptr<MarkerData> marker) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return markerQueue.size() < (std::size_t)maxSize || finished; });
        if (!finished) {
            markerQueue.push(std::move(marker));
        }
        cv.notify_all();
    }

    bool pop(std::shared_ptr<MarkerData>& marker) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return !markerQueue.empty() || finished; });
        if (markerQueue.empty()) {
            return false;
        }
        marker = std::move(markerQueue.front());
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
        cv.wait(lock, [this] { return resultHeap.size() < (std::size_t)maxSize || finished; });
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
        auto marker = std::make_shared<MarkerData>();
        marker->markerIdx = i;

        // DIRECT CALL TO AVOID Main.cpp static variable issues
        if (g_genoType == "PLINK") {
            if (!ptr_gPLINKobj) throw std::runtime_error("ptr_gPLINKobj is null");
            marker->GVec = ptr_gPLINKobj->getOneMarker(
                marker->markerIdx,
                marker->ref,
                marker->alt,
                marker->marker,
                marker->pd,
                marker->chr,
                marker->altFreq,
                marker->altCounts,
                marker->missingRate,
                marker->imputeInfo,
                true,                // output missing indices
                marker->indexForMissing,
                false,               // don't need non-zero indices
                marker->indexForNonZero,
                true                 // isTrueGenotype for PLINK
            );
        } else {
            // Fallback for BGEN or others if supported in future
                throw std::runtime_error("Only PLINK supported in Main2 async mode currently.");
        }

        inputBuffer->push(std::move(marker));
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


static void workerSPAmix(
    GenoInputBuffer* inputBuffer,
    ResultOutputBuffer* outputBuffer,
    SPAmix::SPAmixClass* ptr_gSPAmixobj
) {
    std::shared_ptr<MarkerData> marker;
    while (inputBuffer->pop(marker)) {

        // QC
        double MAF = std::min(marker->altFreq, 1 - marker->altFreq);
        double MAC = 2 * MAF * marker->GVec.size() * (1 - marker->missingRate);
        if ((marker->missingRate > g_missingRate_cutoff) ||
            (MAF < g_marker_minMAF_cutoff) ||
            (MAC < g_marker_minMAC_cutoff)) {
            outputBuffer->push(marker->markerIdx, ""); // skip but keep order
            continue;
        }

        // Impute and possibly flip
        imputeGenoAndFlip(
            marker->GVec,
            marker->altFreq,
            marker->indexForMissing,
            marker->missingRate,
            g_impute_method,
            g_method
        );

        std::string info = marker->chr + ":" + std::to_string(marker->pd) + ":" + marker->ref + ":" + marker->alt;
        std::string resultLine;

        ptr_gSPAmixobj->getMarkerPval(marker->GVec, marker->altFreq);
        arma::vec pvalVec = ptr_gSPAmixobj->getpvalVec();
        arma::vec zScoreVec = ptr_gSPAmixobj->getzScoreVec();

        for (int p = 0; p < g_nPheno; p++) {
            if (p > 0) resultLine += "\n";
            resultLine += "pheno_" + std::to_string(p + 1) + "\t";
            resultLine += marker->marker + "\t";
            resultLine += info + "\t";
            resultLine += std::to_string(marker->altFreq) + "\t";
            resultLine += std::to_string(marker->altCounts) + "\t";
            resultLine += std::to_string(marker->missingRate) + "\t";
            resultLine += std::to_string(pvalVec[p]) + "\t";
            resultLine += std::to_string(zScoreVec[p]);

        }

        outputBuffer->push(marker->markerIdx, resultLine);
    }
}


// [[Rcpp::export]]
void mainMarkerInCPP2(
        Rcpp::List control,
        Rcpp::List objNull
) {
    auto requireString = [](const Rcpp::List& lst, const char* name) {
        if (!lst.containsElementNamed(name)) {
            Rcpp::stop(std::string("Missing required string field: ") + name);
        }
        return Rcpp::as<std::string>(lst[name]);
    };

    auto requireInt = [](const Rcpp::List& lst, const char* name) {
        if (!lst.containsElementNamed(name)) {
            Rcpp::stop(std::string("Missing required integer field: ") + name);
        }
        return Rcpp::as<int>(lst[name]);
    };

    auto requireDouble = [](const Rcpp::List& lst, const char* name) {
        if (!lst.containsElementNamed(name)) {
            Rcpp::stop(std::string("Missing required numeric field: ") + name);
        }
        return Rcpp::as<double>(lst[name]);
    };

    // Pull control parameters
    const std::string outputFile      = requireString(control, "outputFile");
    const int inputBufferSize         = requireInt(control, "inputBufferSize");
    const int outputBufferSize        = requireInt(control, "outputBufferSize");
    const int nWorkers                = requireInt(control, "nWorkers");

    g_impute_method                   = requireString(control, "impute_method");
    g_missingRate_cutoff              = requireDouble(control, "missing_cutoff");
    g_marker_minMAF_cutoff            = requireDouble(control, "min_maf_marker");
    g_marker_minMAC_cutoff            = requireDouble(control, "min_mac_marker");
    g_genoType                        = requireString(control, "genoType");
    g_method                          = requireString(control, "method");
    const double t_SPA_Cutoff         = requireDouble(control, "SPA_Cutoff");

    const std::string t_bimFile       = requireString(control, "bimfile");
    const std::string t_famFile       = requireString(control, "famfile");
    const std::string t_bedFile       = requireString(control, "bedfile");
    const std::string t_AlleleOrder   = requireString(control, "AlleleOrder");

    const std::vector<std::string> t_SampleInModel = Rcpp::as< std::vector<std::string> >(objNull["subjData"]);
    
    // Convert to thread-safe data structures if needed, but arma::mat is generally safe for read-only if not modified
    // However, Rcpp objects (like Rcpp::List) are NOT safe to access from threads.
    // The previous fix used std::vector<OutlierData> in SPAmix. Let's see if user reverted that.
    // The user code passes t_outlierList (Rcpp::List) to the SPAmix constructor.
    // That constructor MUST run in the main thread (here), which it does.
    // So this is safe.
    
    const arma::mat t_resid           = Rcpp::as<arma::mat>(objNull["resid"]);
    const arma::mat t_PCs             = Rcpp::as<arma::mat>(objNull["PCs"]);
    const Rcpp::List t_outlierList    = Rcpp::as<Rcpp::List>(objNull["outLierList"]);

    int M = 0;
    int N = 0;
    std::string header;


    if (g_genoType == "PLINK") {

        if (ptr_gPLINKobj) delete ptr_gPLINKobj;

        ptr_gPLINKobj = new PLINK::PlinkClass(
            t_bimFile, t_famFile, t_bedFile, t_SampleInModel, t_AlleleOrder
        );
        
        M = ptr_gPLINKobj->getM();
        N = ptr_gPLINKobj->getN();
    }

    Rcpp::Rcout << "    Number of subjects: " << N << std::endl;
    Rcpp::Rcout << "    Number of markers: "  << M << std::endl;


    if (g_method == "SPAmix") {
        
        g_nPheno = t_resid.n_cols;
        header = "Pheno\tMarker\tInfo\tAltFreq\tAltCounts\tMissingRate";

        for (int p = 0; p < g_nPheno; ++p) {
            header += "\tPvalue" + std::to_string(p + 1) + "\tZscore" + std::to_string(p + 1);
        }
    }

    // Buffers
    GenoInputBuffer inputBuffer(inputBufferSize);
    ResultOutputBuffer outputBuffer(outputBufferSize);

    // Threads
    std::thread reader(readerThread, &inputBuffer, M);
    std::thread writer(writerThread, &outputBuffer, outputFile, header);
    
    std::vector<std::thread> workers;
    workers.reserve(nWorkers);

    std::vector<SPAmix::SPAmixClass*> threadObjPtrs;

    if (g_method == "SPAmix") {
        
        threadObjPtrs.reserve(nWorkers);

        for (int i = 0; i < nWorkers; ++i) {
            
            SPAmix::SPAmixClass* t_ptr_gSPAmixobj = 
                new SPAmix::SPAmixClass(
                    t_resid, 
                    t_PCs, 
                    N,
                    t_SPA_Cutoff,
                    t_outlierList
            );
            threadObjPtrs.push_back(t_ptr_gSPAmixobj);
            workers.emplace_back(workerSPAmix, &inputBuffer, &outputBuffer, t_ptr_gSPAmixobj);
        }
    } else {
        Rcpp::stop("Only SPAmix is currently supported in Main2 async mode.");
    }

    // Wait for worker completion
    reader.join();
    for (auto& w : workers) w.join();
    
    // Clean up thread-local objects (SPAmixClass instances)
    for (auto ptr : threadObjPtrs) {
        delete ptr;
    }

    outputBuffer.setFinished();
    writer.join();

    Rcpp::Rcout << "Analysis complete. Results written to: " << outputFile << std::endl;
}

