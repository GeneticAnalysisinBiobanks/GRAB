// mtMain.cpp -- Thread pool, per-chunk worker, sequential writer
//
// Functions by role:
//
//   [Inline helpers]
//     numToStr           — format a double as a short string ("NA", "Inf", or %.6g)
//     cctPvalue          — Cauchy combination test aggregate p-value
//     appendCell         — append a tab-separated cell to an ostringstream
//     appendMeta         — append the 8 standard marker meta-columns
//
//   [File-local structs]
//     ThreadContext      — per-worker copies of the active statistical method object
//
//   [Output header]       (static)
//     getHeader           — produce the TSV column header for a given method
//
//   [Worker context]      (static)
//     makeThreadContext   — clone the global method object into a per-thread ThreadContext
//
//   [Engine]              (mtMarker namespace, external linkage)
//     mainMarkerChunksCore — coordinate the PLINK cursor, worker pool, and writer thread

#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <exception>
#include <cstdio>
#include <RcppArmadillo.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <zlib.h>

#include "mtPLINK.h"
#include "mtPOLMM.h"
#include "mtWtCoxG.h"
#include "mtLEAF.h"
#include "mtSPAGRM.h"
#include "mtSAGELD.h"
#include "mtSPAsqr.h"
#include "mtSPACox.h"
#include "mtSPAmix.h"
#include "mtSPAmixPlus.h"


// ---- Global method pointers (one per active statistical method) ----

POLMM::POLMMClass*           ptr_gPOLMMobj      = nullptr;
WtCoxG::WtCoxGClass*         ptr_gWtCoxGobj     = nullptr;
LEAF::LEAFClass*             ptr_gLEAFobj       = nullptr;
SPAGRM::SPAGRMClass*         ptr_gSPAGRMobj     = nullptr;
SAGELD::SAGELDClass*         ptr_gSAGELDobj     = nullptr;
SPAsqr::SPAsqrClass*         ptr_gSPAsqrobj     = nullptr;
SPACox::SPACoxClass*         ptr_gSPACoxobj     = nullptr;
SPAmix::SPAmixClass*         ptr_gSPAmixobj     = nullptr;
SPAmixPlus::SPAmixPlusClass* ptr_gSPAmixPlusobj = nullptr;


namespace {

// ---- Inline helpers ----

// Format a double as a short string: "NA", "Inf", "-Inf", or %.6g.
std::string numToStr(double x) {
  if (std::isnan(x)) return "NA";
  if (std::isinf(x)) return (x > 0 ? "Inf" : "-Inf");
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.6g", x);
  return buf;
}

// Cauchy combination test: aggregate a vector of p-values into one p-value.
double cctPvalue(const std::vector<double>& pvals) {
  std::vector<double> p;
  p.reserve(pvals.size());
  for (double x : pvals)
    if (!std::isnan(x)) p.push_back(x);
  if (p.empty()) return std::numeric_limits<double>::quiet_NaN();
  double tStat = 0.0;
  for (double x : p) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) x = 0.999;
    tStat += std::tan((0.5 - x) * M_PI);
  }
  tStat /= static_cast<double>(p.size());
  if (tStat > 1e15) return (1.0 / tStat) / M_PI;
  return 0.5 - std::atan(tStat) / M_PI;
}

// Append one tab-separated value to an output row buffer.
void appendCell(std::ostringstream& oss, const std::string& v, bool first = false) {
  if (!first) oss << '\t';
  oss << v;
}
void appendCell(std::ostringstream& oss, double v, bool first = false) {
  appendCell(oss, numToStr(v), first);
}

// Append the 8 standard marker meta-columns (Marker, Chr, Pos, Ref, Alt,
// AltCount, AltFreq, MissRate) as the first cells of an output row.
void appendMeta(
  std::ostringstream& oss,
  const std::string& marker,
  const std::string& chr,
  uint32_t pos,
  const std::string& ref,
  const std::string& alt,
  double altCounts,
  double altFreq,
  double missingRate
) {
  appendCell(oss, marker, true);
  appendCell(oss, chr);
  appendCell(oss, std::to_string(pos));
  appendCell(oss, ref);
  appendCell(oss, alt);
  appendCell(oss, altCounts);
  appendCell(oss, altFreq);
  appendCell(oss, missingRate);
}

// ---- File-local structs ----

// Per-worker independent copies of the active method object.
struct ThreadContext {
  std::unique_ptr<POLMM::POLMMClass>           polmm;
  std::unique_ptr<SPACox::SPACoxClass>         spacox;
  std::unique_ptr<SPAmix::SPAmixClass>         spamix;
  std::unique_ptr<SPAmixPlus::SPAmixPlusClass> spamixPlus;
  std::unique_ptr<SPAGRM::SPAGRMClass>         spagrm;
  std::unique_ptr<SAGELD::SAGELDClass>         sageld;
  std::unique_ptr<WtCoxG::WtCoxGClass>         wtcoxg;
  std::unique_ptr<SPAsqr::SPAsqrClass>         spasqr;
  std::unique_ptr<LEAF::LEAFClass>             leaf;
};


// ---- Output header ----

const char* META_HEADER = "Marker\tChr\tPos\tRef\tAlt\tAltCount\tAltFreq\tMissRate";

// Build the TSV column header line for the given method.
std::string getHeader(const std::string& method,
                      const std::string& sageldMethod,
                      const std::vector<double>& taus,
                      int nPheno,
                      int nCluster) {
  std::ostringstream oss;
  oss << META_HEADER;

  if (method == "POLMM") {
    oss << "\tPvalue\tbeta\tseBeta\tzScore";
  } else if (method == "SPACox") {
    oss << "\tPvalue\tzScore";
  } else if (method == "SPAmix" || method == "SPAmixPlus") {
    oss << "\tPheno\tPvalue\tzScore";
  } else if (method == "SPAGRM") {
    oss << "\tzScore\tPvalue\thwepval";
  } else if (method == "SAGELD") {
    if (sageldMethod == "GALLOP")
      oss << "\tMethod\tBeta_G\tBeta_GxE\tSE_G\tSE_GxE\tPvalue_G\tPvalue_GxE\thwepval";
    else
      oss << "\tMethod\tzScore_G\tzScore_GxE\tPvalue_G\tPvalue_GxE\thwepval";
  } else if (method == "WtCoxG") {
    oss << "\tWtCoxG.ext\tWtCoxG.noext\tzscore.ext\tzscore.noext"
        << "\tAF_ref\tAN_ref\tpvalue_bat";
  } else if (method == "SPAsqr") {
    oss << "\thwepval";
    for (double tau : taus) oss << "\tZ_tau" << numToStr(tau);
    for (double tau : taus) oss << "\tP_tau" << numToStr(tau);
    oss << "\tP_CCT";
  } else if (method == "LEAF") {
    oss << "\tmeta.p_ext\tmeta.p_noext";
    for (int i = 1; i <= nCluster; ++i)
      oss << "\tcl" << i << ".p_ext\tcl" << i << ".p_noext\tcl" << i << ".p_batch";
  } else {
    throw std::runtime_error("Unsupported method in mtMarker: " + method);
  }
  return oss.str();
}


// ---- Worker context ----

// Clone the global method object into a per-thread ThreadContext.
// Each method class stores per-marker mutable state (result vectors, scratch),
// so per-thread copies are required for thread safety.  The copy happens once
// per thread, not per chunk, so the cost is amortized.
ThreadContext makeThreadContext(const std::string& method) {
  ThreadContext ctx;
  if (method == "POLMM") {
    if (!ptr_gPOLMMobj) throw std::runtime_error("POLMM object is not initialized.");
    ctx.polmm.reset(new POLMM::POLMMClass(*ptr_gPOLMMobj));
  } else if (method == "SPACox") {
    if (!ptr_gSPACoxobj) throw std::runtime_error("SPACox object is not initialized.");
    ctx.spacox.reset(new SPACox::SPACoxClass(*ptr_gSPACoxobj));
  } else if (method == "SPAmix") {
    if (!ptr_gSPAmixobj) throw std::runtime_error("SPAmix object is not initialized.");
    ctx.spamix.reset(new SPAmix::SPAmixClass(*ptr_gSPAmixobj));
  } else if (method == "SPAmixPlus") {
    if (!ptr_gSPAmixPlusobj) throw std::runtime_error("SPAmixPlus object is not initialized.");
    ctx.spamixPlus.reset(new SPAmixPlus::SPAmixPlusClass(*ptr_gSPAmixPlusobj));
  } else if (method == "SPAGRM") {
    if (!ptr_gSPAGRMobj) throw std::runtime_error("SPAGRM object is not initialized.");
    ctx.spagrm.reset(new SPAGRM::SPAGRMClass(*ptr_gSPAGRMobj));
  } else if (method == "SAGELD") {
    if (!ptr_gSAGELDobj) throw std::runtime_error("SAGELD object is not initialized.");
    ctx.sageld.reset(new SAGELD::SAGELDClass(*ptr_gSAGELDobj));
  } else if (method == "WtCoxG") {
    if (!ptr_gWtCoxGobj) throw std::runtime_error("WtCoxG object is not initialized.");
    ctx.wtcoxg.reset(new WtCoxG::WtCoxGClass(*ptr_gWtCoxGobj));
  } else if (method == "SPAsqr") {
    if (!ptr_gSPAsqrobj) throw std::runtime_error("SPAsqr object is not initialized.");
    ctx.spasqr.reset(new SPAsqr::SPAsqrClass(*ptr_gSPAsqrobj));
  } else if (method == "LEAF") {
    if (!ptr_gLEAFobj) throw std::runtime_error("LEAF object is not initialized.");
    ctx.leaf.reset(new LEAF::LEAFClass(*ptr_gLEAFobj));
  }
  return ctx;
}

} // namespace


// ---- The Engine ----

void mtMarkerEngine(
  const std::string method,
  const std::string bedFile,
  const std::string bimFile,
  const std::string famFile,
  const std::string outputFile,
  const std::vector<std::string>& subjData,
  const std::string AlleleOrder,
  const int nMarkersEachChunk,
  const int nthreads,
  const std::string impute_method,
  const double missing_cutoff,
  const double min_maf_marker,
  const double min_mac_marker,
  const std::string IDsToIncludeFile,
  const std::string RangesToIncludeFile,
  const std::string IDsToExcludeFile,
  const std::string RangesToExcludeFile
) {

  mtPLINK::PlinkData plinkData(bedFile, bimFile, famFile, subjData, AlleleOrder);

  mtPLINK::MarkerFilterConfig filterConfig;
  filterConfig.IDsToIncludeFile    = IDsToIncludeFile;
  filterConfig.RangesToIncludeFile = RangesToIncludeFile;
  filterConfig.IDsToExcludeFile    = IDsToExcludeFile;
  filterConfig.RangesToExcludeFile = RangesToExcludeFile;

  std::vector<mtPLINK::MarkerInfo> markerInfo = plinkData.getFilteredMarkers(filterConfig);
  if (markerInfo.empty()) {
    throw std::runtime_error("No markers remain after PLINK marker filtering.");
  }

  std::vector<std::vector<uint64_t>> chunkIndices = mtPLINK::PlinkData::buildChunks(markerInfo, nMarkersEachChunk);
  const int effective_nthreads = std::min(nthreads, static_cast<int>(chunkIndices.size()));

  Rprintf("Number of subjects in the input file: %u\n", plinkData.nSubjInFile());
  Rprintf("Number of subjects to test: %u\n", plinkData.nSubjUsed());
  Rprintf("Number of markers in the input file: %u\n", plinkData.nMarkers());
  Rprintf("Number of markers to test: %zu\n", markerInfo.size());
  Rprintf("Number of markers in each chunk: %d\n", nMarkersEachChunk);
  Rprintf("Number of chunks for all markers: %zu\n", chunkIndices.size());
  Rprintf("Number of threads: %d\n", effective_nthreads);

  int nPheno = 1;
  if (method == "SPAmix") nPheno = ptr_gSPAmixobj->getNpheno();
  if (method == "SPAmixPlus") nPheno = ptr_gSPAmixPlusobj->getNpheno();
  if (method == "SPAsqr") nPheno = ptr_gSPAsqrobj->get_ntaus();

  const std::string sageldMethod = (method == "SAGELD" && ptr_gSAGELDobj) ? ptr_gSAGELDobj->getMethod() : "SAGELD";
  const std::vector<double> spasqrTaus = (method == "SPAsqr" && ptr_gSPAsqrobj) ? ptr_gSPAsqrobj->getTaus() : std::vector<double>{};
  const int leafNcluster = (method == "LEAF" && ptr_gLEAFobj) ? ptr_gLEAFobj->get_Ncluster() : 0;

  const std::string header = getHeader(method, sageldMethod, spasqrTaus, nPheno, leafNcluster);

  std::vector<std::string> chunkOutput(chunkIndices.size());
  std::vector<char> chunkReady(chunkIndices.size(), 0);
  std::atomic<size_t> nextChunk(0);

  std::exception_ptr workerError = nullptr;
  std::mutex errorMutex;
  std::mutex writeMutex;
  std::mutex logMutex;
  std::condition_variable writeCv;
  bool stopWriter = false;


  std::thread writerThread([&]() {
    const bool useGzip = outputFile.size() >= 3 &&
      outputFile.compare(outputFile.size() - 3, 3, ".gz") == 0;

    gzFile gz = nullptr;
    std::ofstream plainOut;

    if (useGzip) {
      gz = gzopen(outputFile.c_str(), "wb");
      if (!gz) {
        std::lock_guard<std::mutex> lock(errorMutex);
        workerError = std::make_exception_ptr(std::runtime_error("Cannot open gzip output file: " + outputFile));
        stopWriter = true;
        writeCv.notify_all();
        return;
      }
    } else {
      plainOut.open(outputFile.c_str());
      if (!plainOut.is_open()) {
        std::lock_guard<std::mutex> lock(errorMutex);
        workerError = std::make_exception_ptr(std::runtime_error("Cannot open output file: " + outputFile));
        stopWriter = true;
        writeCv.notify_all();
        return;
      }
    }

    auto writeStr = [&](const std::string& s) {
      if (useGzip) gzwrite(gz, s.data(), static_cast<unsigned>(s.size()));
      else plainOut << s;
    };

    writeStr(header + "\n");

    for (size_t i = 0; i < chunkOutput.size(); ++i) {
      std::unique_lock<std::mutex> lk(writeMutex);
      writeCv.wait(lk, [&]() { return chunkReady[i] || stopWriter; });
      if (!chunkReady[i]) break;
      writeStr(chunkOutput[i]);
      {
        std::lock_guard<std::mutex> lg(logMutex);
        std::fprintf(stderr, "[INFO] Writing finished: chunk %zu/%zu\n", i + 1, chunkOutput.size());
        std::fflush(stderr);
      }
    }

    if (useGzip) gzclose(gz);
    else plainOut.close();
  });


  auto workerFn = [&]() {
    try {
      // Per-thread: copy method object (once) + open .bed file handle
      ThreadContext ctx = makeThreadContext(method);
      mtPLINK::PlinkCursor cursor(plinkData);

      while (true) {
        size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= chunkIndices.size()) break;

        const auto& gIdx = chunkIndices[cidx];
        std::ostringstream out;

        if (!gIdx.empty()) cursor.beginSequentialBlock(gIdx.front());

        // WtCoxG per-chunk prep: look up refMap and call updateMarkerInfo
        if (method == "WtCoxG") {
          ctx.wtcoxg->prepareChunk(gIdx);
        }

        // LEAF per-chunk prep: look up per-cluster refMaps
        if (method == "LEAF") {
          ctx.leaf->prepareChunk(gIdx);
        }

        for (size_t i = 0; i < gIdx.size(); ++i) {
          double altFreq = NAN, altCounts = NAN, missingRate = NAN;
          std::vector<uint32_t> indexForMissing;

          // Read genotype from per-thread cursor; metadata from shared plinkData
          arma::vec GVec = cursor.getGenotypes(gIdx[i],
            altFreq, altCounts, missingRate, indexForMissing);

          const std::string& chr    = plinkData.chr(gIdx[i]);
          const std::string& ref    = plinkData.ref(gIdx[i]);
          const std::string& alt    = plinkData.alt(gIdx[i]);
          const std::string& marker = plinkData.markerId(gIdx[i]);
          uint32_t pd               = plinkData.pos(gIdx[i]);

          const double n = static_cast<double>(GVec.size());
          const double maf = std::min(altFreq, 1.0 - altFreq);
          const double mac = 2.0 * maf * n * (1.0 - missingRate);
          bool passQC = !((missingRate > missing_cutoff) || (maf < min_maf_marker) || (mac < min_mac_marker));

          bool flip = false;
          if (passQC) {
            int nMissing = indexForMissing.size();
            double imputeG = 0;
            if (impute_method == "mean"){
              imputeG = 2 * altFreq;
            }

            if (impute_method == "minor"){
              if (altFreq > 0.5){
                imputeG = 2;
              }else{
                imputeG = 0;
              }
            }

            for (int j = 0; j < nMissing; j++){
              uint32_t index = indexForMissing.at(j);
              GVec.at(index) = imputeG;
            }

            if (method != "WtCoxG") {
              if (altFreq > 0.5){
                GVec = 2 - GVec;
                flip = true;
              }
            }
          }

          if (method == "SPAmix" || method == "SPAmixPlus") {
            std::vector<double> p(nPheno, NAN), z(nPheno, NAN);
            if (passQC) {
              if (method == "SPAmix") {
                ctx.spamix->getMarkerPval(GVec, altFreq);
                arma::vec pT = ctx.spamix->getpvalVec(), zT = ctx.spamix->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) { p[j] = pT[j]; z[j] = zT[j]; }
              } else {
                ctx.spamixPlus->getMarkerPval(GVec, altFreq);
                arma::vec pT = ctx.spamixPlus->getpvalVec(), zT = ctx.spamixPlus->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) { p[j] = pT[j]; z[j] = zT[j]; }
              }
            }
            for (int j = 0; j < nPheno; ++j) {
              std::ostringstream row;
              appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
              appendCell(row, "pheno_" + std::to_string(j + 1));
              appendCell(row, p[j]); appendCell(row, z[j]);
              out << row.str() << '\n';
            }
            continue;
          }


          if (method == "POLMM") {
            double Beta = NAN, seBeta = NAN, pval = NAN, zScore = NAN;
            if (passQC) {
              ctx.polmm->getMarkerPval(GVec, Beta, seBeta, pval, altFreq, zScore);
              Beta *= (1.0 - 2.0 * static_cast<double>(flip));
            }
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, pval);
            appendCell(row, Beta); appendCell(row, seBeta); appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SPACox") {
            double pval = NAN, zScore = NAN;
            if (passQC) pval = ctx.spacox->getMarkerPval(GVec, altFreq, zScore);
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, pval); appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SPAGRM") {
            double pval = NAN, zScore = NAN, hwepval = NAN;
            if (passQC) pval = ctx.spagrm->getMarkerPval(GVec, altFreq, zScore, hwepval);
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, zScore);
            appendCell(row, pval); appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SAGELD") {
            double hwepval = NAN;
            double zG = NAN, zGxE = NAN, pG = NAN, pGxE = NAN;
            double bG = NAN, bGxE = NAN, seG = NAN, seGxE = NAN;
            if (passQC) {
              ctx.sageld->getMarkerPval(GVec, altFreq, hwepval);
              arma::vec pT = ctx.sageld->getpvalVec(), zT = ctx.sageld->getzScoreVec();
              arma::vec bT = ctx.sageld->getBetaVec(), seT = ctx.sageld->getseBetaVec();
              zG = zT[0]; zGxE = zT[1]; pG = pT[0]; pGxE = pT[1];
              bG = bT[0] * (1.0 - 2.0 * static_cast<double>(flip));
              bGxE = bT[1] * (1.0 - 2.0 * static_cast<double>(flip));
              seG = seT[0]; seGxE = seT[1];
            }
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, sageldMethod);
            if (sageldMethod == "GALLOP") {
              appendCell(row, bG); appendCell(row, bGxE);
              appendCell(row, seG); appendCell(row, seGxE);
              appendCell(row, pG); appendCell(row, pGxE);
            } else {
              appendCell(row, zG); appendCell(row, zGxE);
              appendCell(row, pG); appendCell(row, pGxE);
            }
            appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }


          if (method == "WtCoxG") {
            double pExt = NAN, pNoext = NAN, zExt = NAN, zNoext = NAN;
            if (passQC) {
              arma::vec pT = ctx.wtcoxg->getpvalVec(GVec, static_cast<int>(i));
              arma::vec zT = ctx.wtcoxg->getZScoreVec();
              pExt = pT[0]; pNoext = pT[1]; zExt = zT[0]; zNoext = zT[1];
            }
            const auto& wr = ctx.wtcoxg->getChunkRefInfo(i);
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, pExt); appendCell(row, pNoext);
            appendCell(row, zExt); appendCell(row, zNoext);
            appendCell(row, wr.AF_ref); appendCell(row, wr.AN_ref);
            appendCell(row, wr.pvalue_bat);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SPAsqr") {
            double hwepval = NAN;
            std::vector<double> z(spasqrTaus.size(), NAN), p(spasqrTaus.size(), NAN);
            if (passQC) {
              arma::vec zT;
              arma::vec pT = ctx.spasqr->getMarkerPval(GVec, altFreq, zT, hwepval);
              for (size_t j = 0; j < p.size(); ++j) { p[j] = pT[j]; z[j] = zT[j]; }
            }
            double pCCT = cctPvalue(p);
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, hwepval);
            for (double v : z) appendCell(row, v);
            for (double v : p) appendCell(row, v);
            appendCell(row, pCCT);
            out << row.str() << '\n';
            continue;
          }


          if (method == "LEAF") {
            std::vector<double> pExt(leafNcluster, NAN), pNoext(leafNcluster, NAN);
            std::vector<double> sExt(leafNcluster, NAN), sNoext(leafNcluster, NAN);

            if (passQC) {
              arma::vec all = ctx.leaf->getMarkerZSP(GVec, static_cast<int>(i));
              const int nOut = 2 * leafNcluster;
              for (int c = 0; c < leafNcluster; ++c) {
                sExt[c] = all[nOut + 2*c]; sNoext[c] = all[nOut + 2*c + 1];
                pExt[c] = all[2*nOut + 2*c]; pNoext[c] = all[2*nOut + 2*c + 1];
              }
            }

            auto metaP = [&](const std::vector<double>& scores, const std::vector<double>& pvals) {
              double sumScore = 0, sumVar = 0;
              boost::math::chi_squared dist(1.0);
              for (size_t c = 0; c < scores.size(); ++c) {
                if (std::isnan(scores[c]) || std::isnan(pvals[c]) || pvals[c] <= 0 || pvals[c] >= 1) continue;
                double chisq = boost::math::quantile(dist, 1.0 - pvals[c]);
                if (chisq < 1e-30) chisq = 1e-30;
                double var = (scores[c] * scores[c]) / chisq;
                if (std::isnan(var)) var = 0;
                sumScore += scores[c]; sumVar += var;
              }
              if (sumVar <= 0) return std::numeric_limits<double>::quiet_NaN();
              double z = sumScore / std::sqrt(sumVar);
              return std::erfc(std::fabs(z) / std::sqrt(2.0));
            };

            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, metaP(sExt, pExt)); appendCell(row, metaP(sNoext, pNoext));
            for (int c = 0; c < leafNcluster; ++c) {
              appendCell(row, pExt[c]); appendCell(row, pNoext[c]);
              appendCell(row, ctx.leaf->getChunkRefInfo(c, i).pvalue_bat);
            }
            out << row.str() << '\n';
            continue;
          }
        }

        {
          std::lock_guard<std::mutex> lk(writeMutex);
          chunkOutput[cidx] = out.str();
          chunkReady[cidx] = 1;
        } {
          std::lock_guard<std::mutex> lg(logMutex);
          std::fprintf(stderr, "[INFO] Calculation finished: chunk %zu/%zu\n", cidx + 1, chunkIndices.size());
          std::fflush(stderr);
        }
        writeCv.notify_all();
      }
    } catch (...) {
      { std::lock_guard<std::mutex> lk(errorMutex); if (!workerError) workerError = std::current_exception(); }
      { std::lock_guard<std::mutex> lk(writeMutex); stopWriter = true; }
      writeCv.notify_all();
    }
  };


  const bool runInline = (effective_nthreads == 1);
  if (runInline) {
    workerFn();
  } else {
    std::vector<std::thread> workers;
    workers.reserve(effective_nthreads);
    for (int t = 0; t < effective_nthreads; ++t) workers.emplace_back(workerFn);
    for (auto& th : workers) th.join();
  }

  { std::lock_guard<std::mutex> lk(writeMutex); stopWriter = true; }
  writeCv.notify_all();
  writerThread.join();

  if (workerError) std::rethrow_exception(workerError);
}
