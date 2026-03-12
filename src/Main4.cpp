// Main4.cpp - chunk-level multi-threaded marker analysis and file writing

#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <limits>
#include <cmath>
#include <exception>

#include <boost/math/distributions/chi_squared.hpp>

#include "PLINK.h"
#include "BGEN.h"
#include "POLMM.h"
#include "SPACox.h"
#include "SPAmix.h"
#include "SPAGRM.h"
#include "SAGELD.h"
#include "WtCoxG.h"
#include "SPAsqr.h"
#include "LEAF.h"
#include "SPAmixPlus.h"
#include "UTIL.h"
#include "Main.h"

namespace {

std::string numToStr(double x) {
  if (std::isnan(x)) return "NA";
  if (std::isinf(x)) return (x > 0 ? "Inf" : "-Inf");
  std::ostringstream oss;
  oss << std::setprecision(15) << x;
  return oss.str();
}

void appendCell(std::ostringstream &oss, const std::string &v, bool first = false) {
  if (!first) oss << '\t';
  oss << v;
}

void appendCell(std::ostringstream &oss, double v, bool first = false) {
  appendCell(oss, numToStr(v), first);
}

double cctPvalue(const std::vector<double> &pvals) {
  std::vector<double> p;
  p.reserve(pvals.size());
  for (double x : pvals) {
    if (!std::isnan(x)) p.push_back(x);
  }
  if (p.empty()) return std::numeric_limits<double>::quiet_NaN();

  double tStat = 0.0;
  for (double x : p) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) x = 0.999;
    tStat += std::tan((0.5 - x) * M_PI);
  }
  tStat /= static_cast<double>(p.size());

  if (tStat > 1e15) {
    return (1.0 / tStat) / M_PI;
  }
  return 0.5 - std::atan(tStat) / M_PI;
}

std::string getHeader(const std::string &method,
                      const std::string &sageldMethod,
                      const std::vector<double> &taus,
                      int nPheno,
                      int nCluster) {
  std::ostringstream oss;

  if (method == "POLMM") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tPvalue\tbeta\tseBeta\tzScore";
  } else if (method == "SPACox") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tPvalue\tzScore";
  } else if (method == "SPAmix" || method == "SPAmixPlus") {
    oss << "Pheno\tMarker\tInfo\tAltFreq\tAltCounts\tMissingRate\tPvalue\tzScore";
  } else if (method == "SPAGRM") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tzScore\tPvalue\thwepval";
  } else if (method == "SAGELD") {
    if (sageldMethod == "GALLOP") {
      oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tMethod\tBeta_G\tBeta_GxE\tSE_G\tSE_GxE\tPvalue_G\tPvalue_GxE\thwepval";
    } else {
      oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tMethod\tzScore_G\tzScore_GxE\tPvalue_G\tPvalue_GxE\thwepval";
    }
  } else if (method == "WtCoxG") {
    oss << "WtCoxG.ext\tWtCoxG.noext\tscore.ext\tscore.noext\tzscore.ext\tzscore.noext"
        << "\tAF_ref\tAN_ref\tTPR\tsigma2\tpvalue_bat\tw.ext\tvar.ratio.w0\tvar.ratio.int\tvar.ratio.ext";
  } else if (method == "SPAsqr") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\thwepval";
    for (double tau : taus) {
      oss << "\tZ_tau" << numToStr(tau);
    }
    for (double tau : taus) {
      oss << "\tP_tau" << numToStr(tau);
    }
    oss << "\tP_CCT";
  } else if (method == "LEAF") {
    oss << "Marker\tInfo\tAltFreq\tMissingRate\tmeta.p_ext\tmeta.p_noext";
    for (int i = 1; i <= nCluster; ++i) {
      oss << "\tcl" << i << ".p_ext\tcl" << i << ".p_noext\tcl" << i << ".p_batch";
    }
  } else {
    Rcpp::stop("Unsupported method in mainMarkerChunksInCPP4: " + method);
  }

  return oss.str();
}

struct WtRow {
  double AF_ref;
  double AN_ref;
  double TPR;
  double sigma2;
  double pvalue_bat;
  double w_ext;
  double var_ratio_w0;
  double var_ratio_int;
  double var_ratio_ext;
};

struct LeafRow {
  double AF_ref;
  double AN_ref;
  double TPR;
  double sigma2;
  double pvalue_bat;
  double w_ext;
  double var_ratio_w0;
  double var_ratio_int;
  double var_ratio_ext;
};

struct ThreadContext {
  std::unique_ptr<POLMM::POLMMClass> polmm;
  std::unique_ptr<SPACox::SPACoxClass> spacox;
  std::unique_ptr<SPAmix::SPAmixClass> spamix;
  std::unique_ptr<SPAmixPlus::SPAmixPlusClass> spamixPlus;
  std::unique_ptr<SPAGRM::SPAGRMClass> spagrm;
  std::unique_ptr<SAGELD::SAGELDClass> sageld;
  std::unique_ptr<WtCoxG::WtCoxGClass> wtcoxg;
  std::unique_ptr<SPAsqr::SPAsqrClass> spasqr;
  std::unique_ptr<LEAF::LEAFClass> leaf;
};

ThreadContext makeThreadContext(const std::string &method) {
  ThreadContext ctx;
  // Each worker gets its own method object copy to avoid shared mutable state.
  if (method == "POLMM") {
    if (!ptr_gPOLMMobj) Rcpp::stop("POLMM object is not initialized.");
    ctx.polmm.reset(new POLMM::POLMMClass(*ptr_gPOLMMobj));
  } else if (method == "SPACox") {
    if (!ptr_gSPACoxobj) Rcpp::stop("SPACox object is not initialized.");
    ctx.spacox.reset(new SPACox::SPACoxClass(*ptr_gSPACoxobj));
  } else if (method == "SPAmix") {
    if (!ptr_gSPAmixobj) Rcpp::stop("SPAmix object is not initialized.");
    ctx.spamix.reset(new SPAmix::SPAmixClass(*ptr_gSPAmixobj));
  } else if (method == "SPAmixPlus") {
    if (!ptr_gSPAmixPlusobj) Rcpp::stop("SPAmixPlus object is not initialized.");
    ctx.spamixPlus.reset(new SPAmixPlus::SPAmixPlusClass(*ptr_gSPAmixPlusobj));
  } else if (method == "SPAGRM") {
    if (!ptr_gSPAGRMobj) Rcpp::stop("SPAGRM object is not initialized.");
    ctx.spagrm.reset(new SPAGRM::SPAGRMClass(*ptr_gSPAGRMobj));
  } else if (method == "SAGELD") {
    if (!ptr_gSAGELDobj) Rcpp::stop("SAGELD object is not initialized.");
    ctx.sageld.reset(new SAGELD::SAGELDClass(*ptr_gSAGELDobj));
  } else if (method == "WtCoxG") {
    if (!ptr_gWtCoxGobj) Rcpp::stop("WtCoxG object is not initialized.");
    ctx.wtcoxg.reset(new WtCoxG::WtCoxGClass(*ptr_gWtCoxGobj));
  } else if (method == "SPAsqr") {
    if (!ptr_gSPAsqrobj) Rcpp::stop("SPAsqr object is not initialized.");
    ctx.spasqr.reset(new SPAsqr::SPAsqrClass(*ptr_gSPAsqrobj));
  } else if (method == "LEAF") {
    if (!ptr_gLEAFobj) Rcpp::stop("LEAF object is not initialized.");
    ctx.leaf.reset(new LEAF::LEAFClass(*ptr_gLEAFobj));
  }
  return ctx;
}

} // namespace

// [[Rcpp::export]]
void mainMarkerChunksInCPP4(
  const std::string t_method,
  const Rcpp::List t_chunkIndexList,
  const std::string t_outputFile,
  const unsigned int t_nThreads,
  const std::string t_impute_method,
  const double t_missing_cutoff,
  const double t_min_maf_marker,
  const double t_min_mac_marker,
  const Rcpp::Nullable<Rcpp::List> t_extraParams = R_NilValue
) {
  if (t_chunkIndexList.size() == 0) {
    Rcpp::stop("t_chunkIndexList is empty.");
  }
  if (t_nThreads < 1) {
    Rcpp::stop("t_nThreads should be >= 1.");
  }

  // Convert chunk markers to plain C++ storage once (thread-safe read later).
  std::vector<std::vector<uint64_t>> chunkMarkers(t_chunkIndexList.size());
  for (int i = 0; i < t_chunkIndexList.size(); ++i) {
    std::vector<int64_t> idx = Rcpp::as<std::vector<int64_t>>(t_chunkIndexList[i]);
    chunkMarkers[i].reserve(idx.size());
    for (int64_t v : idx) {
      chunkMarkers[i].push_back(static_cast<uint64_t>(v));
    }
  }

  std::string sageldMethod = "SAGELD";
  std::vector<double> spasqrTaus;
  std::unordered_map<uint64_t, WtRow> wtMap;
  std::vector<std::unordered_map<uint64_t, LeafRow>> leafMaps;
  int leafNcluster = 0;

  std::string readerGenoType;
  std::string plinkBimFile, plinkFamFile, plinkBedFile;
  std::vector<std::string> readerSampleInModel;
  std::string bgenFile, bgenIndexFile;
  std::vector<std::string> bgenSampleInBgen;
  bool bgenSparseDosage = false;
  bool bgenDropMissing = false;
  std::string readerAlleleOrder;

  if (!t_extraParams.isNull()) {
    Rcpp::List extra(t_extraParams);

    if (extra.containsElementNamed("reader_config") && !Rcpp::as<Rcpp::RObject>(extra["reader_config"]).isNULL()) {
      Rcpp::List rc = Rcpp::as<Rcpp::List>(extra["reader_config"]);
      readerGenoType = Rcpp::as<std::string>(rc["genoType"]);
      readerAlleleOrder = Rcpp::as<std::string>(rc["alleleOrder"]);
      readerSampleInModel = Rcpp::as<std::vector<std::string>>(rc["sampleInModel"]);

      if (readerGenoType == "PLINK") {
        plinkBimFile = Rcpp::as<std::string>(rc["bimFile"]);
        plinkFamFile = Rcpp::as<std::string>(rc["famFile"]);
        plinkBedFile = Rcpp::as<std::string>(rc["bedFile"]);
      } else if (readerGenoType == "BGEN") {
        bgenFile = Rcpp::as<std::string>(rc["bgenFile"]);
        bgenIndexFile = Rcpp::as<std::string>(rc["bgiFile"]);
        bgenSampleInBgen = Rcpp::as<std::vector<std::string>>(rc["sampleInBgen"]);
        if (rc.containsElementNamed("isSparseDosageInBgen")) {
          bgenSparseDosage = Rcpp::as<bool>(rc["isSparseDosageInBgen"]);
        }
        if (rc.containsElementNamed("isDropmissingdosagesInBgen")) {
          bgenDropMissing = Rcpp::as<bool>(rc["isDropmissingdosagesInBgen"]);
        }
      }
    }

    if (extra.containsElementNamed("sageld_method") && !Rcpp::as<Rcpp::RObject>(extra["sageld_method"]).isNULL()) {
      sageldMethod = Rcpp::as<std::string>(extra["sageld_method"]);
    }

    if (extra.containsElementNamed("spasqr_taus") && !Rcpp::as<Rcpp::RObject>(extra["spasqr_taus"]).isNULL()) {
      spasqrTaus = Rcpp::as<std::vector<double>>(extra["spasqr_taus"]);
    }

    if (t_method == "WtCoxG" && extra.containsElementNamed("wtcoxg_merge") && !Rcpp::as<Rcpp::RObject>(extra["wtcoxg_merge"]).isNULL()) {
      Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(extra["wtcoxg_merge"]);
      std::vector<int64_t> genoIndex = Rcpp::as<std::vector<int64_t>>(df["genoIndex"]);
      std::vector<double> AF_ref = Rcpp::as<std::vector<double>>(df["AF_ref"]);
      std::vector<double> AN_ref = Rcpp::as<std::vector<double>>(df["AN_ref"]);
      std::vector<double> TPR = Rcpp::as<std::vector<double>>(df["TPR"]);
      std::vector<double> sigma2 = Rcpp::as<std::vector<double>>(df["sigma2"]);
      std::vector<double> pvalue_bat = Rcpp::as<std::vector<double>>(df["pvalue_bat"]);
      std::vector<double> w_ext = Rcpp::as<std::vector<double>>(df["w.ext"]);
      std::vector<double> var_ratio_w0 = Rcpp::as<std::vector<double>>(df["var.ratio.w0"]);
      std::vector<double> var_ratio_int = Rcpp::as<std::vector<double>>(df["var.ratio.int"]);
      std::vector<double> var_ratio_ext = Rcpp::as<std::vector<double>>(df["var.ratio.ext"]);

      for (size_t i = 0; i < genoIndex.size(); ++i) {
        wtMap[static_cast<uint64_t>(genoIndex[i])] = WtRow{
          AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
          w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
        };
      }
    }

    if (t_method == "LEAF" && extra.containsElementNamed("leaf_subgeno") && !Rcpp::as<Rcpp::RObject>(extra["leaf_subgeno"]).isNULL()) {
      Rcpp::List subList = Rcpp::as<Rcpp::List>(extra["leaf_subgeno"]);
      leafNcluster = subList.size();
      leafMaps.resize(leafNcluster);

      for (int c = 0; c < leafNcluster; ++c) {
        Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(subList[c]);
        std::vector<int64_t> genoIndex = Rcpp::as<std::vector<int64_t>>(df["genoIndex"]);
        std::vector<double> AF_ref = Rcpp::as<std::vector<double>>(df["AF_ref"]);
        std::vector<double> AN_ref = Rcpp::as<std::vector<double>>(df["AN_ref"]);
        std::vector<double> TPR = Rcpp::as<std::vector<double>>(df["TPR"]);
        std::vector<double> sigma2 = Rcpp::as<std::vector<double>>(df["sigma2"]);
        std::vector<double> pvalue_bat = Rcpp::as<std::vector<double>>(df["pvalue_bat"]);
        std::vector<double> w_ext = Rcpp::as<std::vector<double>>(df["w.ext"]);
        std::vector<double> var_ratio_w0 = Rcpp::as<std::vector<double>>(df["var.ratio.w0"]);
        std::vector<double> var_ratio_int = Rcpp::as<std::vector<double>>(df["var.ratio.int"]);
        std::vector<double> var_ratio_ext = Rcpp::as<std::vector<double>>(df["var.ratio.ext"]);

        for (size_t i = 0; i < genoIndex.size(); ++i) {
          leafMaps[c][static_cast<uint64_t>(genoIndex[i])] = LeafRow{
            AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
            w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
          };
        }
      }
    }

    if (extra.containsElementNamed("leaf_ncluster") && !Rcpp::as<Rcpp::RObject>(extra["leaf_ncluster"]).isNULL()) {
      leafNcluster = Rcpp::as<int>(extra["leaf_ncluster"]);
    }
  }

  if (readerGenoType.empty()) {
    Rcpp::stop("reader_config is required in t_extraParams for mainMarkerChunksInCPP4.");
  }

  int nPheno = 1;
  if (t_method == "SPAmix") nPheno = ptr_gSPAmixobj->getNpheno();
  if (t_method == "SPAmixPlus") nPheno = ptr_gSPAmixPlusobj->getNpheno();
  if (t_method == "SPAsqr") nPheno = ptr_gSPAsqrobj->get_ntaus();

  const std::string header = getHeader(t_method, sageldMethod, spasqrTaus, nPheno, leafNcluster);

  std::vector<std::string> chunkOutput(chunkMarkers.size());
  std::vector<char> chunkReady(chunkMarkers.size(), 0);
  std::atomic<size_t> nextChunk(0);

  std::exception_ptr workerError = nullptr;
  std::mutex errorMutex;
  std::mutex writeMutex;
  std::condition_variable writeCv;
  bool stopWriter = false;

  std::thread writerThread([&]() {
    std::ofstream out(t_outputFile.c_str());
    if (!out.is_open()) {
      std::lock_guard<std::mutex> lock(errorMutex);
      workerError = std::make_exception_ptr(std::runtime_error("Cannot open output file: " + t_outputFile));
      stopWriter = true;
      writeCv.notify_all();
      return;
    }

    out << header << '\n';

    // Preserve deterministic output order by flushing chunks strictly by index.
    for (size_t i = 0; i < chunkOutput.size(); ++i) {
      std::unique_lock<std::mutex> lk(writeMutex);
      writeCv.wait(lk, [&]() { return chunkReady[i] || stopWriter; });
      if (!chunkReady[i]) {
        break;
      }
      out << chunkOutput[i];
      Rcpp::Rcout << "[INFO] Completed chunks: " << (i + 1) << "/" << chunkOutput.size() << std::endl;
    }
    out.close();
  });

  auto workerFn = [&]() {
    try {
      ThreadContext ctx = makeThreadContext(t_method);

      std::unique_ptr<PLINK::PlinkClass> plinkReader;
      std::unique_ptr<BGEN::BgenClass> bgenReader;

      if (readerGenoType == "PLINK") {
        plinkReader.reset(new PLINK::PlinkClass(
          plinkBimFile,
          plinkFamFile,
          plinkBedFile,
          readerSampleInModel,
          readerAlleleOrder
        ));
      } else if (readerGenoType == "BGEN") {
        bgenReader.reset(new BGEN::BgenClass(
          bgenFile,
          bgenIndexFile,
          bgenSampleInBgen,
          readerSampleInModel,
          bgenSparseDosage,
          bgenDropMissing,
          readerAlleleOrder
        ));
      } else {
        Rcpp::stop("Unsupported reader genoType in reader_config: " + readerGenoType);
      }

      while (true) {
        // Dynamic scheduling keeps workers busy when chunk runtimes are uneven.
        size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= chunkMarkers.size()) break;

        const auto &gIdx = chunkMarkers[cidx];
        std::ostringstream out;

        std::vector<WtRow> wtChunkRows;
        std::vector<std::vector<LeafRow>> leafChunkRows;

        if (t_method == "WtCoxG") {
          wtChunkRows.reserve(gIdx.size());
          std::vector<double> AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext;
          AF_ref.reserve(gIdx.size());
          AN_ref.reserve(gIdx.size());
          TPR.reserve(gIdx.size());
          sigma2.reserve(gIdx.size());
          pvalue_bat.reserve(gIdx.size());
          w_ext.reserve(gIdx.size());
          var_w0.reserve(gIdx.size());
          var_int.reserve(gIdx.size());
          var_ext.reserve(gIdx.size());

          for (auto gi : gIdx) {
            auto it = wtMap.find(gi);
            if (it == wtMap.end()) {
              wtChunkRows.push_back(WtRow{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
            } else {
              wtChunkRows.push_back(it->second);
            }
          }

          for (const auto &r : wtChunkRows) {
            AF_ref.push_back(r.AF_ref);
            AN_ref.push_back(r.AN_ref);
            TPR.push_back(r.TPR);
            sigma2.push_back(r.sigma2);
            pvalue_bat.push_back(r.pvalue_bat);
            w_ext.push_back(r.w_ext);
            var_w0.push_back(r.var_ratio_w0);
            var_int.push_back(r.var_ratio_int);
            var_ext.push_back(r.var_ratio_ext);
          }

          Rcpp::DataFrame chunkDf = Rcpp::DataFrame::create(
            Rcpp::Named("AF_ref") = AF_ref,
            Rcpp::Named("AN_ref") = AN_ref,
            Rcpp::Named("TPR") = TPR,
            Rcpp::Named("sigma2") = sigma2,
            Rcpp::Named("pvalue_bat") = pvalue_bat,
            Rcpp::Named("w.ext") = w_ext,
            Rcpp::Named("var.ratio.w0") = var_w0,
            Rcpp::Named("var.ratio.int") = var_int,
            Rcpp::Named("var.ratio.ext") = var_ext
          );
          ctx.wtcoxg->updateMarkerInfo(chunkDf);
        }

        if (t_method == "LEAF") {
          leafChunkRows.resize(leafNcluster);
          for (int c = 0; c < leafNcluster; ++c) {
            auto &rows = leafChunkRows[c];
            rows.reserve(gIdx.size());

            std::vector<double> AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext;
            AF_ref.reserve(gIdx.size());
            AN_ref.reserve(gIdx.size());
            TPR.reserve(gIdx.size());
            sigma2.reserve(gIdx.size());
            pvalue_bat.reserve(gIdx.size());
            w_ext.reserve(gIdx.size());
            var_w0.reserve(gIdx.size());
            var_int.reserve(gIdx.size());
            var_ext.reserve(gIdx.size());

            for (auto gi : gIdx) {
              auto it = leafMaps[c].find(gi);
              if (it == leafMaps[c].end()) {
                rows.push_back(LeafRow{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
              } else {
                rows.push_back(it->second);
              }
            }

            for (const auto &r : rows) {
              AF_ref.push_back(r.AF_ref);
              AN_ref.push_back(r.AN_ref);
              TPR.push_back(r.TPR);
              sigma2.push_back(r.sigma2);
              pvalue_bat.push_back(r.pvalue_bat);
              w_ext.push_back(r.w_ext);
              var_w0.push_back(r.var_ratio_w0);
              var_int.push_back(r.var_ratio_int);
              var_ext.push_back(r.var_ratio_ext);
            }

            Rcpp::DataFrame chunkDf = Rcpp::DataFrame::create(
              Rcpp::Named("AF_ref") = AF_ref,
              Rcpp::Named("AN_ref") = AN_ref,
              Rcpp::Named("TPR") = TPR,
              Rcpp::Named("sigma2") = sigma2,
              Rcpp::Named("pvalue_bat") = pvalue_bat,
              Rcpp::Named("w.ext") = w_ext,
              Rcpp::Named("var.ratio.w0") = var_w0,
              Rcpp::Named("var.ratio.int") = var_int,
              Rcpp::Named("var.ratio.ext") = var_ext
            );
            ctx.leaf->updateMarkerInfo(c, chunkDf);
          }
        }

        for (size_t i = 0; i < gIdx.size(); ++i) {
          double altFreq = NAN, altCounts = NAN, missingRate = NAN, imputeInfo = NAN;
          std::vector<uint32_t> indexForMissing, indexForNonZero;
          std::string chr, ref, alt, marker;
          uint32_t pd = 0;

          arma::vec GVec;
          if (readerGenoType == "PLINK") {
            GVec = plinkReader->getOneMarker(
              gIdx[i],
              ref,
              alt,
              marker,
              pd,
              chr,
              altFreq,
              altCounts,
              missingRate,
              imputeInfo,
              true,
              indexForMissing,
              false,
              indexForNonZero,
              true
            );
          } else {
            bool isBoolRead = false;
            GVec = bgenReader->getOneMarker(
              gIdx[i],
              ref,
              alt,
              marker,
              pd,
              chr,
              altFreq,
              altCounts,
              missingRate,
              imputeInfo,
              true,
              indexForMissing,
              false,
              indexForNonZero,
              isBoolRead
            );
          }

          const std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;
          const double n = static_cast<double>(GVec.size());
          const double maf = std::min(altFreq, 1.0 - altFreq);
          const double mac = 2.0 * maf * n * (1.0 - missingRate);

          // Apply marker-level missing/MAF/MAC filters before method-specific work.
          bool passQC = !(
            (missingRate > t_missing_cutoff) ||
            (maf < t_min_maf_marker) ||
            (mac < t_min_mac_marker)
          );

          // Imputation may also flip allele coding to keep ALT as the minor allele.
          bool flip = false;
          if (passQC) {
            flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, t_impute_method, t_method);
          }

          if (t_method == "SPAmix" || t_method == "SPAmixPlus") {
            std::vector<double> p(nPheno, NAN), z(nPheno, NAN);
            if (passQC) {
              if (t_method == "SPAmix") {
                ctx.spamix->getMarkerPval(GVec, altFreq);
                arma::vec pTmp = ctx.spamix->getpvalVec();
                arma::vec zTmp = ctx.spamix->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) {
                  p[j] = pTmp[j];
                  z[j] = zTmp[j];
                }
              } else {
                ctx.spamixPlus->getMarkerPval(GVec, altFreq);
                arma::vec pTmp = ctx.spamixPlus->getpvalVec();
                arma::vec zTmp = ctx.spamixPlus->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) {
                  p[j] = pTmp[j];
                  z[j] = zTmp[j];
                }
              }
            }

            for (int j = 0; j < nPheno; ++j) {
              std::ostringstream row;
              appendCell(row, "pheno_" + std::to_string(j + 1), true);
              appendCell(row, marker);
              appendCell(row, info);
              appendCell(row, altFreq);
              appendCell(row, altCounts);
              appendCell(row, missingRate);
              appendCell(row, p[j]);
              appendCell(row, z[j]);
              out << row.str() << '\n';
            }
            continue;
          }

          if (t_method == "POLMM") {
            double Beta = NAN, seBeta = NAN, pval = NAN, zScore = NAN;
            if (passQC) {
              ctx.polmm->getMarkerPval(GVec, Beta, seBeta, pval, altFreq, zScore);
              Beta = Beta * (1.0 - 2.0 * static_cast<double>(flip));
            }
            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, pval);
            appendCell(row, Beta);
            appendCell(row, seBeta);
            appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }

          if (t_method == "SPACox") {
            double pval = NAN, zScore = NAN;
            if (passQC) {
              pval = ctx.spacox->getMarkerPval(GVec, altFreq, zScore);
            }
            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, pval);
            appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }

          if (t_method == "SPAGRM") {
            double pval = NAN, zScore = NAN, hwepval = NAN;
            if (passQC) {
              pval = ctx.spagrm->getMarkerPval(GVec, altFreq, zScore, hwepval);
            }
            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, zScore);
            appendCell(row, pval);
            appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }

          if (t_method == "SAGELD") {
            double hwepval = NAN;
            double zG = NAN, zGxE = NAN, pG = NAN, pGxE = NAN;
            double bG = NAN, bGxE = NAN, seG = NAN, seGxE = NAN;
            if (passQC) {
              ctx.sageld->getMarkerPval(GVec, altFreq, hwepval);
              arma::vec pTmp = ctx.sageld->getpvalVec();
              arma::vec zTmp = ctx.sageld->getzScoreVec();
              arma::vec bTmp = ctx.sageld->getBetaVec();
              arma::vec seTmp = ctx.sageld->getseBetaVec();
              zG = zTmp[0];
              zGxE = zTmp[1];
              pG = pTmp[0];
              pGxE = pTmp[1];
              bG = bTmp[0] * (1.0 - 2.0 * static_cast<double>(flip));
              bGxE = bTmp[1] * (1.0 - 2.0 * static_cast<double>(flip));
              seG = seTmp[0];
              seGxE = seTmp[1];
            }

            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, sageldMethod);
            if (sageldMethod == "GALLOP") {
              appendCell(row, bG);
              appendCell(row, bGxE);
              appendCell(row, seG);
              appendCell(row, seGxE);
              appendCell(row, pG);
              appendCell(row, pGxE);
            } else {
              appendCell(row, zG);
              appendCell(row, zGxE);
              appendCell(row, pG);
              appendCell(row, pGxE);
            }
            appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }

          if (t_method == "WtCoxG") {
            double pExt = NAN, pNoext = NAN, sExt = NAN, sNoext = NAN, zExt = NAN, zNoext = NAN;
            if (passQC) {
              arma::vec pTmp = ctx.wtcoxg->getpvalVec(GVec, static_cast<int>(i));
              arma::vec sTmp = ctx.wtcoxg->getScoreVec();
              arma::vec zTmp = ctx.wtcoxg->getZScoreVec();
              pExt = pTmp[0];
              pNoext = pTmp[1];
              sExt = sTmp[0];
              sNoext = sTmp[1];
              zExt = zTmp[0];
              zNoext = zTmp[1];
            }

            const WtRow wr = wtChunkRows[i];
            std::ostringstream row;
            appendCell(row, pExt, true);
            appendCell(row, pNoext);
            appendCell(row, sExt);
            appendCell(row, sNoext);
            appendCell(row, zExt);
            appendCell(row, zNoext);
            appendCell(row, wr.AF_ref);
            appendCell(row, wr.AN_ref);
            appendCell(row, wr.TPR);
            appendCell(row, wr.sigma2);
            appendCell(row, wr.pvalue_bat);
            appendCell(row, wr.w_ext);
            appendCell(row, wr.var_ratio_w0);
            appendCell(row, wr.var_ratio_int);
            appendCell(row, wr.var_ratio_ext);
            out << row.str() << '\n';
            continue;
          }

          if (t_method == "SPAsqr") {
            double hwepval = NAN;
            std::vector<double> z(spasqrTaus.size(), NAN), p(spasqrTaus.size(), NAN);
            if (passQC) {
              arma::vec zTmp;
              arma::vec pTmp = ctx.spasqr->getMarkerPval(GVec, altFreq, zTmp, hwepval);
              for (size_t j = 0; j < p.size(); ++j) {
                p[j] = pTmp[j];
                z[j] = zTmp[j];
              }
            }
            double pCCT = cctPvalue(p);

            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, hwepval);
            for (double v : z) appendCell(row, v);
            for (double v : p) appendCell(row, v);
            appendCell(row, pCCT);
            out << row.str() << '\n';
            continue;
          }

          if (t_method == "LEAF") {
            std::vector<double> pExt(leafNcluster, NAN), pNoext(leafNcluster, NAN);
            std::vector<double> sExt(leafNcluster, NAN), sNoext(leafNcluster, NAN);

            if (passQC) {
              arma::vec all = ctx.leaf->getMarkerZSP(GVec, static_cast<int>(i));
              const int nOut = 2 * leafNcluster;
              for (int c = 0; c < leafNcluster; ++c) {
                int extIdx = 2 * c;
                int noextIdx = 2 * c + 1;
                sExt[c] = all[nOut + extIdx];
                sNoext[c] = all[nOut + noextIdx];
                pExt[c] = all[2 * nOut + extIdx];
                pNoext[c] = all[2 * nOut + noextIdx];
              }
            }

            auto metaP = [&](const std::vector<double> &scores, const std::vector<double> &pvals) {
              double sumScore = 0.0;
              double sumVar = 0.0;
              boost::math::chi_squared dist(1.0);

              for (size_t c = 0; c < scores.size(); ++c) {
                double sc = scores[c];
                double pv = pvals[c];
                if (std::isnan(sc) || std::isnan(pv) || pv <= 0.0 || pv >= 1.0) continue;
                double chisq = boost::math::quantile(dist, 1.0 - pv);
                if (chisq < 1e-30) chisq = 1e-30;
                double var = (sc * sc) / chisq;
                if (std::isnan(var)) var = 0.0;
                sumScore += sc;
                sumVar += var;
              }

              if (sumVar <= 0.0) return std::numeric_limits<double>::quiet_NaN();
              double z = sumScore / std::sqrt(sumVar);
              return std::erfc(std::fabs(z) / std::sqrt(2.0));
            };

            double metaExt = metaP(sExt, pExt);
            double metaNoext = metaP(sNoext, pNoext);

            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, missingRate);
            appendCell(row, metaExt);
            appendCell(row, metaNoext);

            for (int c = 0; c < leafNcluster; ++c) {
              appendCell(row, pExt[c]);
              appendCell(row, pNoext[c]);
              appendCell(row, leafChunkRows[c][i].pvalue_bat);
            }
            out << row.str() << '\n';
            continue;
          }
        }

        {
          // Hand completed chunk text to the single writer thread.
          std::lock_guard<std::mutex> lk(writeMutex);
          chunkOutput[cidx] = out.str();
          chunkReady[cidx] = 1;
        }
        writeCv.notify_all();
      }

    } catch (...) {
      {
        std::lock_guard<std::mutex> lock(errorMutex);
        if (!workerError) workerError = std::current_exception();
      }
      {
        std::lock_guard<std::mutex> lk(writeMutex);
        stopWriter = true;
      }
      writeCv.notify_all();
    }
  };

  const unsigned int nThreads = std::min<unsigned int>(t_nThreads, static_cast<unsigned int>(chunkMarkers.size()));
  std::vector<std::thread> workers;
  workers.reserve(nThreads);
  for (unsigned int t = 0; t < nThreads; ++t) {
    workers.emplace_back(workerFn);
  }
  for (auto &th : workers) {
    th.join();
  }

  {
    std::lock_guard<std::mutex> lk(writeMutex);
    stopWriter = true;
  }
  writeCv.notify_all();
  writerThread.join();

  if (workerError) {
    std::rethrow_exception(workerError);
  }
}
