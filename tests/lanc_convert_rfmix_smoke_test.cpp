// lanc_convert_rfmix_smoke_test.cpp — standalone smoke test for P3's
// convertRfmixToLanc (src/localplus/lanc_convert_rfmix.{hpp,cpp}).
//
// Hand-writes a tiny plain-text VCF (5 samples, chr20, 2 windows, one
// biallelic pair per window plus one multiallelic site that must be
// skipped) and a matching tiny RFMix .msp.tsv (K=3, 2 windows, including
// an unassigned (-1) call and a missing genotype to exercise both the
// ancestry-unassigned (0xF) path and the allele-missing plane), runs
// convertRfmixToLanc, then re-opens the output via LancData/LancCursor and
// asserts the decoded dosage_k/hapcount_k against hand-computed values for
// every (marker, sample, k) cell. A second pass exercises --keep.
//
// A third, focused case (runBoundaryCase) regresses the 0-based [spos, epos)
// window-membership fix (dev-notes/methods/recode_rfmix/03_merged_and_fixes.md
// item #4): a synthetic msp of 3 contiguous windows tiling half-open in
// 1-based RFMix coordinates, with query SNVs placed EXACTLY on window spos
// boundaries (including POS==spos0 of the first window). It asserts (a) the
// first-window first SNP is not dropped, (b) each boundary SNP lands in the
// LATER window that 1-based [spos, epos) contains it, and (c) the per-window
// (and total) assigned SNV counts equal RFMix's own `n snps` column.
//
// Build (from repo root, after `make -j`):
//   g++ -std=c++17 -O2 -Isrc -Ithird_party/eigen-5.0.0 \
//       -Ithird_party/zstd-1.5.7/lib -Ithird_party/htslib-1.23.1 \
//       tests/lanc_convert_rfmix_smoke_test.cpp \
//       build/localplus/lanc_io.o build/localplus/lanc_convert_rfmix.o \
//       build/geno_factory/variant_filter.o build/io/subject_filter.o \
//       build/htslib/*.o build/htslib/htscodecs/*.o \
//       build/zstd/common/*.o build/zstd/compress/*.o build/zstd/decompress/*.o \
//       build/zlib/*.o build/libdeflate/*.o build/libdeflate/x86/*.o build/libdeflate/arm/*.o \
//       -lpthread -lstdc++fs -o <bin>
//
// Prints a PASS/FAIL summary; exits non-zero on any failure.

#include "localplus/lanc_convert_rfmix.hpp"
#include "localplus/lanc_io.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static int g_failures = 0;
static int g_checks = 0;

static void check(bool cond, const std::string &msg) {
    ++g_checks;
    if (!cond) {
        ++g_failures;
        std::fprintf(stderr, "  FAIL: %s\n", msg.c_str());
    }
}

static bool feq(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    return a == b;
}

static void writeFile(const std::string &path, const std::string &content) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot create fixture file: " + path);
    f << content;
}

// ── Ground truth: K=3, N=5 samples (S1..S5), chr20, two windows ─────────
//
// RFMix .msp spos/epos are 1-based; convertRfmixToLanc stores them 0-based
// (minus one) so membership is the half-open 0-based [spos, epos):
//   Window1: 1-based [1,501)    -> 0-based [0,500):  calls (S1..S5, hap0/hap1)
//     S1:(0,0) S2:(1,2) S3:(2,2) S4:(0,1) S5:(-1,0)
//   Window2: 1-based [501,1001) -> 0-based [500,1000): calls
//     S1:(1,1) S2:(0,0) S3:(1,0) S4:(2,-1) S5:(2,2)
//
// Variants (0-based pos): v1=100 (window1), vMulti=200 (window1, ALT="G,T",
// must be skipped), v2=300 (window1), v3=600 (window2), v4=800 (window2).
static const int K = 3;
static const uint32_t N = 5;
static const int8_t W1[5][2] = {{0, 0}, {1, 2}, {2, 2}, {0, 1}, {-1, 0}};
static const int8_t W2[5][2] = {{1, 1}, {0, 0}, {1, 0}, {2, -1}, {2, 2}};

// GT per variant: {a0, a1, miss0, miss1} per sample.
struct GtCell { int a0, a1, m0, m1; };
static const GtCell V1[5] = {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}, {1, 1, 0, 0}, {0, 0, 1, 0}};
static const GtCell V2[5] = {{1, 1, 0, 0}, {0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}, {1, 1, 0, 0}};
static const GtCell V3[5] = {{0, 0, 0, 0}, {1, 1, 0, 0}, {0, 1, 0, 0}, {1, 0, 0, 1}, {0, 1, 0, 0}};
static const GtCell V4[5] = {{1, 0, 0, 0}, {0, 0, 0, 0}, {1, 1, 0, 0}, {0, 1, 0, 0}, {1, 1, 0, 0}};

static void expectedCell(
    const int8_t (*win)[2],
    const GtCell *gt,
    uint32_t s,
    int k,
    double &dos,
    double &hap
) {
    const int c0 = win[s][0];
    const int c1 = win[s][1];
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    dos = 0.0;
    hap = 0.0;
    if (c0 == k) {
        hap += 1.0;
        dos += gt[s].m0 ? NaN : static_cast<double>(gt[s].a0);
    }
    if (c1 == k) {
        hap += 1.0;
        dos += gt[s].m1 ? NaN : static_cast<double>(gt[s].a1);
    }
}

static void writeVcfFixture(const std::string &path) {
    std::string s;
    s += "##fileformat=VCFv4.2\n";
    s += "##contig=<ID=chr20,length=1000000>\n";
    s += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    s += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\n";
    s += "chr20\t101\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1\t.|0\n";
    s += "chr20\t201\trsX\tA\tG,T\t.\tPASS\t.\tGT\t0|1\t1|2\t0|0\t1|1\t2|0\n"; // multiallelic: must be skipped
    s += "chr20\t301\trs2\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|1\t1|0\t0|0\t1|1\n";
    s += "chr20\t601\trs3\tG\tA\t.\tPASS\t.\tGT\t0|0\t1|1\t0|1\t1|.\t0|1\n";
    s += "chr20\t801\trs4\tT\tC\t.\tPASS\t.\tGT\t1|0\t0|0\t1|1\t0|1\t1|1\n";
    writeFile(path, s);
}

static void writeMspFixture(const std::string &path) {
    std::string s;
    s += "#Subpopulation order/codes: 0=AFR\t1=EUR\t2=EAS\n";
    s += "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1\tS3.0\tS3.1\tS4.0\tS4.1\tS5.0\tS5.1\n";
    // 1-based RFMix spos/epos: [1,501) and [501,1001) -> 0-based [0,500) / [500,1000).
    s += "chr20\t1\t501\t0.0\t0.5\t50\t0\t0\t1\t2\t2\t2\t0\t1\t-1\t0\n";
    s += "chr20\t501\t1001\t0.5\t1.0\t50\t1\t1\t0\t0\t1\t0\t2\t-1\t2\t2\n";
    writeFile(path, s);
}

// ── Focused boundary case (regresses the 0-based [spos, epos) fix) ─────────
//
// K=3, N=2 (T1,T2), chr21.  Three contiguous windows tiling half-open in
// 1-based RFMix coordinates (epos_i == spos_{i+1}); each sets a distinct
// constant ancestry so the decoded ancestry identifies the window:
//   W_a: 1-based [1,11)  -> 0-based [0,10)   n snps = 3   ancestry 0
//   W_b: 1-based [11,21) -> 0-based [10,20)  n snps = 2   ancestry 1
//   W_c: 1-based [21,31) -> 0-based [20,30)  n snps = 4   ancestry 2
static void writeBoundaryMspFixture(const std::string &path) {
    std::string s;
    s += "#Subpopulation order/codes: 0=A0\t1=A1\t2=A2\n";
    s += "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tT1.0\tT1.1\tT2.0\tT2.1\n";
    s += "chr21\t1\t11\t0.00\t0.10\t3\t0\t0\t0\t0\n";
    s += "chr21\t11\t21\t0.10\t0.20\t2\t1\t1\t1\t1\n";
    s += "chr21\t21\t31\t0.20\t0.30\t4\t2\t2\t2\t2\n";
    writeFile(path, s);
}

// 9 biallelic SNVs.  1-based POS (pos0 = POS-1) placed so each window's
// assigned SNV count equals its synthetic `n snps`, with three SNVs sitting
// EXACTLY on window spos boundaries: POS=1 is spos0 of the first window,
// POS=11 is spos1, POS=21 is spos2.
//   POS  1, 3, 5    (pos0  0, 2, 4)              -> W_a (anc0), 3 sites
//   POS 11,15       (pos0 10,14)                 -> W_b (anc1), 2 sites
//   POS 21,23,25,27 (pos0 20,22,24,26)           -> W_c (anc2), 4 sites
static void writeBoundaryVcfFixture(const std::string &path) {
    std::string s;
    s += "##fileformat=VCFv4.2\n";
    s += "##contig=<ID=chr21,length=1000000>\n";
    s += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    s += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\tT2\n";
    const int pos[9] = {1, 3, 5, 11, 15, 21, 23, 25, 27};
    for (int i = 0; i < 9; ++i) {
        s += "chr21\t" + std::to_string(pos[i]) + "\tb" + std::to_string(i + 1) +
             "\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\n";
    }
    writeFile(path, s);
}

static void runBoundaryCase(const std::string &baseDir) {
    const std::string dir = baseDir + "/boundary";
    std::error_code ec;
    fs::create_directories(dir, ec);
    const std::string mspPrefix = dir + "/msp";
    const std::string genoPrefix = dir + "/geno";
    const std::string outPrefix = dir + "/out";
    writeBoundaryMspFixture(mspPrefix + ".chr21.msp.tsv");
    writeBoundaryVcfFixture(genoPrefix + ".chr21.vcf");

    try {
        convertRfmixToLanc(genoPrefix, /*expectBcf=*/false, mspPrefix, outPrefix, /*compressionLevel=*/3);
    } catch (const std::exception &e) {
        check(false, std::string("boundary convertRfmixToLanc threw: ") + e.what());
        return;
    }

    constexpr int nWin = 3;
    const int expNSnps[nWin] = {3, 2, 4}; // W_a, W_b, W_c synthetic `n snps`
    const int winAnc[nWin] = {0, 1, 2};   // each window's distinct constant ancestry
    constexpr int nMk = 9;
    const int expTotal = 9;
    // Query SNVs in .bim (marker) order and the window each must land in under
    // the corrected 0-based [spos, epos) rule.
    const int expWin[nMk] = {0, 0, 0, 1, 1, 2, 2, 2, 2};
    const int expPos1based[nMk] = {1, 3, 5, 11, 15, 21, 23, 25, 27};

    // (a) No SNV dropped: .bim has all 9 rows, and the FIRST row is POS==1
    //     (spos0 of the first window) — the first-window first SNP is kept.
    {
        std::ifstream bim(outPrefix + ".chr21.bim");
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(bim, line))
            if (!line.empty()) lines.push_back(line);
        check(lines.size() == static_cast<size_t>(expTotal),
              "boundary: bim has 9 rows (no SNV dropped)");
        if (!lines.empty())
            check(lines[0] == "chr21\tb1\t0\t1\tA\tG",
                  "boundary: first bim row is POS==1 (first-window first SNP NOT dropped)");
    }

    const uint32_t Nb = 2;
    const int Kb = 3;
    const uint32_t nWords = (Nb + 63) / 64;
    std::vector<uint64_t> maskAll(nWords, 0);
    for (uint32_t i = 0; i < Nb; ++i) maskAll[i / 64] |= (uint64_t(1) << (i % 64));

    LancData data(outPrefix, maskAll, Nb, Nb);
    check(data.nAncestries() == Kb, "boundary: nAncestries == 3");
    check(data.nMarkers() == static_cast<uint32_t>(expTotal), "boundary: nMarkers == 9 (no SNV dropped)");
    if (data.nMarkers() != static_cast<uint32_t>(expTotal)) return;

    auto cur = data.makeCursor();
    cur->beginSequentialBlock(0);

    int assignedByWin[nWin] = {0, 0, 0};
    for (int m = 0; m < nMk; ++m) {
        Eigen::MatrixXd dosMat(Nb, Kb), hapMat(Nb, Kb);
        cur->getAllAncestries(static_cast<uint64_t>(m), dosMat, hapMat);
        // Each window sets both haplotypes of both samples to one ancestry, so
        // exactly one column k has hapcount 2 for every sample: that k is the
        // window's ancestry, which identifies the assigned window.
        int decodedAnc = -1;
        for (int k = 0; k < Kb; ++k)
            if (hapMat(0, k) == 2.0 && hapMat(1, k) == 2.0) { decodedAnc = k; break; }
        check(decodedAnc >= 0,
              std::string("boundary: marker ") + std::to_string(m) + " decodes to a constant ancestry");
        if (decodedAnc < 0) continue;
        // (b) Each SNV — including the boundary SNVs (POS 1, 11, 21) — is
        //     assigned to the window whose 1-based [spos, epos) contains it.
        const int wantAnc = winAnc[expWin[m]];
        check(decodedAnc == wantAnc,
              std::string("boundary: marker ") + std::to_string(m) + " (POS=" +
                  std::to_string(expPos1based[m]) + ") -> window " + std::to_string(expWin[m]) +
                  " (ancestry " + std::to_string(wantAnc) + "), decoded ancestry " +
                  std::to_string(decodedAnc));
        for (int w = 0; w < nWin; ++w)
            if (winAnc[w] == decodedAnc) { ++assignedByWin[w]; break; }
    }

    // (c) Per-window assigned SNV count == synthetic `n snps` (per-window + total).
    int assignedTotal = 0;
    for (int w = 0; w < nWin; ++w) {
        assignedTotal += assignedByWin[w];
        check(assignedByWin[w] == expNSnps[w],
              std::string("boundary: window ") + std::to_string(w) + " assigned " +
                  std::to_string(assignedByWin[w]) + " SNVs == n snps " + std::to_string(expNSnps[w]));
    }
    check(assignedTotal == expTotal,
          std::string("boundary: total assigned ") + std::to_string(assignedTotal) +
              " == 9 (sum of n snps)");

    std::fprintf(stderr,
                 "[boundary] per-window assigned vs n snps: W_a %d/%d, W_b %d/%d, W_c %d/%d (total %d/%d)\n",
                 assignedByWin[0], expNSnps[0], assignedByWin[1], expNSnps[1], assignedByWin[2], expNSnps[2],
                 assignedTotal, expTotal);
}

int main(int argc, char **argv) {
    std::string baseDir = (argc > 1) ? argv[1] : "./lanc_conv_smoke_tmp";
    std::error_code ec;
    fs::remove_all(baseDir, ec);
    fs::create_directories(baseDir, ec);
    std::fprintf(stderr, "[lanc_convert_rfmix_smoke] work dir: %s\n", baseDir.c_str());

    const std::string genoPrefix = baseDir + "/geno";
    const std::string mspPrefix = baseDir + "/msp";
    const std::string outPrefix = baseDir + "/out";

    writeVcfFixture(genoPrefix + ".chr20.vcf");
    writeMspFixture(mspPrefix + ".chr20.msp.tsv");

    // ── Pass 1: full cohort, no --keep/--remove ──────────────────────────
    try {
        convertRfmixToLanc(genoPrefix, /*expectBcf=*/false, mspPrefix, outPrefix, /*compressionLevel=*/3);
    } catch (const std::exception &e) {
        check(false, std::string("convertRfmixToLanc threw: ") + e.what());
        std::fprintf(stderr, "\n[lanc_convert_rfmix_smoke] %d checks run, %d failures\n", g_checks, g_failures);
        return 1;
    }

    check(fs::exists(outPrefix + ".chr20.lanc"), "out.chr20.lanc exists (token stripped from 'chr20' -> '20')");
    check(fs::exists(outPrefix + ".chr20.bim"), "out.chr20.bim exists");
    check(fs::exists(outPrefix + ".fam"), "out.fam exists");

    // ── .bim: 4 kept markers (multiallelic vMulti skipped), correct columns ──
    {
        std::ifstream bim(outPrefix + ".chr20.bim");
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(bim, line))
            if (!line.empty()) lines.push_back(line);
        check(lines.size() == 4, "bim has 4 rows (multiallelic site skipped)");
        if (lines.size() == 4) {
            check(lines[0] == "chr20\trs1\t0\t101\tA\tG", "bim row 1 (rs1) matches");
            check(lines[1] == "chr20\trs2\t0\t301\tC\tT", "bim row 2 (rs2) matches");
            check(lines[2] == "chr20\trs3\t0\t601\tG\tA", "bim row 3 (rs3) matches");
            check(lines[3] == "chr20\trs4\t0\t801\tT\tC", "bim row 4 (rs4) matches");
        }
    }

    // ── .fam: 5 samples, original order ──────────────────────────────────
    {
        std::ifstream fam(outPrefix + ".fam");
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(fam, line))
            if (!line.empty()) lines.push_back(line);
        check(lines.size() == 5, "fam has 5 rows");
        const char *expIID[5] = {"S1", "S2", "S3", "S4", "S5"};
        for (int i = 0; i < 5 && i < static_cast<int>(lines.size()); ++i) {
            std::string exp = std::string(expIID[i]) + "\t" + expIID[i] + "\t0\t0\t0\t-9";
            check(lines[static_cast<size_t>(i)] == exp, std::string("fam row ") + std::to_string(i) + " == " + exp);
        }
    }

    // ── Decode via LancData/LancCursor and check every (marker,sample,k) ──
    {
        const uint32_t nWords = (N + 63) / 64;
        std::vector<uint64_t> maskAll(nWords, 0);
        for (uint32_t i = 0; i < N; ++i) maskAll[i / 64] |= (uint64_t(1) << (i % 64));

        LancData data(outPrefix, maskAll, N, N);
        check(data.nAncestries() == K, "nAncestries == 3");
        check(data.nMarkers() == 4, "nMarkers == 4");
        check(data.nSubjInFile() == N, "nSubjInFile == 5");
        check(data.hasNoMissing() == false, "hasNoMissing == false (fixture has 2 missing genotypes)");

        auto cur = data.makeCursor();
        cur->beginSequentialBlock(0);

        const int8_t(*wins[4])[2] = {W1, W1, W2, W2};
        const GtCell *gts[4] = {V1, V2, V3, V4};
        const char *names[4] = {"rs1", "rs2", "rs3", "rs4"};

        for (int m = 0; m < 4; ++m) {
            Eigen::MatrixXd dosMat(N, K), hapMat(N, K);
            cur->getAllAncestries(static_cast<uint64_t>(m), dosMat, hapMat);
            for (uint32_t s = 0; s < N; ++s) {
                for (int k = 0; k < K; ++k) {
                    double ed, eh;
                    expectedCell(wins[m], gts[m], s, k, ed, eh);
                    bool ok = feq(dosMat(static_cast<int>(s), k), ed) && feq(hapMat(static_cast<int>(s), k), eh);
                    check(ok, std::string("marker ") + names[m] + " sample S" + std::to_string(s + 1) + " k=" +
                                  std::to_string(k) + " (got dos=" + std::to_string(dosMat(static_cast<int>(s), k)) +
                                  " hap=" + std::to_string(hapMat(static_cast<int>(s), k)) + ", want dos=" +
                                  std::to_string(ed) + " hap=" + std::to_string(eh) + ")");
                }
            }
        }
    }

    // ── Pass 2: --keep S5,S1,S3 (out of order in the file) -> fam/decode
    //    must preserve ORIGINAL BCF/msp order (S1, S3, S5), not keep-file order.
    {
        const std::string keepFile = baseDir + "/keep.txt";
        writeFile(keepFile, "S5\nS1\nS3\n");
        const std::string outPrefix2 = baseDir + "/outkeep";

        convertRfmixToLanc(genoPrefix, /*expectBcf=*/false, mspPrefix, outPrefix2, /*compressionLevel=*/3, keepFile);

        std::ifstream fam(outPrefix2 + ".fam");
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(fam, line))
            if (!line.empty()) lines.push_back(line);
        check(lines.size() == 3, "--keep fam has 3 rows");
        const char *expIID[3] = {"S1", "S3", "S5"}; // original order, not keep-file order
        for (int i = 0; i < 3 && i < static_cast<int>(lines.size()); ++i) {
            std::string exp = std::string(expIID[i]) + "\t" + expIID[i] + "\t0\t0\t0\t-9";
            check(lines[static_cast<size_t>(i)] == exp,
                  std::string("--keep fam row ") + std::to_string(i) + " == " + exp + " (order preserved)");
        }

        const uint32_t nUsed = 3;
        const uint32_t nWords = (nUsed + 63) / 64;
        std::vector<uint64_t> maskAll(nWords, 0);
        for (uint32_t i = 0; i < nUsed; ++i) maskAll[i / 64] |= (uint64_t(1) << (i % 64));
        LancData data2(outPrefix2, maskAll, nUsed, nUsed);

        auto cur2 = data2.makeCursor();
        Eigen::MatrixXd dosMat(nUsed, K), hapMat(nUsed, K);
        cur2->getAllAncestries(0, dosMat, hapMat); // marker rs1
        // kept subject order is S1(idx0), S3(idx2), S5(idx4) in the original W1/V1 arrays
        const uint32_t origIdx[3] = {0, 2, 4};
        for (uint32_t di = 0; di < nUsed; ++di) {
            for (int k = 0; k < K; ++k) {
                double ed, eh;
                expectedCell(W1, V1, origIdx[di], k, ed, eh);
                bool ok = feq(dosMat(static_cast<int>(di), k), ed) && feq(hapMat(static_cast<int>(di), k), eh);
                check(ok, std::string("--keep marker rs1 row ") + std::to_string(di) + " k=" + std::to_string(k));
            }
        }
    }

    // ── Focused boundary case: 3 contiguous windows, boundary-placed SNVs ──
    runBoundaryCase(baseDir);

    std::fprintf(stderr, "\n[lanc_convert_rfmix_smoke] %d checks run, %d failures\n", g_checks, g_failures);
    if (g_failures == 0) {
        std::printf("PASS: all %d checks passed\n", g_checks);
        return 0;
    }
    std::printf("FAIL: %d of %d checks failed\n", g_failures, g_checks);
    return 1;
}
