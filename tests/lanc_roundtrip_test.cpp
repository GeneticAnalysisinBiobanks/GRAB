// lanc_roundtrip_test.cpp — self-contained P1 gate for the .lanc v1 format.
//
// Exercises LancWriter / LancData / LancCursor over two synthetic configurations
// (NO_MISSING and WITH-missing), multiple chromosomes, N=300, M spanning several
// blocks (blockLen=512), K=4 with unassigned (0xF) nibbles and missing bits.
//
// Checks:
//   1. Deterministic PRNG-generated per-(marker,individual,hap) source arrays.
//   2. getAllAncestries decode == expected dosage_k/hapcount_k (formula), both
//      sequentially and via random jumps.
//   3. getAdmixGenotypes(marker,k) == column k of getAllAncestries; altFreq ==
//      sum(dos)/sum(hap) (guarded).
//   4. Subject filtering (every other subject) yields the correct N_used rows.
//   5. Multi-chr prefix stream == chr1 markers followed by chr2 markers.
//   6. Determinism: encoding twice yields byte-identical .lanc files.
//
// Build (example):
//   g++ -std=c++17 -O2 -Isrc -Ithird_party/eigen-5.0.0 -Ithird_party/zstd-1.5.7/lib
//       tests/lanc_roundtrip_test.cpp src/localplus/lanc_io.cpp
//       src/geno_factory/variant_filter.cpp build/zstd/(common,compress,decompress)/*.o
//       -lstdc++fs -o <bin>
//
// Prints a PASS/FAIL summary; exits non-zero on any failure.

#include "localplus/lanc_io.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ── failure tracking ────────────────────────────────────────────────────
static int g_failures = 0;
static int g_checks = 0;

static void check(bool cond, const std::string &msg) {
    ++g_checks;
    if (!cond) {
        ++g_failures;
        std::fprintf(stderr, "  FAIL: %s\n", msg.c_str());
    }
}

// NaN==NaN treated as equal; otherwise exact value equality.
static bool feq(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    return a == b;
}

// ── deterministic PRNG (splitmix64, fixed seed) ─────────────────────────
struct Splitmix64 {
    uint64_t state;
    explicit Splitmix64(uint64_t seed) : state(seed) {}
    uint64_t next() {
        uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};

// ── synthetic dataset ───────────────────────────────────────────────────
struct ChrArrays {
    std::string chrom;
    uint32_t M;
    std::vector<uint8_t> anc0, anc1, alt0, alt1, miss0, miss1; // each size M*N
};

struct Dataset {
    bool withMissing;
    uint32_t N, K;
    uint16_t blockLen;
    std::vector<ChrArrays> chrs;
};

static Dataset generateDataset(bool withMissing, uint64_t seed) {
    Dataset ds;
    ds.withMissing = withMissing;
    ds.N = 300;
    ds.K = 4;
    ds.blockLen = 512;

    Splitmix64 rng(seed);
    struct ChrSpec { const char *chrom; uint32_t M; };
    const ChrSpec specs[2] = {{"1", 1300}, {"2", 1000}};

    auto drawAnc = [&](uint32_t K) -> uint8_t {
        uint64_t r = rng.next();
        if ((r & 0xFF) < 26) return LANC_ANC_UNASSIGNED; // ~10% unassigned
        return static_cast<uint8_t>((r >> 8) % K);
    };
    auto drawBit = [&]() -> uint8_t { return static_cast<uint8_t>(rng.next() & 1); };
    auto drawMiss = [&]() -> uint8_t { return static_cast<uint8_t>((rng.next() & 0xFF) < 13 ? 1 : 0); }; // ~5%

    for (int c = 0; c < 2; ++c) {
        ChrArrays ca;
        ca.chrom = specs[c].chrom;
        ca.M = specs[c].M;
        const size_t n = static_cast<size_t>(ca.M) * ds.N;
        ca.anc0.resize(n);
        ca.anc1.resize(n);
        ca.alt0.resize(n);
        ca.alt1.resize(n);
        ca.miss0.assign(n, 0);
        ca.miss1.assign(n, 0);
        for (size_t idx = 0; idx < n; ++idx) {
            ca.anc0[idx] = drawAnc(ds.K);
            ca.anc1[idx] = drawAnc(ds.K);
            ca.alt0[idx] = drawBit();
            ca.alt1[idx] = drawBit();
            if (withMissing) {
                ca.miss0[idx] = drawMiss();
                ca.miss1[idx] = drawMiss();
            }
        }
        ds.chrs.push_back(std::move(ca));
    }
    return ds;
}

// ── file writers ────────────────────────────────────────────────────────
static void writeFam(const std::string &path, uint32_t N) {
    std::ofstream f(path);
    for (uint32_t i = 0; i < N; ++i)
        f << "F" << i << " IND" << i << " 0 0 0 -9\n";
}

static void writeBim(const std::string &path, const ChrArrays &ca) {
    std::ofstream f(path);
    for (uint32_t m = 0; m < ca.M; ++m)
        f << ca.chrom << " snp_" << ca.chrom << "_" << m << " 0 " << (m + 1) << " A G\n";
}

static void writeDataset(const std::string &prefix, const Dataset &ds) {
    for (const auto &ca : ds.chrs) {
        std::string lancPath = prefix + ".chr" + ca.chrom + ".lanc";
        std::string bimPath = prefix + ".chr" + ca.chrom + ".bim";
        LancWriter w(lancPath, static_cast<uint8_t>(ds.K), ds.N, ds.blockLen, 3);
        for (uint32_t m = 0; m < ca.M; ++m) {
            const size_t off = static_cast<size_t>(m) * ds.N;
            const uint8_t *miss0 = ds.withMissing ? &ca.miss0[off] : nullptr;
            const uint8_t *miss1 = ds.withMissing ? &ca.miss1[off] : nullptr;
            w.addMarker(&ca.anc0[off], &ca.anc1[off], &ca.alt0[off], &ca.alt1[off], miss0, miss1);
        }
        w.close();
        writeBim(bimPath, ca);
    }
    writeFam(prefix + ".fam", ds.N);
}

// ── expected decode formula ─────────────────────────────────────────────
static void expectedCell(
    const ChrArrays &ca,
    uint32_t N,
    uint32_t m,
    uint32_t i,
    int k,
    double &dos,
    double &hap
) {
    const size_t idx = static_cast<size_t>(m) * N + i;
    const int c0 = ca.anc0[idx];
    const int c1 = ca.anc1[idx];
    const uint8_t a0 = ca.alt0[idx] & 1;
    const uint8_t a1 = ca.alt1[idx] & 1;
    const uint8_t m0 = ca.miss0[idx] & 1;
    const uint8_t m1 = ca.miss1[idx] & 1;
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    hap = static_cast<double>(c0 == k) + static_cast<double>(c1 == k);
    dos = (c0 == k ? (m0 ? NaN : static_cast<double>(a0)) : 0.0) +
          (c1 == k ? (m1 ? NaN : static_cast<double>(a1)) : 0.0);
}

// Locate the chromosome + local marker for a global (concatenated) index.
static const ChrArrays &locate(const Dataset &ds, uint32_t globalMarker, uint32_t &localM) {
    uint32_t acc = 0;
    for (const auto &ca : ds.chrs) {
        if (globalMarker < acc + ca.M) {
            localM = globalMarker - acc;
            return ca;
        }
        acc += ca.M;
    }
    localM = 0;
    return ds.chrs.back();
}

// ── byte-identity of two files ──────────────────────────────────────────
static bool filesIdentical(const std::string &a, const std::string &b) {
    std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
    if (!fa || !fb) return false;
    std::vector<char> ba((std::istreambuf_iterator<char>(fa)), std::istreambuf_iterator<char>());
    std::vector<char> bb((std::istreambuf_iterator<char>(fb)), std::istreambuf_iterator<char>());
    return ba == bb;
}

// ── verification for one dataset/prefix ─────────────────────────────────
static void verifyDataset(const std::string &prefix, const Dataset &ds, const std::string &label) {
    const uint32_t N = ds.N;
    const int K = static_cast<int>(ds.K);
    const uint32_t totalM = ds.chrs[0].M + ds.chrs[1].M;

    // usedMask: all subjects included.
    const uint32_t nWords = (N + 63) / 64;
    std::vector<uint64_t> maskAll(nWords, 0);
    for (uint32_t i = 0; i < N; ++i) maskAll[i / 64] |= (uint64_t(1) << (i % 64));

    LancData data(prefix, maskAll, N, N);

    check(data.nAncestries() == K, label + ": nAncestries == K");
    check(data.nMarkers() == totalM, label + ": nMarkers == totalM");
    check(data.nSubjInFile() == N, label + ": nSubjInFile == N");
    check(data.nSubjUsed() == N, label + ": nSubjUsed == N");
    check(data.hasNoMissing() == !ds.withMissing, label + ": hasNoMissing flag matches config");

    // (5) multi-chr prefix stream = chr1 markers then chr2 markers.
    {
        bool streamOk = true;
        for (uint32_t g = 0; g < totalM; ++g) {
            uint32_t localM;
            const ChrArrays &ca = locate(ds, g, localM);
            std::string expChrom = ca.chrom;
            std::string expId = "snp_" + ca.chrom + "_" + std::to_string(localM);
            if (std::string(data.chr(g)) != expChrom) streamOk = false;
            if (std::string(data.markerId(g)) != expId) streamOk = false;
            if (data.pos(g) != localM + 1) streamOk = false;
        }
        check(streamOk, label + ": concatenated marker stream == chr1 then chr2");
    }

    // (2)+(3) full decode via getAllAncestries + getAdmixGenotypes.
    auto verifyMarker = [&](LancCursor &cur, uint32_t g) {
        uint32_t localM;
        const ChrArrays &ca = locate(ds, g, localM);

        Eigen::MatrixXd dosMat(N, K), hapMat(N, K);
        cur.getAllAncestries(g, dosMat, hapMat);

        for (uint32_t i = 0; i < N; ++i) {
            for (int k = 0; k < K; ++k) {
                double ed, eh;
                expectedCell(ca, N, localM, i, k, ed, eh);
                if (!feq(dosMat(i, k), ed) || !feq(hapMat(i, k), eh)) {
                    check(false, label + ": getAllAncestries mismatch at g=" + std::to_string(g) +
                                     " i=" + std::to_string(i) + " k=" + std::to_string(k));
                    return; // avoid flooding
                }
            }
        }

        // getAdmixGenotypes per ancestry.
        for (int k = 0; k < K; ++k) {
            Eigen::VectorXd dosCol(N), hapCol(N);
            double af = cur.getAdmixGenotypes(g, k, dosCol, hapCol);
            bool colOk = true;
            for (uint32_t i = 0; i < N; ++i)
                if (!feq(dosCol(i), dosMat(i, k)) || !feq(hapCol(i), hapMat(i, k))) colOk = false;
            check(colOk, label + ": getAdmixGenotypes col == getAllAncestries col (g=" +
                             std::to_string(g) + ", k=" + std::to_string(k) + ")");
            double dosSum = dosCol.sum();
            double hapSum = hapCol.sum();
            double expAf = (hapSum > 0.0) ? dosSum / hapSum : 0.0;
            check(feq(af, expAf), label + ": altFreq == sum(dos)/sum(hap) (g=" +
                                      std::to_string(g) + ", k=" + std::to_string(k) + ")");
        }
    };

    // Sequential pass over every marker.
    {
        auto cur = data.makeCursor();
        cur->beginSequentialBlock(0);
        for (uint32_t g = 0; g < totalM; ++g) verifyMarker(*cur, g);
    }

    // Random-jump pass (deterministic permutation).
    {
        auto cur = data.makeCursor();
        std::vector<uint32_t> order(totalM);
        std::iota(order.begin(), order.end(), 0u);
        // Deterministic shuffle via splitmix64.
        Splitmix64 rng(0xABCDEF0123456789ULL);
        for (uint32_t i = totalM; i > 1; --i) {
            uint32_t j = static_cast<uint32_t>(rng.next() % i);
            std::swap(order[i - 1], order[j]);
        }
        for (uint32_t g : order) verifyMarker(*cur, g);
    }

    // (4) subject filtering: every other subject.
    {
        std::vector<uint32_t> selected;
        std::vector<uint64_t> maskEven(nWords, 0);
        for (uint32_t i = 0; i < N; i += 2) {
            maskEven[i / 64] |= (uint64_t(1) << (i % 64));
            selected.push_back(i);
        }
        const uint32_t nUsed = static_cast<uint32_t>(selected.size());

        LancData dataSub(prefix, maskEven, N, nUsed);
        check(dataSub.nSubjUsed() == nUsed, label + ": subset nSubjUsed");
        auto cur = dataSub.makeCursor();

        // Check a spread of markers across both chromosomes.
        std::vector<uint32_t> probes = {0, 1, 511, 512, 513, ds.chrs[0].M - 1,
                                        ds.chrs[0].M, ds.chrs[0].M + 512, totalM - 1};
        for (uint32_t g : probes) {
            if (g >= totalM) continue;
            uint32_t localM;
            const ChrArrays &ca = locate(ds, g, localM);
            Eigen::MatrixXd dosMat(nUsed, K), hapMat(nUsed, K);
            cur->getAllAncestries(g, dosMat, hapMat);
            bool ok = true;
            for (uint32_t di = 0; di < nUsed; ++di) {
                uint32_t s = selected[di];
                for (int k = 0; k < K; ++k) {
                    double ed, eh;
                    expectedCell(ca, N, localM, s, k, ed, eh);
                    if (!feq(dosMat(di, k), ed) || !feq(hapMat(di, k), eh)) ok = false;
                }
            }
            check(ok, label + ": subset decode matches selected rows (g=" + std::to_string(g) + ")");
        }
    }
}

int main(int argc, char **argv) {
    std::string baseDir = (argc > 1) ? argv[1] : "./lanc_rt_tmp";
    std::error_code ec;
    fs::create_directories(baseDir, ec);

    std::fprintf(stderr, "[lanc_roundtrip] work dir: %s\n", baseDir.c_str());

    struct Cfg { const char *name; bool withMissing; uint64_t seed; };
    const Cfg cfgs[2] = {
        {"NO_MISSING", false, 0x1111111111111111ULL},
        {"WITH_MISSING", true, 0x2222222222222222ULL},
    };

    for (const auto &cfg : cfgs) {
        std::fprintf(stderr, "[lanc_roundtrip] === config %s ===\n", cfg.name);
        Dataset ds = generateDataset(cfg.withMissing, cfg.seed);

        std::string prefixA = baseDir + "/" + cfg.name + "_A";
        std::string prefixB = baseDir + "/" + cfg.name + "_B";
        writeDataset(prefixA, ds);
        writeDataset(prefixB, ds);

        // (6) determinism: A vs B byte-identical per chromosome file.
        for (const auto &ca : ds.chrs) {
            std::string fa = prefixA + ".chr" + ca.chrom + ".lanc";
            std::string fb = prefixB + ".chr" + ca.chrom + ".lanc";
            check(filesIdentical(fa, fb),
                  std::string(cfg.name) + ": determinism chr" + ca.chrom + " byte-identical");
        }

        verifyDataset(prefixA, ds, cfg.name);
    }

    std::fprintf(stderr, "\n[lanc_roundtrip] %d checks run, %d failures\n", g_checks, g_failures);
    if (g_failures == 0) {
        std::printf("PASS: all %d checks passed\n", g_checks);
        return 0;
    }
    std::printf("FAIL: %d of %d checks failed\n", g_failures, g_checks);
    return 1;
}
