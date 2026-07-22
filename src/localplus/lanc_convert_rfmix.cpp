// lanc_convert_rfmix.cpp — RFMix .msp.tsv + phased query BCF/VCF -> .lanc
//
// Multi-chromosome refactor of the single-file .abed converter that this
// module supersedes (see dev-notes/methods/recode_rfmix/
// 01_lanc_format_and_recoding.md section 6).  The MSP parser, chromosome
// comparator, and window cursor are copied unmodified from that converter;
// the per-ancestry `for k` collapse loop is deleted outright and replaced
// with direct per-haplotype plane emission via LancWriter::addMarker.
//
// Control flow:
//   1. Glob {mspPrefix}*.msp.tsv; parse each (K, sample IDs, windows);
//      verify each file is single-chromosome and derive its token (the
//      chromosome string with a leading "chr"/"CHR" stripped).
//   2. Verify all msp files agree on K and on the sample-ID list + order
//      (the query cohort is shared across chromosomes).
//   3. Glob {genoPrefix}*.bcf (or *.vcf.gz/*.vcf.zst/*.vcf); for each,
//      cross-check format vs expectBcf and determine its chromosome from
//      the first record (or, for an empty file, a single header contig);
//      match to an msp entry by token.
//   4. Error if any chromosome has an msp file without a matching geno
//      file, or vice versa. Process matched chromosomes in sorted order.
//   5. Per chromosome: stream biallelic SNVs; for each variant whose 0-based
//      POS falls in an msp window (monotone WindowCursor), build per-kept-
//      sample per-haplotype ancestry/allele/missing arrays and call
//      LancWriter::addMarker. Emit a matching .bim line.
//   6. Write one shared .fam after all chromosomes (filtered query samples,
//      original BCF/msp order).
#include "localplus/lanc_convert_rfmix.hpp"
#include "localplus/lanc_io.hpp"
#include "io/subject_filter.hpp"
#include "util/logging.hpp"

#include <zlib.h>
extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// ── simple line reader for gzFile (copied from the old .abed converter) ────

bool mspReadLine(
    gzFile gz,
    std::string &line
) {
    line.clear();
    char buf[131072];
    while (gzgets(gz, buf, sizeof(buf)) != nullptr) {
        size_t len = std::strlen(buf);
        line.append(buf, len);
        if (len > 0 && buf[len - 1] == '\n') {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
                line.pop_back();
            return true;
        }
    }
    return !line.empty();
}

std::vector<std::string> mspSplitTabs(const std::string &s) {
    std::vector<std::string> t;
    size_t start = 0;
    while (start <= s.size()) {
        auto pos = s.find('\t', start);
        if (pos == std::string::npos) {
            t.push_back(s.substr(start));
            break;
        }
        t.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return t;
}

// ── MSP data structures (copied from the old .abed converter) ──────────────

struct MspWindow {
    std::string chrom;
    int32_t spos;              // 0-based, inclusive
    int32_t epos;              // 0-based, exclusive
    std::vector<int8_t> calls; // length 2*nSamples; each is ancestry 0..K-1, or -1
};

struct MspData {
    int K = 0;
    uint32_t nSamples = 0;
    std::vector<std::string> sampleIDs;
    std::vector<MspWindow> windows;
};

// ── Parse MSP file (copied from the old .abed converter) ───────────────────

MspData parseMsp(const std::string &mspFile) {
    gzFile gz = gzopen(mspFile.c_str(), "rb");
    if (!gz) throw std::runtime_error("Cannot open MSP file: " + mspFile);

    MspData msp;
    std::string line;

    // Line 1: #Subpopulation order/codes: 0=X\t1=Y\t...
    if (!mspReadLine(gz, line)) {
        gzclose(gz);
        throw std::runtime_error("MSP file is empty: " + mspFile);
    }
    {
        auto sep = line.find(": ");
        if (sep == std::string::npos) {
            gzclose(gz);
            throw std::runtime_error("MSP line 1 missing ': ': " + line);
        }
        auto rest = line.substr(sep + 2);
        msp.K = 1;
        for (char c : rest)
            if (c == '\t') ++msp.K;
    }
    if (msp.K < 1) {
        gzclose(gz);
        throw std::runtime_error("MSP: could not detect K from line 1: " + mspFile);
    }
    if (msp.K > 15) {
        gzclose(gz);
        throw std::runtime_error("MSP: K=" + std::to_string(msp.K) +
                                 " exceeds the .lanc limit of 15 ancestries: " + mspFile);
    }

    // Line 2: #chm\tspos\tepos\tsgpos\tegpos\tn snps\tIID0.0\tIID0.1\tIID1.0\tIID1.1\t...
    if (!mspReadLine(gz, line)) {
        gzclose(gz);
        throw std::runtime_error("MSP file has only 1 line: " + mspFile);
    }
    {
        auto cols = mspSplitTabs(line);
        if (cols.size() < 7) {
            gzclose(gz);
            throw std::runtime_error("MSP header line 2 has fewer than 7 columns: " + mspFile);
        }
        size_t nHapCols = cols.size() - 6;
        if (nHapCols == 0 || nHapCols % 2 != 0) {
            gzclose(gz);
            throw std::runtime_error("MSP header line 2: expected even number of subject columns (got " +
                                     std::to_string(nHapCols) + "): " + mspFile);
        }
        msp.nSamples = static_cast<uint32_t>(nHapCols / 2);
        msp.sampleIDs.resize(msp.nSamples);
        for (uint32_t s = 0; s < msp.nSamples; ++s) {
            const auto &hapID = cols[6 + s * 2]; // e.g. "SAMPLE1.0"
            auto dot = hapID.rfind('.');
            msp.sampleIDs[s] = (dot != std::string::npos) ? hapID.substr(0, dot) : hapID;
        }
    }

    // Data rows
    while (mspReadLine(gz, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cols = mspSplitTabs(line);
        size_t expected = 6 + 2 * msp.nSamples;
        if (cols.size() < expected) {
            gzclose(gz);
            throw std::runtime_error("MSP data row has " + std::to_string(cols.size()) + " columns, expected " +
                                     std::to_string(expected) + ": " + mspFile);
        }
        MspWindow w;
        w.chrom = cols[0];
        w.spos = std::stoi(cols[1]);
        w.epos = std::stoi(cols[2]);
        w.calls.resize(2 * msp.nSamples);
        for (uint32_t h = 0; h < 2 * msp.nSamples; ++h) {
            int call = std::stoi(cols[6 + h]);
            w.calls[h] = static_cast<int8_t>(call >= 0 && call < msp.K ? call : -1);
        }
        msp.windows.push_back(std::move(w));
    }

    gzclose(gz);

    if (msp.windows.empty()) throw std::runtime_error("MSP file has no data rows: " + mspFile);

    infoMsg("MSP %s: K=%d, %u subjects, %zu windows", mspFile.c_str(), msp.K, msp.nSamples, msp.windows.size());
    return msp;
}

// ── Chromosome token / comparison helpers ────────────────────────────────
// stripChrPrefix / chromLessThan are copied from the old .abed converter (and
// match the identical helper in lanc_io.cpp, which is file-local there) so
// the output filenames are discoverable by LancData's {prefix}.chr*.lanc
// glob and sort in the same numeric-aware order the reader expects.

std::string stripChrPrefix(const std::string &s) {
    if (s.size() >= 3 && (s[0] == 'c' || s[0] == 'C') && (s[1] == 'h' || s[1] == 'H') &&
        (s[2] == 'r' || s[2] == 'R'))
        return s.substr(3);
    return s;
}

bool chromLessThan(
    const std::string &a,
    const std::string &b
) {
    if (a == b) return false;
    std::string sa = stripChrPrefix(a), sb = stripChrPrefix(b);
    bool aNum = !sa.empty() && sa.find_first_not_of("0123456789") == std::string::npos;
    bool bNum = !sb.empty() && sb.find_first_not_of("0123456789") == std::string::npos;
    if (aNum && bNum) return std::stoll(sa) < std::stoll(sb);
    if (aNum != bNum) return aNum;
    return sa < sb;
}

// ── Window lookup cursor (copied from the old .abed converter) ─────────────

struct WindowCursor {
    const std::vector<MspWindow> &windows;
    size_t idx = 0;

    int findWindow(
        const std::string &chrom,
        int32_t pos0
    ) {
        while (idx < windows.size()) {
            const auto &w = windows[idx];
            if (chromLessThan(w.chrom, chrom) || (w.chrom == chrom && w.epos <= pos0)) {
                ++idx;
                continue;
            }
            break;
        }
        if (idx >= windows.size()) return -1;
        const auto &w = windows[idx];
        if (w.chrom == chrom && w.spos <= pos0 && pos0 < w.epos) return static_cast<int>(idx);
        return -1;
    }
};

// ── Directory glob: {prefix}<anything>{suffix}, any suffix in `suffixes` ──

std::vector<std::string> globFiles(
    const std::string &prefix,
    const std::vector<std::string> &suffixes
) {
    namespace fs = std::filesystem;
    fs::path prefixPath(prefix);
    fs::path parentDir = prefixPath.has_parent_path() ? prefixPath.parent_path() : fs::path(".");
    const std::string base = prefixPath.filename().string();

    std::vector<std::string> out;
    std::error_code ec;
    if (!fs::exists(parentDir, ec)) return out;
    for (const auto &entry : fs::directory_iterator(parentDir, ec)) {
        if (!entry.is_regular_file()) continue;
        const std::string fn = entry.path().filename().string();
        if (fn.size() < base.size() || fn.compare(0, base.size(), base) != 0) continue;
        for (const auto &suf : suffixes) {
            if (fn.size() >= suf.size() && fn.compare(fn.size() - suf.size(), suf.size(), suf) == 0) {
                out.push_back(entry.path().string());
                break;
            }
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

// ── Format cross-check, matching the old .abed converter's processVcf ──────

void checkGenoFormat(
    htsFile *fp,
    bool expectBcf,
    const std::string &path
) {
    const htsFormat *fmt = hts_get_format(fp);
    const enum htsExactFormat detected = fmt ? fmt->format : unknown_format;
    if (!expectBcf && detected == bcf)
        throw std::runtime_error(path + " appears to be a BCF2 file. Try --bcf instead of --vcf.");
    if (expectBcf && detected != bcf)
        throw std::runtime_error(path + " does not appear to be a BCF2 file. Try --vcf instead of --bcf.");
}

// Determine a genotype file's chromosome from its first data record (or, for
// an empty file, its single header contig if unambiguous).
std::string detectGenoChrom(
    const std::string &path,
    bool expectBcf
) {
    const char *kindLabel = expectBcf ? "BCF" : "VCF";
    htsFile *fp = hts_open(path.c_str(), "r");
    if (!fp) throw std::runtime_error(std::string("Cannot open ") + kindLabel + ": " + path);
    checkGenoFormat(fp, expectBcf, path);

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        throw std::runtime_error(std::string("Cannot read ") + kindLabel + " header: " + path);
    }
    bcf1_t *rec = bcf_init();

    std::string chrom;
    if (bcf_read(fp, hdr, rec) == 0) {
        chrom = bcf_hdr_id2name(hdr, rec->rid);
    } else if (hdr->n[BCF_DT_CTG] == 1) {
        chrom = bcf_hdr_id2name(hdr, 0);
    } else {
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        throw std::runtime_error("Cannot determine chromosome for empty genotype file (" +
                                 std::to_string(hdr->n[BCF_DT_CTG]) + " header contigs, no records): " + path);
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    return chrom;
}

// ── Per-chromosome stream: MSP windows + genotypes -> LancWriter planes ──
// Direct refactor of the old .abed converter's processVcf (nthreads<=1 branch).
// The per-ancestry `for k` dosage/hapcount collapse loop is deleted; each
// marker instead emits one per-haplotype ancestry/allele/missing triple per
// kept sample straight into LancWriter::addMarker.
uint32_t processChrToLanc(
    const std::string &genoPath,
    bool expectBcf,
    const MspData &msp,
    LancWriter &writer,
    std::ofstream &bimOut,
    const std::vector<uint32_t> &keptIndices,
    uint64_t &nMissingOut,
    int nthreads
) {
    const char *kindLabel = expectBcf ? "BCF" : "VCF";

    htsFile *fp = hts_open(genoPath.c_str(), "r");
    if (!fp) throw std::runtime_error(std::string("Cannot open ") + kindLabel + ": " + genoPath);
    checkGenoFormat(fp, expectBcf, genoPath);

    if (nthreads > 1) {
        int nDecompressThreads = std::min(nthreads, 4);
        hts_set_threads(fp, nDecompressThreads);
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        throw std::runtime_error(std::string("Cannot read ") + kindLabel + " header: " + genoPath);
    }

    bcf1_t *rec = bcf_init();

    const int nSamplesVcf = bcf_hdr_nsamples(hdr);
    if (static_cast<uint32_t>(nSamplesVcf) != msp.nSamples) {
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        throw std::runtime_error("Subject count mismatch: " + genoPath + " has " +
                                 std::to_string(nSamplesVcf) + " subjects, MSP has " +
                                 std::to_string(msp.nSamples));
    }
    for (int s = 0; s < nSamplesVcf; ++s) {
        const char *vcfID = hdr->samples[s];
        if (msp.sampleIDs[static_cast<uint32_t>(s)] != vcfID) {
            bcf_destroy(rec);
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("Subject name mismatch at index " + std::to_string(s) + " in " + genoPath +
                                     ": genotype file has '" + vcfID + "', MSP has '" +
                                     msp.sampleIDs[static_cast<uint32_t>(s)] + "'");
        }
    }

    int32_t *gtArr = nullptr;
    int ngtArr = 0;

    const uint32_t Nall = msp.nSamples;
    const uint32_t Nout = keptIndices.empty() ? Nall : static_cast<uint32_t>(keptIndices.size());

    WindowCursor cursor{msp.windows};
    uint32_t nMatch = 0;
    uint64_t nMissing = 0;
    bool phasingChecked = false;

    std::vector<uint8_t> anc0(Nout), anc1(Nout), alt0(Nout), alt1(Nout), miss0(Nout), miss1(Nout);

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        if (rec->n_allele != 2) continue; // skip non-biallelic sites

        const char *chromStr = bcf_hdr_id2name(hdr, rec->rid);
        const int32_t pos0 = rec->pos;

        const int winIdx = cursor.findWindow(chromStr, pos0);
        if (winIdx < 0) continue;

        bcf_get_genotypes(hdr, rec, &gtArr, &ngtArr);
        if (!phasingChecked && ngtArr >= 2) {
            if (!bcf_gt_is_phased(gtArr[1]))
                throw std::runtime_error("--make-lanc: " + genoPath +
                                         " genotypes are unphased. A phased VCF/BCF (genotypes separated "
                                         "by '|') is required.");
            phasingChecked = true;
        }

        const MspWindow &w = msp.windows[static_cast<size_t>(winIdx)];
        for (uint32_t j = 0; j < Nout; ++j) {
            const uint32_t s = keptIndices.empty() ? j : keptIndices[j];

            const int8_t call0 = w.calls[2 * s + 0];
            const int8_t call1 = w.calls[2 * s + 1];
            anc0[j] = (call0 < 0) ? LANC_ANC_UNASSIGNED : static_cast<uint8_t>(call0);
            anc1[j] = (call1 < 0) ? LANC_ANC_UNASSIGNED : static_cast<uint8_t>(call1);

            const int32_t gt0 = gtArr[2 * s + 0];
            const int32_t gt1 = gtArr[2 * s + 1];

            if (bcf_gt_is_missing(gt0)) {
                miss0[j] = 1;
                alt0[j] = 0;
                ++nMissing;
            } else {
                const int a0 = bcf_gt_allele(gt0);
                if (a0 < 0 || a0 > 1)
                    throw std::runtime_error("--make-lanc: unexpected allele code at biallelic site " +
                                             std::string(chromStr) + ":" + std::to_string(pos0 + 1) + " in " +
                                             genoPath);
                miss0[j] = 0;
                alt0[j] = static_cast<uint8_t>(a0);
            }
            if (bcf_gt_is_missing(gt1)) {
                miss1[j] = 1;
                alt1[j] = 0;
                ++nMissing;
            } else {
                const int a1 = bcf_gt_allele(gt1);
                if (a1 < 0 || a1 > 1)
                    throw std::runtime_error("--make-lanc: unexpected allele code at biallelic site " +
                                             std::string(chromStr) + ":" + std::to_string(pos0 + 1) + " in " +
                                             genoPath);
                miss1[j] = 0;
                alt1[j] = static_cast<uint8_t>(a1);
            }
        }

        writer.addMarker(anc0.data(), anc1.data(), alt0.data(), alt1.data(), miss0.data(), miss1.data());

        bimOut << chromStr << "\t"
               << (rec->d.id && rec->d.id[0] != '.' ? rec->d.id
                                                    : std::string(chromStr) + ":" + std::to_string(pos0 + 1))
               << "\t0\t" << (pos0 + 1) << "\t" << rec->d.allele[0] << "\t" << rec->d.allele[1] << "\n";

        ++nMatch;
        if (nMatch % 10000 == 0) infoMsg("  Wrote %u markers (%s)", nMatch, genoPath.c_str());
    }

    std::free(gtArr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    nMissingOut = nMissing;
    return nMatch;
}

} // namespace

// ── Public API ────────────────────────────────────────────────────────────

void convertRfmixToLanc(
    const std::string &genoPrefix,
    bool expectBcf,
    const std::string &mspPrefix,
    const std::string &outPrefix,
    int compressionLevel,
    const std::string &keepFile,
    const std::string &removeFile,
    int nthreads
) {
    if (nthreads < 1) nthreads = 1;

    // ── 1. Glob + parse MSP files; derive each file's chromosome token ──
    auto mspPaths = globFiles(mspPrefix, {".msp.tsv"});
    if (mspPaths.empty())
        throw std::runtime_error("convertRfmixToLanc: no files matched " + mspPrefix + "*.msp.tsv");

    struct ChrEntry {
        std::string token;    // chromosome string with leading chr/CHR stripped
        std::string mspChrom; // literal chromosome as it appears in the msp file
        std::string mspPath;
        MspData msp;
        std::string genoPath;
        std::string genoChrom;
    };
    std::vector<ChrEntry> entries;
    entries.reserve(mspPaths.size());

    for (const auto &mp : mspPaths) {
        MspData msp = parseMsp(mp);
        const std::string &chrom0 = msp.windows.front().chrom;
        for (const auto &w : msp.windows) {
            if (w.chrom != chrom0)
                throw std::runtime_error("convertRfmixToLanc: " + mp + " mixes chromosomes ('" + chrom0 +
                                         "' and '" + w.chrom +
                                         "'); each --rfmix-msp file must be single-chromosome.");
        }
        ChrEntry e;
        e.mspChrom = chrom0;
        e.token = stripChrPrefix(chrom0);
        e.mspPath = mp;
        e.msp = std::move(msp);
        entries.push_back(std::move(e));
    }

    // ── 2. Cross-file K / sample-ID (list + order) consistency ──────────
    const int K = entries.front().msp.K;
    const std::vector<std::string> &sampleIDs0 = entries.front().msp.sampleIDs;
    for (size_t i = 1; i < entries.size(); ++i) {
        const MspData &m = entries[i].msp;
        if (m.K != K)
            throw std::runtime_error("convertRfmixToLanc: K mismatch between " + entries.front().mspPath + " (K=" +
                                     std::to_string(K) + ") and " + entries[i].mspPath + " (K=" +
                                     std::to_string(m.K) + ")");
        if (m.sampleIDs.size() != sampleIDs0.size())
            throw std::runtime_error("convertRfmixToLanc: subject-count mismatch between " +
                                     entries.front().mspPath + " (" + std::to_string(sampleIDs0.size()) +
                                     ") and " + entries[i].mspPath + " (" + std::to_string(m.sampleIDs.size()) +
                                     ")");
        for (size_t s = 0; s < sampleIDs0.size(); ++s) {
            if (m.sampleIDs[s] != sampleIDs0[s])
                throw std::runtime_error("convertRfmixToLanc: subject order/ID mismatch between " +
                                         entries.front().mspPath + " and " + entries[i].mspPath + " at index " +
                                         std::to_string(s) + ": '" + sampleIDs0[s] + "' vs '" + m.sampleIDs[s] +
                                         "'");
        }
    }

    // ── 3. Glob geno files; match each to an msp entry by token ─────────
    const std::vector<std::string> genoSuffixes =
        expectBcf ? std::vector<std::string>{".bcf"}
                  : std::vector<std::string>{".vcf.gz", ".vcf.zst", ".vcf"};
    auto genoPaths = globFiles(genoPrefix, genoSuffixes);
    if (genoPaths.empty())
        throw std::runtime_error("convertRfmixToLanc: no files matched " + genoPrefix + "* (" +
                                 (expectBcf ? "*.bcf" : "*.vcf.gz / *.vcf.zst / *.vcf") + ")");

    std::unordered_map<std::string, size_t> tokenToEntry;
    for (size_t i = 0; i < entries.size(); ++i) {
        if (!tokenToEntry.emplace(entries[i].token, i).second)
            throw std::runtime_error("convertRfmixToLanc: duplicate chromosome token '" + entries[i].token +
                                     "' across --rfmix-msp files (" + entries[i].mspPath + ")");
    }

    for (const auto &gp : genoPaths) {
        const std::string gchrom = detectGenoChrom(gp, expectBcf);
        const std::string tok = stripChrPrefix(gchrom);
        auto it = tokenToEntry.find(tok);
        if (it == tokenToEntry.end()) {
            infoMsg("[WARN] convertRfmixToLanc: %s (chromosome '%s') has no matching --rfmix-msp file; skipped",
                    gp.c_str(), gchrom.c_str());
            continue;
        }
        ChrEntry &e = entries[it->second];
        if (!e.genoPath.empty())
            throw std::runtime_error("convertRfmixToLanc: multiple genotype files match chromosome token '" +
                                     tok + "': " + e.genoPath + " and " + gp);
        e.genoPath = gp;
        e.genoChrom = gchrom;
    }

    for (const auto &e : entries) {
        if (e.genoPath.empty())
            throw std::runtime_error("convertRfmixToLanc: chromosome '" + e.mspChrom + "' (msp file " +
                                     e.mspPath + ") has no matching genotype file under prefix " + genoPrefix);
    }

    // ── 4. Sort chromosomes (numeric-aware, matches LancData / .bim order) ──
    std::sort(entries.begin(), entries.end(),
              [](const ChrEntry &a, const ChrEntry &b) { return chromLessThan(a.token, b.token); });

    // ── 5. --keep / --remove (shared cohort; order preserved) ───────────
    auto keptIndices = buildKeptIndices(sampleIDs0, keepFile, removeFile);
    if (keptIndices.empty() && (!keepFile.empty() || !removeFile.empty()))
        throw std::runtime_error("convertRfmixToLanc: no subjects remain after --keep/--remove filters");
    const uint32_t nKept =
        keptIndices.empty() ? static_cast<uint32_t>(sampleIDs0.size()) : static_cast<uint32_t>(keptIndices.size());

    // ── 6. Write the shared .fam once (filtered, original BCF/msp order) ──
    {
        std::ofstream fam(outPrefix + ".fam");
        if (!fam) throw std::runtime_error("Cannot create " + outPrefix + ".fam");
        for (uint32_t j = 0; j < nKept; ++j) {
            const uint32_t s = keptIndices.empty() ? j : keptIndices[j];
            const std::string &sid = sampleIDs0[s];
            fam << sid << "\t" << sid << "\t0\t0\t0\t-9\n";
        }
    }

    // ── 7. Per chromosome: stream genotypes + msp windows -> .lanc/.bim ──
    infoMsg("convertRfmixToLanc: %zu chromosome(s), K=%d, %u subjects (%u kept)", entries.size(), K,
            static_cast<uint32_t>(sampleIDs0.size()), nKept);

    uint64_t totalMissing = 0;
    uint64_t totalMarkers = 0;
    for (const auto &e : entries) {
        const std::string lancPath = outPrefix + ".chr" + e.token + ".lanc";
        const std::string bimPath = outPrefix + ".chr" + e.token + ".bim";

        infoMsg("Chromosome %s: msp=%s geno=%s -> %s", e.token.c_str(), e.mspPath.c_str(), e.genoPath.c_str(),
                lancPath.c_str());

        LancWriter writer(lancPath, static_cast<uint8_t>(K), nKept, LANC_DEFAULT_BLOCKLEN, compressionLevel);
        std::ofstream bimOut(bimPath);
        if (!bimOut) throw std::runtime_error("Cannot create " + bimPath);

        uint64_t nMissing = 0;
        const uint32_t nMatch =
            processChrToLanc(e.genoPath, expectBcf, e.msp, writer, bimOut, keptIndices, nMissing, nthreads);

        writer.close();
        bimOut.close();

        if (nMatch == 0)
            throw std::runtime_error("convertRfmixToLanc: no " + std::string(expectBcf ? "BCF" : "VCF") +
                                     " variants in " + e.genoPath + " fell inside MSP windows for chromosome " +
                                     e.token + ". Check chromosome naming and coordinate systems.");

        totalMissing += nMissing;
        totalMarkers += nMatch;
    }

    if (totalMissing > 0)
        infoMsg("[INFO] convertRfmixToLanc: %llu missing genotypes encountered across all chromosomes",
                static_cast<unsigned long long>(totalMissing));

    infoMsg("Conversion complete: %s.chr*.lanc (%llu markers total x %u subjects x %d ancestries, %zu chromosomes)",
            outPrefix.c_str(), static_cast<unsigned long long>(totalMarkers), nKept, K, entries.size());
}
