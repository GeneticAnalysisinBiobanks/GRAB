// abed_convert_msp.cpp — Convert phased VCF/BCF + MSP ancestry windows to .abed
//
// Algorithm:
//  1. Parse MSP: detect K from header, collect sample IDs, load all windows.
//  2. Write .fam from MSP sample IDs.
//  3. Two-pass over VCF (each pass opens the file independently via htslib):
//       Pass 1: count biallelic SNPs that fall inside an MSP window.
//       Pass 2: write .bim + 2*K .abed tracks per marker.
//  4. For each variant at vcf-side 0-based position p:
//       find MSP window with spos <= p < epos and matching chrom.
//       Per-ancestry k, per-sample s:
//         dosage[k][s]   = sum_{h=0,1} { (call[2s+h]==k) * alt_allele(gt[2s+h]) }
//         hapcount[k][s] = sum_{h=0,1} { call[2s+h]==k }
//       Tracks are written in order: anc0_dosage, anc0_hapcount, anc1_dosage, ...

#include "localplus/abed_convert_msp.hpp"
#include "localplus/abed_io.hpp"
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
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// ── simple line reader for gzFile ──────────────────────────────────────────

static bool mspReadLine(
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

static std::vector<std::string> mspSplitTabs(const std::string &s) {
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

// ── MSP data structures ─────────────────────────────────────────────────────

struct MspWindow {
    std::string chrom;
    int32_t spos;              // 0-based, inclusive
    int32_t epos;              // 0-based, exclusive
    std::vector<int8_t> calls; // length 2*nSamples; each is ancestry 0..K-1
};

struct MspData {
    int K;
    uint32_t nSamples;
    std::vector<std::string> sampleIDs;
    std::vector<MspWindow> windows;
};

// ── Parse MSP file ──────────────────────────────────────────────────────────

static MspData parseMsp(const std::string &mspFile) {
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
        // K = number of tab-separated entries in rest
        msp.K = 1;
        for (char c : rest)
            if (c == '\t') ++msp.K;
    }
    if (msp.K < 1) {
        gzclose(gz);
        throw std::runtime_error("MSP: could not detect K from line 1");
    }

    // Line 2: #chm\tspos\tepos\tsgpos\tegpos\tn snps\tIID0.0\tIID0.1\tIID1.0\tIID1.1\t...
    if (!mspReadLine(gz, line)) {
        gzclose(gz);
        throw std::runtime_error("MSP file has only 1 line: " + mspFile);
    }
    {
        // Strip leading "#" if present
        auto cols = mspSplitTabs(line);
        if (cols.size() < 7) {
            gzclose(gz);
            throw std::runtime_error("MSP header line 2 has fewer than 7 columns");
        }
        // columns 6+ are haplotype IDs: IID0.0 IID0.1 IID1.0 IID1.1 ...
        // Number of sample columns must be even
        size_t nHapCols = cols.size() - 6;
        if (nHapCols == 0 || nHapCols % 2 != 0) {
            gzclose(gz);
            throw std::runtime_error("MSP header line 2: expected even number of subject columns (got " +
                                     std::to_string(nHapCols) + ")");
        }
        msp.nSamples = static_cast<uint32_t>(nHapCols / 2);
        msp.sampleIDs.resize(msp.nSamples);
        for (uint32_t s = 0; s < msp.nSamples; ++s) {
            const auto &hapID = cols[6 + s * 2]; // e.g. "SAMPLE1.0"
            // Strip last .digit suffix
            auto dot = hapID.rfind('.');
            if (dot != std::string::npos) {
                msp.sampleIDs[s] = hapID.substr(0, dot);
            } else {
                msp.sampleIDs[s] = hapID;
            }
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
                                     std::to_string(expected));
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

    infoMsg("MSP: K=%d, %u subjects, %zu windows", msp.K, msp.nSamples, msp.windows.size());
    return msp;
}

// ── Window lookup helper ─────────────────────────────────────────────────────
// Returns index into msp.windows for the window covering (chrom, pos0based),
// or -1 if not found.  Uses a cursor that advances monotonically.

// ── Chromosome comparison helper ─────────────────────────────────────────────
// Compares chromosome names numerically when possible (chr1 < chr2 < ... < chr10)
// instead of lexicographically (where chr10 < chr2).
static bool chromLessThan(
    const std::string &a,
    const std::string &b
) {
    if (a == b) return false;
    // Strip leading "chr" or "Chr" or "CHR" prefix
    auto stripChr = [](const std::string &s) -> std::string {
        if (s.size() >= 3 && (s[0] == 'c' || s[0] == 'C') && (s[1] == 'h' || s[1] == 'H') &&
            (s[2] == 'r' || s[2] == 'R'))
            return s.substr(3);
        return s;
    };
    std::string sa = stripChr(a), sb = stripChr(b);
    // Try numeric comparison
    bool aNum = !sa.empty() && sa.find_first_not_of("0123456789") == std::string::npos;
    bool bNum = !sb.empty() && sb.find_first_not_of("0123456789") == std::string::npos;
    if (aNum && bNum) return std::stoi(sa) < std::stoi(sb);
    // Numeric before non-numeric (e.g., chr22 < chrX)
    if (aNum != bNum) return aNum;
    return sa < sb;
}

struct WindowCursor {
    const std::vector<MspWindow> &windows;
    size_t idx = 0;

    int findWindow(
        const std::string &chrom,
        int32_t pos0
    ) {
        // Advance past windows that end before or are on earlier chroms
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
        return -1; // gap between windows on the same chrom, or different chrom
    }

};

// ── Core logic: iterate VCF, match to MSP, write tracks ────────────────────

// countFlag: if true, only count matching markers; writer/bimOut may be null
// keptIndices: 0-based sample indices to include (empty = use all, for pass 1)
// nMissingOut: if non-null, receives count of missing genotypes encountered
// ── Core logic: single-pass VCF → .abed + .bim (with optional parallelism) ────

// keptIndices: 0-based sample indices to include (empty = use all)
// nthreads: number of compute threads (>=2 enables batch parallel encoding)
static uint32_t processVcf(
    const std::string &vcfFile,
    const MspData &msp,
    AbedWriter &writer,
    std::ofstream &bimOut,
    const std::vector<uint32_t> &keptIndices,
    uint64_t *nMissingOut = nullptr,
    int nthreads = 1
) {
    htsFile *fp = hts_open(vcfFile.c_str(), "r");
    if (!fp) throw std::runtime_error("Cannot open VCF: " + vcfFile);

    if (nthreads > 1) {
        int nDecompressThreads = std::min(nthreads, 4);
        hts_set_threads(fp, nDecompressThreads);
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        throw std::runtime_error("Cannot read VCF header: " + vcfFile);
    }

    bcf1_t *rec = bcf_init();

    int nSamplesVcf = bcf_hdr_nsamples(hdr);
    if ((uint32_t)nSamplesVcf != msp.nSamples) {
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        throw std::runtime_error("Subject count mismatch: VCF has " + std::to_string(nSamplesVcf) +
                                 " subjects, MSP has " + std::to_string(msp.nSamples));
    }
    for (int s = 0; s < nSamplesVcf; ++s) {
        const char *vcfID = hdr->samples[s];
        if (msp.sampleIDs[(uint32_t)s] != vcfID) {
            bcf_destroy(rec);
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("Subject name mismatch at index " + std::to_string(s) + ": VCF has '" + vcfID +
                                     "', MSP has '" + msp.sampleIDs[(uint32_t)s] + "'");
        }
    }

    int32_t *gtArr = nullptr;
    int ngtArr = 0;

    const int K = msp.K;
    const uint32_t Nall = msp.nSamples;
    const uint32_t Nout = keptIndices.empty() ? Nall : static_cast<uint32_t>(keptIndices.size());
    const uint64_t bpt = (static_cast<uint64_t>(Nout) + 3) / 4;

    WindowCursor cursor{msp.windows};
    uint32_t nMatch = 0;
    uint64_t nMissing = 0;

    // ── Sequential — single-threaded write ──────────────────────────
    if (nthreads <= 1) {
        bool phasingChecked = false;
        while (bcf_read(fp, hdr, rec) == 0) {
            bcf_unpack(rec, BCF_UN_ALL);
            if (rec->n_allele != 2) continue;

            const char *chromStr = bcf_hdr_id2name(hdr, rec->rid);
            int32_t pos0 = rec->pos;

            int winIdx = cursor.findWindow(chromStr, pos0);
            if (winIdx < 0) continue;

            bimOut << chromStr << "\t"
                   << (rec->d.id && rec->d.id[0] != '.' ? rec->d.id
                                                        : std::string(chromStr) + ":" + std::to_string(pos0 + 1))
                   << "\t0\t" << (pos0 + 1) << "\t" << rec->d.allele[0] << "\t" << rec->d.allele[1] << "\n";

            bcf_get_genotypes(hdr, rec, &gtArr, &ngtArr);
            if (!phasingChecked && ngtArr >= 2) {
                if (!bcf_gt_is_phased(gtArr[1]))
                    throw std::runtime_error("--make-abfile: VCF genotypes are unphased. "
                                             "A phased VCF (genotypes separated by '|') is required.");
                phasingChecked = true;
            }

            const auto &w = msp.windows[(size_t)winIdx];
            std::vector<int8_t> dosage(Nout), hapcount(Nout);
            for (int k = 0; k < K; ++k) {
                for (uint32_t j = 0; j < Nout; ++j) {
                    uint32_t s = keptIndices.empty() ? j : keptIndices[j];
                    int dos = 0, hap = 0;
                    for (int h = 0; h < 2; ++h) {
                        int8_t call = w.calls[2 * s + h];
                        if (call == (int8_t)k) {
                            ++hap;
                            int32_t gt = gtArr[2 * s + h];
                            if (bcf_gt_is_missing(gt)) {
                                ++nMissing;
                            } else {
                                dos += bcf_gt_allele(gt);
                            }
                        }
                    }
                    dosage[j] = static_cast<int8_t>(dos >= 0 && dos <= 2 ? dos : -1);
                    hapcount[j] = static_cast<int8_t>(hap >= 0 && hap <= 2 ? hap : -1);
                }
                auto pd = AbedWriter::encode(dosage);
                auto ph = AbedWriter::encode(hapcount);
                writer.writeTrack(pd.data(), bpt);
                writer.writeTrack(ph.data(), bpt);
            }
            ++nMatch;
            if (nMatch % 10000 == 0) infoMsg("  Wrote %u markers", nMatch);
        }
    }

    // ── Parallel — batch I/O then parallel compute then ordered write ─
    else {
        struct MarkerSlot {
            std::string bimLine;
            int winIdx = -1;
            std::vector<int32_t> gt;
            std::vector<std::vector<uint8_t> > tracks;
            uint64_t nMissing = 0;
        };

        const int batchCap = std::max(nthreads, 64);
        std::vector<MarkerSlot> batch(static_cast<size_t>(batchCap));
        for (auto &slot : batch) {
            slot.gt.reserve(2 * Nall);
            slot.tracks.resize(2 * K);
        }
        int batchUsed = 0;

        auto processSlice = [&](int from, int to) {
            std::vector<int8_t> dosage(Nout), hapcount(Nout);
            for (int i = from; i < to; ++i) {
                MarkerSlot &slot = batch[i];
                slot.nMissing = 0;
                const auto &w = msp.windows[(size_t)slot.winIdx];
                for (int k = 0; k < K; ++k) {
                    for (uint32_t j = 0; j < Nout; ++j) {
                        uint32_t s = keptIndices.empty() ? j : keptIndices[j];
                        int dos = 0, hap = 0;
                        for (int h = 0; h < 2; ++h) {
                            int8_t call = w.calls[2 * s + h];
                            if (call == (int8_t)k) {
                                ++hap;
                                int32_t gtv = slot.gt[2 * s + h];
                                if (bcf_gt_is_missing(gtv)) {
                                    ++slot.nMissing;
                                } else {
                                    dos += bcf_gt_allele(gtv);
                                }
                            }
                        }
                        dosage[j] = static_cast<int8_t>(dos >= 0 && dos <= 2 ? dos : -1);
                        hapcount[j] = static_cast<int8_t>(hap >= 0 && hap <= 2 ? hap : -1);
                    }
                    slot.tracks[2 * k] = AbedWriter::encode(dosage);
                    slot.tracks[2 * k + 1] = AbedWriter::encode(hapcount);
                }
            }
        };

        auto flushBatch = [&]() {
            if (batchUsed == 0) return;
            int nt = std::min(nthreads, batchUsed);
            int sz = (batchUsed + nt - 1) / nt;
            std::vector<std::thread> workers;
            workers.reserve(static_cast<size_t>(nt - 1));
            for (int t = 0; t < nt - 1; ++t) {
                int from = t * sz;
                int to = std::min(from + sz, batchUsed);
                workers.emplace_back(processSlice, from, to);
            }
            processSlice((nt - 1) * sz, batchUsed);
            for (auto &w : workers)
                w.join();
            for (int i = 0; i < batchUsed; ++i) {
                bimOut << batch[i].bimLine;
                for (const auto &track : batch[i].tracks)
                    writer.writeTrack(track.data(), bpt);
                nMissing += batch[i].nMissing;
            }
            batchUsed = 0;
        };

        bool phasingChecked = false;
        while (bcf_read(fp, hdr, rec) == 0) {
            bcf_unpack(rec, BCF_UN_ALL);
            if (rec->n_allele != 2) continue;

            const char *chromStr = bcf_hdr_id2name(hdr, rec->rid);
            int32_t pos0 = rec->pos;

            int winIdx = cursor.findWindow(chromStr, pos0);
            if (winIdx < 0) continue;

            bcf_get_genotypes(hdr, rec, &gtArr, &ngtArr);
            if (!phasingChecked && ngtArr >= 2) {
                if (!bcf_gt_is_phased(gtArr[1]))
                    throw std::runtime_error("--make-abfile: VCF genotypes are unphased. "
                                             "A phased VCF (genotypes separated by '|') is required.");
                phasingChecked = true;
            }

            MarkerSlot &slot = batch[batchUsed];
            {
                const char *id_cstr = (rec->d.id && rec->d.id[0] != '.') ? rec->d.id : nullptr;
                slot.bimLine.clear();
                slot.bimLine += chromStr;
                slot.bimLine += '\t';
                if (id_cstr)
                    slot.bimLine += id_cstr;
                else {
                    slot.bimLine += chromStr;
                    slot.bimLine += ':';
                    slot.bimLine += std::to_string(pos0 + 1);
                }
                slot.bimLine += "\t0\t";
                slot.bimLine += std::to_string(pos0 + 1);
                slot.bimLine += '\t';
                slot.bimLine += rec->d.allele[0];
                slot.bimLine += '\t';
                slot.bimLine += rec->d.allele[1];
                slot.bimLine += '\n';
            }
            slot.winIdx = winIdx;
            slot.gt.assign(gtArr, gtArr + 2 * static_cast<int>(Nall));
            ++batchUsed;
            ++nMatch;

            if (batchUsed == batchCap) {
                flushBatch();
                if (nMatch % 50000 < static_cast<uint32_t>(batchCap)) infoMsg("  Wrote %u markers", nMatch);
            }
        }
        flushBatch();
        if (nMatch > 0 && nMatch % 50000 >= static_cast<uint32_t>(batchCap)) infoMsg("  Wrote %u markers", nMatch);
    }

    std::free(gtArr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    if (nMissingOut) *nMissingOut = nMissing;
    return nMatch;
}

// ── Public API ──────────────────────────────────────────────────────────────

void convertVcfMspToAbed(
    const std::string &vcfFile,
    const std::string &mspFile,
    const std::string &outPrefix,
    const std::string &keepFile,
    const std::string &removeFile,
    int nthreads
) {
    if (nthreads < 1) nthreads = 1;
    // 1. Parse MSP
    MspData msp = parseMsp(mspFile);

    // 2. Apply --keep / --remove subject filtering
    auto keptIndices = buildKeptIndices(msp.sampleIDs, keepFile, removeFile);
    if (keptIndices.empty() && (!keepFile.empty() || !removeFile.empty()))
        throw std::runtime_error("convertVcfMspToAbed: no subjects remain after --keep/--remove filters");

    // Build filtered sample ID list for .fam output
    const uint32_t nKept = keptIndices.empty() ? msp.nSamples : static_cast<uint32_t>(keptIndices.size());
    std::vector<std::string> keptSampleIDs(nKept);
    if (keptIndices.empty()) {
        keptSampleIDs = msp.sampleIDs;
    } else {
        for (uint32_t j = 0; j < nKept; ++j)
            keptSampleIDs[j] = msp.sampleIDs[keptIndices[j]];
    }

    // 3. Write .fam (filtered)
    {
        std::ofstream fam(outPrefix + ".fam");
        if (!fam) throw std::runtime_error("Cannot create " + outPrefix + ".fam");
        for (const auto &sid : keptSampleIDs)
            fam << sid << "\t" << sid << "\t0\t0\t0\t-9\n";
    }

    // 4. Single-pass: write .abed + .bim
    //    The VCF→abed path never writes 0b01 (missing) codes — missing genotypes
    //    are silently imputed to ref (dosage=0 → 0b00).  Safe to set NO_MISSING.
    infoMsg("Writing .abed and .bim (%d thread%s)...", nthreads, nthreads > 1 ? "s" : "");
    AbedWriter writer(outPrefix + ".abed", static_cast<uint8_t>(msp.K), nKept,
                      /*noMissing=*/ true, nthreads);

    std::ofstream bimOut(outPrefix + ".bim");
    if (!bimOut) throw std::runtime_error("Cannot create " + outPrefix + ".bim");

    uint64_t nMissing = 0;
    uint32_t nMarkers = processVcf(vcfFile, msp, writer, bimOut, keptIndices, &nMissing, nthreads);

    writer.close();
    bimOut.close();

    if (nMarkers == 0)
        throw std::runtime_error("No VCF variants fell inside MSP windows. "
                                 "Check chromosome naming and coordinate systems.");

    if (nMissing > 0)
        infoMsg("[INFO] --make-abfile: %llu missing genotypes encountered (imputed to ref)",
                (unsigned long long)nMissing);

    infoMsg("Conversion complete: %s.abed (%u markers x %u subjects x %d ancestries)", outPrefix.c_str(), nMarkers,
            nKept, msp.K);
}
