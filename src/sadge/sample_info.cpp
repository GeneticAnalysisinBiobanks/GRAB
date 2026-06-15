// sample_info.cpp — SADGE pedigree parsing + family-structure model.
#include "sadge/sample_info.hpp"

#include "util/text_scanner.hpp"

#include <fstream>
#include <stdexcept>
#include <unordered_set>

namespace sadge {

SampleInfo buildSampleInfo(
    const std::string &famFile,
    const std::vector<std::string> &genoIIDs
) {
    std::unordered_set<std::string> genoSet;
    genoSet.reserve(genoIIDs.size() * 2);
    for (const auto &id : genoIIDs)
        genoSet.insert(id);

    // ── Read the 4-column .fam and apply the R filters ────────────────────
    // Keep row iff IID ∈ geno and every non-"0" parent ∈ geno
    // (pipeline1.R filters FatherID/MotherID ∈ {"0"} ∪ geno).
    struct Row {
        std::string fid, iid, fa, mo;
    };
    std::vector<Row> kept;
    {
        std::ifstream in(famFile);
        if (!in.is_open())
            throw std::runtime_error("SADGE: cannot open --sadge-fam file: " + famFile);
        std::string line;
        size_t lineNo = 0;
        while (std::getline(in, line)) {
            ++lineNo;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            text::TokenScanner tok(line);
            Row r;
            r.fid = tok.next();
            r.iid = tok.next();
            r.fa = tok.next();
            r.mo = tok.next();
            if (r.iid.empty() || r.mo.empty())
                throw std::runtime_error(
                    "SADGE: --sadge-fam needs >=4 columns (FID IID FatherID MotherID) on line " +
                    std::to_string(lineNo));
            if (!genoSet.count(r.iid)) continue;
            const bool faOk = (r.fa == "0") || genoSet.count(r.fa) > 0;
            const bool moOk = (r.mo == "0") || genoSet.count(r.mo) > 0;
            if (!faOk || !moOk) continue;
            kept.push_back(std::move(r));
        }
    }

    // ── Group kept rows by FID, preserving first-appearance order ─────────
    std::unordered_map<std::string, int> famIndexOf;
    std::vector<std::vector<int> > famRows; // family index → kept row indices
    famIndexOf.reserve(kept.size());
    for (int i = 0; i < static_cast<int>(kept.size()); ++i) {
        auto it = famIndexOf.find(kept[i].fid);
        int fi;
        if (it == famIndexOf.end()) {
            fi = static_cast<int>(famRows.size());
            famIndexOf.emplace(kept[i].fid, fi);
            famRows.emplace_back();
        } else {
            fi = it->second;
        }
        famRows[fi].push_back(i);
    }

    // ── Pedigree subject set: valid-family siblings + their genotyped parents.
    std::unordered_set<std::string> pedSet;
    for (const auto &idxs : famRows) {
        if (idxs.size() != 1 && idxs.size() != 2) continue; // n_sib ∈ {1,2}
        for (int i : idxs) {
            pedSet.insert(kept[i].iid);
            if (kept[i].fa != "0") pedSet.insert(kept[i].fa);
            if (kept[i].mo != "0") pedSet.insert(kept[i].mo);
        }
    }

    SampleInfo si;
    // Decode-row order = genotype-file order among pedigree subjects.  This
    // matches the row order SADGE's own GenoCursor returns for the pedigree
    // usedMask built from the same genoIIDs.
    si.pedRowByIID.reserve(pedSet.size() * 2);
    for (const auto &iid : genoIIDs) {
        if (pedSet.count(iid)) {
            si.pedRowByIID.emplace(iid, static_cast<int32_t>(si.pedigreeIIDs.size()));
            si.pedigreeIIDs.push_back(iid);
        }
    }

    auto pedRow = [&](const std::string &iid) -> int32_t {
        auto it = si.pedRowByIID.find(iid);
        return (it == si.pedRowByIID.end()) ? -1 : it->second;
    };

    // ── Build per-sibling records (family index = stable grouping key) ────
    si.sibIndexByIID.reserve(kept.size() * 2);
    for (int fi = 0; fi < static_cast<int>(famRows.size()); ++fi) {
        const auto &idxs = famRows[fi];
        const int nSib = static_cast<int>(idxs.size());
        if (nSib != 1 && nSib != 2) continue;
        for (int k = 0; k < nSib; ++k) {
            const Row &row = kept[idxs[k]];
            const int nPar = (row.fa != "0" ? 1 : 0) + (row.mo != "0" ? 1 : 0);
            SibInfo s;
            s.iid = row.iid;
            s.fid = fi;
            s.genoNPar = nPar;
            s.genoNSib = nSib;
            s.fs1 = famStruct(nPar, nSib);
            s.pedSelf = pedRow(row.iid);
            s.pedCoSib = (nSib == 2) ? pedRow(kept[idxs[1 - k]].iid) : -1;
            s.pedPar1 = -1;
            s.pedPar2 = -1;
            // Normalize: the lone parent of a 1-parent family lands in pedPar1.
            if (row.fa != "0" && row.mo != "0") {
                s.pedPar1 = pedRow(row.fa);
                s.pedPar2 = pedRow(row.mo);
            } else if (row.fa != "0") {
                s.pedPar1 = pedRow(row.fa);
            } else if (row.mo != "0") {
                s.pedPar1 = pedRow(row.mo);
            }
            si.sibIndexByIID.emplace(s.iid, static_cast<int>(si.siblings.size()));
            si.siblings.push_back(std::move(s));
        }
    }

    return si;
}

} // namespace sadge
