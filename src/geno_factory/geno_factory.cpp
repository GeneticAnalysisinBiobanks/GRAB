// geno_factory.cpp — Multi-format genotype data factory
//
// parseGenoIIDs():  extract sample IIDs from any supported format
// makeGenoData():   construct the right GenoMeta backend

#include "geno_factory/bgen.hpp"
#include "geno_factory/geno_data.hpp"
#include "geno_factory/pgen.hpp"
#include "geno_factory/plink.hpp"
#include "geno_factory/vcf.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/text_scanner.hpp"

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}

// ── parsePsamIIDs ────────────────────────────────────────────────────
// .psam is tab-delimited, header starts with #FID or #IID.
// Without a header, columns are FID IID PAT MAT SEX PHENO(s) (same as .fam).
static std::vector<std::string> parsePsamIIDs(const std::string &psamFile) {
    std::ifstream ifs(psamFile);
    if (!ifs) throw std::runtime_error("Cannot open " + psamFile);
    std::vector<std::string> iids;
    std::string line;
    int iidCol = 1; // default: column 1 (0-indexed) like .fam
    bool headerSeen = false;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '#' && line.size() > 1 && line[1] == '#') continue; // ## comment
        text::TokenScanner ts(line);
        std::vector<std::string> toks;
        while (!ts.atEnd()) {
            auto sv = ts.nextView();
            if (!sv.empty()) toks.emplace_back(sv);
        }
        if (toks.empty()) continue;
        if (!headerSeen && (toks[0] == "#FID" || toks[0] == "FID" || toks[0] == "#IID" || toks[0] == "IID")) {
            headerSeen = true;
            if (toks[0] == "#IID" || toks[0] == "IID")
                iidCol = 0;
            else
                iidCol = 1;
            continue;
        }
        headerSeen = true;
        if (static_cast<int>(toks.size()) <= iidCol) throw std::runtime_error("Too few columns in " + psamFile);
        iids.push_back(toks[iidCol]);
    }
    return iids;
}

// ── parseVcfSampleIDs ────────────────────────────────────────────────
static std::vector<std::string> parseVcfSampleIDs(const std::string &vcfFile) {
    htsFile *fp = hts_open(vcfFile.c_str(), "r");
    if (!fp) throw std::runtime_error("Cannot open " + vcfFile);
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        throw std::runtime_error("Cannot read header: " + vcfFile);
    }
    int n = bcf_hdr_nsamples(hdr);
    std::vector<std::string> ids;
    ids.reserve(n);
    for (int i = 0; i < n; ++i)
        ids.emplace_back(hdr->samples[i]);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    return ids;
}

// ── parseBgenSampleIDs ───────────────────────────────────────────────
static std::vector<std::string> parseBgenSampleIDs(const std::string &bgenFile) {
    std::ifstream ifs(bgenFile, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open " + bgenFile);
    // Read offset (first 4 bytes)
    uint32_t offset = 0;
    ifs.read(reinterpret_cast<char *>(&offset), 4);
    // Read header
    uint32_t headerLen = 0;
    ifs.read(reinterpret_cast<char *>(&headerLen), 4);
    uint32_t nVariants = 0, nSamples = 0;
    ifs.read(reinterpret_cast<char *>(&nVariants), 4);
    ifs.read(reinterpret_cast<char *>(&nSamples), 4);
    // Skip magic + free area, read flags
    ifs.seekg(4 + headerLen); // past offset + header
    uint32_t flags = 0;
    ifs.seekg(4 + headerLen - 4);
    ifs.read(reinterpret_cast<char *>(&flags), 4);
    bool hasSampleIds = (flags & 0x80000000u) != 0;
    // Position after header
    ifs.seekg(4 + headerLen);
    std::vector<std::string> ids;
    if (hasSampleIds) {
        uint32_t blockLen = 0, nIds = 0;
        ifs.read(reinterpret_cast<char *>(&blockLen), 4);
        ifs.read(reinterpret_cast<char *>(&nIds), 4);
        ids.reserve(nIds);
        for (uint32_t i = 0; i < nIds; ++i) {
            uint16_t idLen = 0;
            ifs.read(reinterpret_cast<char *>(&idLen), 2);
            std::string id(idLen, '\0');
            ifs.read(id.data(), idLen);
            ids.push_back(std::move(id));
        }
    } else {
        // No sample IDs in BGEN — try companion .sample file (Oxford format)
        std::string sampleFile = bgenFile;
        auto dotPos = sampleFile.rfind('.');
        if (dotPos != std::string::npos) sampleFile = sampleFile.substr(0, dotPos);
        sampleFile += ".sample";
        std::ifstream sfs(sampleFile);
        if (sfs) {
            // Oxford .sample: line 1 = header, line 2 = type row, data lines = FID IID ...
            std::string line;
            std::getline(sfs, line); // header (e.g. "ID_1 ID_2 missing")
            std::getline(sfs, line); // type row (e.g. "0 0 0")
            ids.reserve(nSamples);
            while (std::getline(sfs, line)) {
                if (line.empty()) continue;
                // IID is the second whitespace-delimited token
                text::TokenScanner ts(line);
                std::string fid = ts.next();
                std::string iid = ts.next();
                if (iid.empty()) iid = fid; // single-column fallback
                ids.push_back(std::move(iid));
            }
            if (ids.size() != nSamples)
                throw std::runtime_error(sampleFile + ": expected " + std::to_string(nSamples) +
                                         " subjects but found " + std::to_string(ids.size()));
            infoMsg("BGEN: read %u subject IDs from %s", nSamples, sampleFile.c_str());
        } else {
            // No .sample file either — error out
            throw std::runtime_error("BGEN file '" + bgenFile +
                                     "' has no sample identifier block "
                                     "and no companion .sample file found (tried '" +
                                     sampleFile + "'). Provide a .sample file alongside the BGEN file.");
        }
    }
    return ids;
}

// ── parseGenoIIDs ────────────────────────────────────────────────────

std::vector<std::string> parseGenoIIDs(const GenoSpec &spec) {
    switch (spec.format) {
    case GenoFormat::Plink:
        return parseFamIIDs(spec.path + ".fam");
    case GenoFormat::Pgen:
        return parsePsamIIDs(spec.path + ".psam");
    case GenoFormat::Vcf:
        return parseVcfSampleIDs(spec.path);
    case GenoFormat::Bgen:
        return parseBgenSampleIDs(spec.path);
    }
    throw std::runtime_error("Unknown genotype format");
}

// ── makeGenoData ─────────────────────────────────────────────────────

std::unique_ptr<GenoMeta> makeGenoData(
    const GenoSpec &spec,
    const std::vector<uint64_t> &usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    int nMarkersEachChunk
) {
    switch (spec.format) {
    case GenoFormat::Plink:
        infoMsg("Loading PLINK data: %s", spec.path.c_str());
        return std::make_unique<PlinkData>(spec.path + ".bed", spec.path + ".bim", spec.path + ".fam", usedMask,
                                           nSamplesInFile, nUsed, spec.extractFile, std::string{}, spec.excludeFile,
                                           std::string{}, spec.chrFilter, nMarkersEachChunk);
    case GenoFormat::Pgen:
        infoMsg("Loading PGEN data: %s", spec.path.c_str());
        return std::make_unique<PgenData>(spec.path + ".pgen", spec.path + ".pvar", usedMask, nSamplesInFile, nUsed,
                                          spec.chrFilter, nMarkersEachChunk);
    case GenoFormat::Vcf:
        infoMsg("Loading VCF/BCF data: %s", spec.path.c_str());
        return std::make_unique<VcfData>(spec.path, usedMask, nSamplesInFile, nUsed, spec.chrFilter,
                                         nMarkersEachChunk);
    case GenoFormat::Bgen:
        infoMsg("Loading BGEN data: %s", spec.path.c_str());
        return std::make_unique<BgenData>(spec.path, usedMask, nSamplesInFile, nUsed, spec.chrFilter,
                                          nMarkersEachChunk);
    }
    throw std::runtime_error("Unknown genotype format");
}
