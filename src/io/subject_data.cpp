// subject_data.cpp — Unified per-subject data loader and .fam alignment

#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <unordered_map>
#include <unordered_set>


// ══════════════════════════════════════════════════════════════════════
// parseFamIIDs — extract IID (column 2) from .fam
// ══════════════════════════════════════════════════════════════════════

std::vector<std::string> parseFamIIDs(const std::string& famFile) {
  std::ifstream in(famFile);
  if (!in.is_open())
    throw std::runtime_error("Cannot open PLINK .fam file: " + famFile);

  std::vector<std::string> iids;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    text::TokenScanner tok(line);
    /*FID*/ tok.next();
    std::string iid = tok.next();
    if (iid.empty())
      throw std::runtime_error(".fam: missing IID on line " +
                               std::to_string(iids.size() + 1));
    iids.push_back(std::move(iid));
  }
  return iids;
}


// ══════════════════════════════════════════════════════════════════════
// SubjectData construction
// ══════════════════════════════════════════════════════════════════════

SubjectData::SubjectData(std::vector<std::string> famIIDs)
  : m_nFam(static_cast<uint32_t>(famIIDs.size())),
    m_famIIDs(std::move(famIIDs)
) {}


// ══════════════════════════════════════════════════════════════════════
// Generic IID + numeric-columns parser
// ══════════════════════════════════════════════════════════════════════

SubjectData::RawFile SubjectData::parseIIDFile(
    const std::string& filename, int expectCols
) {
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs)
    throw std::runtime_error("Cannot open file: " + filename);

  RawFile rf;
  rf.nCols = 0;
  rf.iids.reserve(65536);
  rf.vals.reserve(65536);

  uint32_t lineNo = 0;
  std::string line;

  // ── Phase 1: Skip ##-comments, detect plink2-style header ──────────
  // Header: first non-## line starting with #FID, FID, #IID, or IID.
  // If no header found, the first data line uses legacy layout
  // (col 0 = IID, remaining = numeric values).

  bool hasHeader = false;
  int  iidCol = -1;
  std::vector<int> valCols;
  int  nHeaderTokens = 0;

  while (std::getline(ifs, line)) {
    ++lineNo;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    // Skip ## extra-comment lines (plink2 convention)
    if (line.size() >= 2 && line[0] == '#' && line[1] == '#') continue;

    if (line[0] == '#') {
      // Single-# line: check for header tokens #FID / #IID
      std::vector<std::string> tokens;
      { std::istringstream ss(line); std::string t;
        while (ss >> t) tokens.push_back(t); }
      const auto& t0 = tokens[0];
      if (t0 == "#FID" || t0 == "#IID") {
        hasHeader = true;
        for (int c = 0; c < static_cast<int>(tokens.size()); ++c) {
          const auto& h = tokens[c];
          if (h == "IID" || h == "#IID")      iidCol = c;
          else if (h == "#FID" || h == "FID") { /* skip */ }
          else { valCols.push_back(c); rf.colNames.push_back(h); }
        }
        if (iidCol < 0)
          throw std::runtime_error(filename +
              ": header found but no IID/#IID column");
        nHeaderTokens = static_cast<int>(tokens.size());
        rf.nCols = static_cast<int>(valCols.size());
        if (expectCols >= 0 && rf.nCols != expectCols)
          throw std::runtime_error(
              filename + ": expected " + std::to_string(expectCols) +
              " value columns, got " + std::to_string(rf.nCols));
        break;
      }
      // Not a header — skip as comment
      continue;
    }

    // Non-# line: check for FID / IID header without # prefix
    {
      const char* p0 = line.c_str();
      const char* e0 = p0 + line.size();
      while (p0 < e0 && (*p0 == ' ' || *p0 == '\t')) ++p0;
      const char* ts = p0;
      while (p0 < e0 && *p0 != ' ' && *p0 != '\t') ++p0;
      std::string firstTok(ts, p0);

      if (firstTok == "FID" || firstTok == "IID") {
        hasHeader = true;
        std::vector<std::string> tokens;
        { std::istringstream ss(line); std::string t;
          while (ss >> t) tokens.push_back(t); }
        for (int c = 0; c < static_cast<int>(tokens.size()); ++c) {
          const auto& h = tokens[c];
          if (h == "IID")      iidCol = c;
          else if (h == "FID") { /* skip */ }
          else { valCols.push_back(c); rf.colNames.push_back(h); }
        }
        if (iidCol < 0)
          throw std::runtime_error(filename +
              ": header found but no IID column");
        nHeaderTokens = static_cast<int>(tokens.size());
        rf.nCols = static_cast<int>(valCols.size());
        if (expectCols >= 0 && rf.nCols != expectCols)
          throw std::runtime_error(
              filename + ": expected " + std::to_string(expectCols) +
              " value columns, got " + std::to_string(rf.nCols));
        break;
      }
    }

    // Not a header — first data line (legacy mode)
    break;
  }

  // ── Phase 2: Parse data lines ──────────────────────────────────────

  if (hasHeader) {
    // Header mode: column positions known from Phase 1
    while (std::getline(ifs, line)) {
      ++lineNo;
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty()) continue;

      std::vector<std::string> tokens;
      { std::istringstream rss(line); std::string t;
        while (rss >> t) tokens.push_back(t); }

      if (static_cast<int>(tokens.size()) != nHeaderTokens)
        throw std::runtime_error(
            filename + " line " + std::to_string(lineNo) +
            ": expected " + std::to_string(nHeaderTokens) +
            " columns, got " + std::to_string(tokens.size()));

      rf.iids.push_back(tokens[iidCol]);
      for (int ci : valCols) {
        char* ep;
        double v = std::strtod(tokens[ci].c_str(), &ep);
        if (ep == tokens[ci].c_str())
          throw std::runtime_error(
              filename + " line " + std::to_string(lineNo) +
              ": non-numeric value in column " + std::to_string(ci + 1));
        rf.vals.push_back(v);
      }
    }
  } else {
    // No header found.  The file must be a pure numeric matrix
    // (all columns are values, no IID column).
    // `line` / lineNo hold the first data line from Phase 1.

    auto parseNumericRow = [&](const std::string& ln, uint32_t lno) {
      const char*       p   = ln.c_str();
      const char* const end = p + ln.size();
      int colsThisRow = 0;
      while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t')) ++p;
        if (p >= end) break;
        char* next;
        double v = std::strtod(p, &next);
        if (next == p)
          throw std::runtime_error(
              filename + " line " + std::to_string(lno) +
              ": non-numeric value");
        rf.vals.push_back(v);
        ++colsThisRow;
        p = next;
      }

      uint32_t rowIdx = static_cast<uint32_t>(rf.iids.size());
      rf.iids.push_back(std::to_string(rowIdx));  // synthetic row index

      if (rowIdx == 0) {
        if (expectCols >= 0 && colsThisRow != expectCols)
          throw std::runtime_error(
              filename + ": expected " + std::to_string(expectCols) +
              " numeric columns, got " + std::to_string(colsThisRow));
        rf.nCols = colsThisRow;
      } else {
        if (colsThisRow != rf.nCols)
          throw std::runtime_error(
              filename + " line " + std::to_string(lno) +
              ": expected " + std::to_string(rf.nCols) +
              " numeric columns, got " + std::to_string(colsThisRow));
      }
    };

    // Verify first data line is numeric; error out if not.
    if (!line.empty()) {
      const char* p = line.c_str();
      const char* end = p + line.size();
      while (p < end && (*p == ' ' || *p == '\t')) ++p;
      if (p < end) {
        char* numEnd;
        std::strtod(p, &numEnd);
        bool isNumeric = (numEnd != p) &&
                         (numEnd >= end || *numEnd == ' ' || *numEnd == '\t');
        if (!isNumeric)
          throw std::runtime_error(
              filename + ": no header (IID/#IID) found and first data value "
              "is not numeric. Please add a header line with an IID column.");
      }
      rf.noIID = true;
      parseNumericRow(line, lineNo);
    }

    // Parse remaining lines
    while (std::getline(ifs, line)) {
      ++lineNo;
      if (text::skipLine(line)) continue;
      parseNumericRow(line, lineNo);
    }
  }

  if (rf.iids.empty())
    throw std::runtime_error(filename + ": no data rows");

  return rf;
}


// ══════════════════════════════════════════════════════════════════════
// Residual loaders
// ══════════════════════════════════════════════════════════════════════

void SubjectData::loadResidOne(const std::string& filename) {
  if (m_hasRawResid)
    throw std::runtime_error("loadResid* already called");
  m_rawResid    = parseIIDFile(filename, -1);  // auto-detect columns
  m_residType   = ResidType::One;
  m_hasRawResid = true;
  m_residColNames = m_rawResid.colNames;
}

void SubjectData::loadResidWtCoxG(const std::string& filename) {
  if (m_hasRawResid)
    throw std::runtime_error("loadResid* already called");
  m_rawResid    = parseIIDFile(filename, 3);
  m_residType   = ResidType::WtCoxG;
  m_hasRawResid = true;
}

void SubjectData::loadResidSPAsqr(const std::string& filename) {
  if (m_hasRawResid)
    throw std::runtime_error("loadResid* already called");
  m_rawResid    = parseIIDFile(filename, -1);  // auto-detect K
  m_residType   = ResidType::SPAsqr;
  m_hasRawResid = true;
}


// ══════════════════════════════════════════════════════════════════════
// Optional per-subject file loaders
// ══════════════════════════════════════════════════════════════════════

void SubjectData::loadCovar(const std::string& filename) {
  if (m_hasRawCovar)
    throw std::runtime_error("loadCovar already called");

  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs)
    throw std::runtime_error("Cannot open covariate file: " + filename);

  // ── Read header (first non-blank line) ─────────────────────────────
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) break;
  }
  if (line.empty())
    throw std::runtime_error("Empty covariate file: " + filename);

  std::vector<std::string> headers;
  { std::istringstream hss(line); std::string tok;
    while (hss >> tok) headers.push_back(tok); }

  // ── Find IID column and covariate columns ──────────────────────────
  // Skip standard plink2 metadata columns: FID, IID, SID, PAT, MAT, SEX, PHENO*
  int iidCol = -1;
  std::vector<int> covarCols;

  for (int c = 0; c < static_cast<int>(headers.size()); ++c) {
    const auto& h = headers[c];
    if (h == "IID" || h == "#IID") {
      iidCol = c;
    } else if (h == "#FID" || h == "FID" || h == "SID" ||
               h == "PAT"  || h == "MAT" || h == "SEX") {
      // skip
    } else if (h.size() >= 5 && h.substr(0, 5) == "PHENO") {
      // skip phenotype columns
    } else {
      covarCols.push_back(c);
    }
  }

  if (iidCol < 0) {
    // No IID column — check if the first line was actually a data row
    // (pure numeric matrix with no header).  Re-parse from the beginning.
    ifs.clear();
    ifs.seekg(0);

    RawFile rf;
    rf.noIID = true;
    rf.iids.reserve(65536);
    rf.vals.reserve(65536);

    uint32_t lineNo2 = 0;
    std::string ln;
    while (std::getline(ifs, ln)) {
      ++lineNo2;
      if (!ln.empty() && ln.back() == '\r') ln.pop_back();
      if (ln.empty()) continue;

      const char*       p   = ln.c_str();
      const char* const end = p + ln.size();
      int colsThisRow = 0;
      while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t')) ++p;
        if (p >= end) break;
        char* next;
        double v = std::strtod(p, &next);
        if (next == p)
          throw std::runtime_error(
              filename + " line " + std::to_string(lineNo2) +
              ": non-numeric value in numeric matrix");
        rf.vals.push_back(v);
        ++colsThisRow;
        p = next;
      }

      uint32_t rowIdx = static_cast<uint32_t>(rf.iids.size());
      rf.iids.push_back(std::to_string(rowIdx));  // synthetic

      if (rowIdx == 0) {
        rf.nCols = colsThisRow;
      } else if (colsThisRow != rf.nCols) {
        throw std::runtime_error(
            filename + " line " + std::to_string(lineNo2) +
            ": expected " + std::to_string(rf.nCols) +
            " columns, got " + std::to_string(colsThisRow));
      }
    }
    if (rf.iids.empty())
      throw std::runtime_error(filename + ": no data rows");

    m_rawCovar    = std::move(rf);
    m_hasRawCovar = true;
    return;
  }
  if (covarCols.empty())
    throw std::runtime_error(filename + ": no covariate columns found in header");

  // ── Parse data rows ────────────────────────────────────────────────
  RawFile rf;
  rf.nCols = static_cast<int>(covarCols.size());
  for (int ci : covarCols)
    rf.colNames.push_back(headers[ci]);
  rf.iids.reserve(65536);
  rf.vals.reserve(65536);

  const int nExpectedToks = static_cast<int>(headers.size());
  uint32_t lineNo = 1; // header was line 1
  while (std::getline(ifs, line)) {
    ++lineNo;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    std::vector<std::string> tokens;
    { std::istringstream rss(line); std::string tok;
      while (rss >> tok) tokens.push_back(tok); }

    if (static_cast<int>(tokens.size()) != nExpectedToks)
      throw std::runtime_error(
          filename + " line " + std::to_string(lineNo) +
          ": expected " + std::to_string(nExpectedToks) +
          " columns, got " + std::to_string(tokens.size()));

    rf.iids.push_back(tokens[iidCol]);
    for (int ci : covarCols) {
      char* end;
      double v = std::strtod(tokens[ci].c_str(), &end);
      if (end == tokens[ci].c_str())
        throw std::runtime_error(
            filename + " line " + std::to_string(lineNo) +
            ": non-numeric covariate value: " + tokens[ci]);
      rf.vals.push_back(v);
    }
  }

  if (rf.iids.empty())
    throw std::runtime_error(filename + ": no data rows");

  m_rawCovar    = std::move(rf);
  m_hasRawCovar = true;
}

void SubjectData::loadPhenoFile(const std::string& filename) {
  if (m_hasRawPheno)
    throw std::runtime_error("loadPhenoFile already called");
  m_rawPheno    = parseIIDFile(filename, -1);
  m_hasRawPheno = true;
}

#if 0  // ── loadEigenVecs removed; replaced by loadPhenoFile ────────────────
void SubjectData::loadEigenVecs(const std::string& filename) {

  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs)
    throw std::runtime_error("Cannot open eigenvec file: " + filename);

  // ── Read header ────────────────────────────────────────────────────
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) break;
  }
  if (line.empty())
    throw std::runtime_error("Empty eigenvec file: " + filename);

  std::vector<std::string> headers;
  { std::istringstream hss(line); std::string tok;
    while (hss >> tok) headers.push_back(tok); }

  // ── Find IID column and PC columns ─────────────────────────────────
  int iidCol = -1;
  std::vector<int> pcCols;

  for (int c = 0; c < static_cast<int>(headers.size()); ++c) {
    const auto& h = headers[c];
    if (h == "IID" || h == "#IID") {
      iidCol = c;
    } else if (h == "#FID" || h == "FID" || h == "SID") {
      // skip
    } else if (h.size() >= 2 && (h[0] == 'P' || h[0] == 'p') &&
               (h[1] == 'C' || h[1] == 'c')) {
      pcCols.push_back(c);
    }
  }

  if (iidCol < 0) {
    // No IID column — treat as pure numeric matrix (no header).
    // Re-parse from beginning.
    ifs.clear();
    ifs.seekg(0);

    RawFile rf;
    rf.noIID = true;
    rf.iids.reserve(65536);
    rf.vals.reserve(65536);

    uint32_t lineNo2 = 0;
    std::string ln;
    while (std::getline(ifs, ln)) {
      ++lineNo2;
      if (!ln.empty() && ln.back() == '\r') ln.pop_back();
      if (ln.empty()) continue;

      const char*       p   = ln.c_str();
      const char* const end = p + ln.size();
      int colsThisRow = 0;
      while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t')) ++p;
        if (p >= end) break;
        char* next;
        double v = std::strtod(p, &next);
        if (next == p)
          throw std::runtime_error(
              filename + " line " + std::to_string(lineNo2) +
              ": non-numeric value in numeric matrix");
        rf.vals.push_back(v);
        ++colsThisRow;
        p = next;
      }

      uint32_t rowIdx = static_cast<uint32_t>(rf.iids.size());
      rf.iids.push_back(std::to_string(rowIdx));  // synthetic

      if (rowIdx == 0) {
        rf.nCols = colsThisRow;
      } else if (colsThisRow != rf.nCols) {
        throw std::runtime_error(
            filename + " line " + std::to_string(lineNo2) +
            ": expected " + std::to_string(rf.nCols) +
            " columns, got " + std::to_string(colsThisRow));
      }
    }
    if (rf.iids.empty())
      throw std::runtime_error(filename + ": no data rows");

    m_rawEigen    = std::move(rf);
    m_hasRawEigen = true;
    return;
  }
  if (pcCols.empty())
    throw std::runtime_error(filename + ": no PC columns found in header");

  // ── Parse data rows ────────────────────────────────────────────────
  RawFile rf;
  rf.nCols = static_cast<int>(pcCols.size());
  rf.iids.reserve(65536);
  rf.vals.reserve(65536);

  const int nExpectedToks = static_cast<int>(headers.size());
  uint32_t lineNo = 1;
  while (std::getline(ifs, line)) {
    ++lineNo;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    std::vector<std::string> tokens;
    { std::istringstream rss(line); std::string tok;
      while (rss >> tok) tokens.push_back(tok); }

    if (static_cast<int>(tokens.size()) != nExpectedToks)
      throw std::runtime_error(
          filename + " line " + std::to_string(lineNo) +
          ": expected " + std::to_string(nExpectedToks) +
          " columns, got " + std::to_string(tokens.size()));

    rf.iids.push_back(tokens[iidCol]);
    for (int ci : pcCols) {
      char* end;
      double v = std::strtod(tokens[ci].c_str(), &end);
      if (end == tokens[ci].c_str())
        throw std::runtime_error(
            filename + " line " + std::to_string(lineNo) +
            ": non-numeric PC value: " + tokens[ci]);
      rf.vals.push_back(v);
    }
  }

  if (rf.iids.empty())
    throw std::runtime_error(filename + ": no data rows");

  m_rawEigen    = std::move(rf);
  m_hasRawEigen = true;
}
#endif  // loadEigenVecs removed


// ══════════════════════════════════════════════════════════════════════
// finalize — intersect, build bitmask, reorder
// ══════════════════════════════════════════════════════════════════════

void SubjectData::finalize() {
  if (m_finalized)
    throw std::runtime_error("SubjectData::finalize already called");
  if (!m_hasRawResid && !m_hasRawCovar && !m_hasRawPheno)
    throw std::runtime_error("SubjectData: no data files loaded before finalize");

  // ── Handle noIID files: assign .fam IIDs if row count matches ────
  auto assignFamIIDs = [&](RawFile& rf, const char* label) {
    if (!rf.noIID) return;
    if (static_cast<uint32_t>(rf.iids.size()) != m_nFam)
      throw std::runtime_error(
          std::string(label) + ": pure numeric matrix has " +
          std::to_string(rf.iids.size()) + " rows but .fam has " +
          std::to_string(m_nFam) + " — cannot infer subject identity");
    rf.iids = m_famIIDs;  // assume same order as .fam
  };
  if (m_hasRawResid)  assignFamIIDs(m_rawResid, "--null-resid");
  if (m_hasRawCovar)  assignFamIIDs(m_rawCovar, "--covar");
  if (m_hasRawPheno)  assignFamIIDs(m_rawPheno, "--pheno");

  // Log assumed column layout for noIID files
  if (m_hasRawResid && m_rawResid.noIID) {
    const int nc = m_rawResid.nCols;
    switch (m_residType) {
    case ResidType::One:
      infoMsg("--null-resid: no IID column; assuming .fam order, %d residual column(s)", nc);
      break;
    case ResidType::WtCoxG:
      infoMsg("--null-resid: no IID column; assuming .fam order, columns: RESID WEIGHT INDICATOR");
      break;
    case ResidType::SPAsqr:
      infoMsg("--null-resid: no IID column; assuming .fam order, columns: R_tau1 .. R_tau%d", nc);
      break;
    default: break;
    }
  }
  if (m_hasRawCovar && m_rawCovar.noIID)
    infoMsg("--covar: no IID column; assuming .fam order, %d covariate column(s)", m_rawCovar.nCols);
  if (m_hasRawPheno && m_rawPheno.noIID)
    infoMsg("--pheno: no IID column; assuming .fam order, %d column(s)", m_rawPheno.nCols);

  // ── Build IID index maps for each loaded file ──────────────────────
  std::unordered_map<std::string, uint32_t> residMap, covarMap, phenoMap;
  if (m_hasRawResid)  residMap  = text::buildIIDMap(m_rawResid.iids);
  if (m_hasRawCovar)  covarMap  = text::buildIIDMap(m_rawCovar.iids);
  if (m_hasRawPheno)  phenoMap  = text::buildIIDMap(m_rawPheno.iids);

  // ── Intersect with .fam, build bitmask ─────────────────────────────
  const uint32_t nWords = (m_nFam + 63) / 64;
  m_usedMask.assign(nWords, 0ULL);

  // First pass: determine which .fam subjects are in the intersection
  std::vector<uint32_t> usedFamIndices;
  usedFamIndices.reserve(m_nFam);

  for (uint32_t f = 0; f < m_nFam; ++f) {
    const auto& iid = m_famIIDs[f];
    if (m_hasRawResid  && residMap.find(iid)  == residMap.end())  continue;
    if (m_hasRawCovar  && covarMap.find(iid)  == covarMap.end())  continue;
    if (m_hasRawPheno  && phenoMap.find(iid)  == phenoMap.end())  continue;
    // Subject passes all loaded-file checks
    m_usedMask[f / 64] |= 1ULL << (f % 64);
    usedFamIndices.push_back(f);
  }
  m_nUsed = static_cast<uint32_t>(usedFamIndices.size());

  if (m_nUsed == 0)
    throw std::runtime_error("SubjectData: no subjects remain after intersection with .fam");

  // ── Allocate and fill dense arrays in .fam order ───────────────────
  const Eigen::Index N = static_cast<Eigen::Index>(m_nUsed);

  // Helper: get the row index in a raw file for the given .fam IID
  // (lookup is guaranteed to succeed because we already checked membership).
  auto rowIn = [](const std::unordered_map<std::string, uint32_t>& m,
                  const std::string& iid) -> uint32_t {
    return m.find(iid)->second;
  };

  // ── Residuals ──────────────────────────────────────────────────────
  if (m_hasRawResid) {
    const int nc = m_rawResid.nCols;
    const double* src = m_rawResid.vals.data();

    switch (m_residType) {
      case ResidType::One: {
        if (nc == 1) {
          m_residuals.resize(N);
          for (Eigen::Index i = 0; i < N; ++i) {
            uint32_t r = rowIn(residMap, m_famIIDs[usedFamIndices[i]]);
            m_residuals[i] = src[r];
          }
        } else {
          m_residMatrix.resize(N, nc);
          for (Eigen::Index i = 0; i < N; ++i) {
            uint32_t r = rowIn(residMap, m_famIIDs[usedFamIndices[i]]);
            for (int c = 0; c < nc; ++c)
              m_residMatrix(i, c) = src[r * nc + c];
          }
          m_residuals = m_residMatrix.col(0);
        }
        m_nResidOneCols = nc;
        break;
      }
      case ResidType::WtCoxG: {
        m_residuals.resize(N);
        m_weights.resize(N);
        m_indicator.resize(N);
        for (Eigen::Index i = 0; i < N; ++i) {
          uint32_t r = rowIn(residMap, m_famIIDs[usedFamIndices[i]]);
          m_residuals[i] = src[r * 3 + 0];
          m_weights[i]   = src[r * 3 + 1];
          m_indicator[i] = src[r * 3 + 2];
        }
        break;
      }
      case ResidType::SPAsqr: {
        m_residMatrix.resize(N, nc);
        for (Eigen::Index i = 0; i < N; ++i) {
          uint32_t r = rowIn(residMap, m_famIIDs[usedFamIndices[i]]);
          for (int c = 0; c < nc; ++c)
            m_residMatrix(i, c) = src[r * nc + c];
        }
        break;
      }
      default: break;
    }
  }

  // ── Design matrix ──────────────────────────────────────────────────
  if (m_hasRawCovar) {
    const int nc = m_rawCovar.nCols;
    const double* src = m_rawCovar.vals.data();
    m_covar.resize(N, nc);
    for (Eigen::Index i = 0; i < N; ++i) {
      uint32_t r = rowIn(covarMap, m_famIIDs[usedFamIndices[i]]);
      for (int c = 0; c < nc; ++c)
        m_covar(i, c) = src[r * nc + c];
    }
  }

  // ── Phenotype data ─────────────────────────────────────────────────
  if (m_hasRawPheno) {
    const int nc = m_rawPheno.nCols;
    const double* src = m_rawPheno.vals.data();
    m_phenoData.resize(N, nc);
    for (Eigen::Index i = 0; i < N; ++i) {
      uint32_t r = rowIn(phenoMap, m_famIIDs[usedFamIndices[i]]);
      for (int c = 0; c < nc; ++c)
        m_phenoData(i, c) = src[r * nc + c];
    }
    m_phenoColNames = m_rawPheno.colNames;
    for (int c = 0; c < static_cast<int>(m_phenoColNames.size()); ++c)
      m_phenoColMap[m_phenoColNames[c]] = c;
  }

  // ── Build column name maps for covar ───────────────────────────────
  if (m_hasRawCovar) {
    m_covarColNames = m_rawCovar.colNames;
    for (int c = 0; c < static_cast<int>(m_covarColNames.size()); ++c)
      m_covarColMap[m_covarColNames[c]] = c;
  }

  // ── Free raw storage ───────────────────────────────────────────────
  m_rawResid  = RawFile{};
  m_rawCovar  = RawFile{};
  m_rawPheno  = RawFile{};

  m_finalized = true;
}


// ══════════════════════════════════════════════════════════════════════
// Named column accessors (valid after finalize)
// ══════════════════════════════════════════════════════════════════════

bool SubjectData::hasColumn(const std::string& name) const {
  return m_covarColMap.count(name) || m_phenoColMap.count(name);
}

Eigen::VectorXd SubjectData::getColumn(const std::string& name) const {
  auto it = m_covarColMap.find(name);
  if (it != m_covarColMap.end())
    return m_covar.col(it->second);
  auto it2 = m_phenoColMap.find(name);
  if (it2 != m_phenoColMap.end())
    return m_phenoData.col(it2->second);
  throw std::runtime_error(
      "SubjectData::getColumn: column '" + name + "' not found in covar or pheno data");
}

Eigen::MatrixXd SubjectData::getColumns(const std::vector<std::string>& names) const {
  const Eigen::Index N = static_cast<Eigen::Index>(m_nUsed);
  Eigen::MatrixXd mat(N, static_cast<Eigen::Index>(names.size()));
  for (size_t i = 0; i < names.size(); ++i)
    mat.col(static_cast<Eigen::Index>(i)) = getColumn(names[i]);
  return mat;
}

void SubjectData::fillColumnsInto(
    const std::vector<std::string>& names,
    Eigen::Ref<Eigen::MatrixXd> target) const {
  for (size_t i = 0; i < names.size(); ++i) {
    auto it = m_covarColMap.find(names[i]);
    if (it != m_covarColMap.end()) {
      target.col(static_cast<Eigen::Index>(i)) = m_covar.col(it->second);
      continue;
    }
    auto it2 = m_phenoColMap.find(names[i]);
    if (it2 != m_phenoColMap.end()) {
      target.col(static_cast<Eigen::Index>(i)) = m_phenoData.col(it2->second);
      continue;
    }
    throw std::runtime_error(
        "SubjectData::fillColumnsInto: column '" + names[i] + "' not found");
  }
}


// ══════════════════════════════════════════════════════════════════════
// setResidWeightIndicator — post-finalize setter for computed regression
// ══════════════════════════════════════════════════════════════════════

void SubjectData::setResidWeightIndicator(
    Eigen::VectorXd resid, Eigen::VectorXd weight, Eigen::VectorXd ind) {
  if (!m_finalized)
    throw std::runtime_error("setResidWeightIndicator: finalize must be called first");
  const auto N = static_cast<Eigen::Index>(m_nUsed);
  if (resid.size() != N || weight.size() != N || ind.size() != N)
    throw std::runtime_error(
        "setResidWeightIndicator: vector length (" +
        std::to_string(resid.size()) + ") does not match nUsed (" +
        std::to_string(m_nUsed) + ")");
  m_residuals  = std::move(resid);
  m_weights    = std::move(weight);
  m_indicator  = std::move(ind);
}


// ══════════════════════════════════════════════════════════════════════
// initFromMask — direct initialization from pre-computed bitmask
// ══════════════════════════════════════════════════════════════════════

void SubjectData::initFromMask(
    std::vector<uint64_t> mask, uint32_t nUsed,
    Eigen::VectorXd resid, Eigen::VectorXd weight, Eigen::VectorXd ind) {
  const auto N = static_cast<Eigen::Index>(nUsed);
  if (resid.size() != N || weight.size() != N || ind.size() != N)
    throw std::runtime_error(
        "initFromMask: vector length (" + std::to_string(resid.size()) +
        ") does not match nUsed (" + std::to_string(nUsed) + ")");
  m_usedMask  = std::move(mask);
  m_nUsed     = nUsed;
  m_residuals = std::move(resid);
  m_weights   = std::move(weight);
  m_indicator = std::move(ind);
  m_finalized = true;
}


// ══════════════════════════════════════════════════════════════════════
// usedIIDs — used IIDs in .fam order
// ══════════════════════════════════════════════════════════════════════

std::vector<std::string> SubjectData::usedIIDs() const {
  std::vector<std::string> out;
  out.reserve(m_nUsed);
  for (uint32_t f = 0; f < m_nFam; ++f) {
    if (m_usedMask[f / 64] & (1ULL << (f % 64)))
      out.push_back(m_famIIDs[f]);
  }
  return out;
}
