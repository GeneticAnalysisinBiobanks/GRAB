// subject_data.cpp — Unified per-subject data loader and .fam alignment

#include "io/subject_data.hpp"

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

    // .fam has 6 columns: FID IID ...  — extract column 2 (IID).
    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };

    /*FID*/ nextTok();
    std::string iid = nextTok();
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
    m_famIIDs(std::move(famIIDs))
{}


// ══════════════════════════════════════════════════════════════════════
// Generic IID + numeric-columns parser
// ══════════════════════════════════════════════════════════════════════

SubjectData::RawFile SubjectData::parseIIDFile(
    const std::string& filename, int expectCols)
{
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
    // Legacy mode: col 0 = IID, remaining cols = numeric values.
    // `line` / lineNo hold the first data line from Phase 1.
    auto parseLegacy = [&](const std::string& ln, uint32_t lno) {
      const char*       p   = ln.c_str();
      const char* const end = p + ln.size();
      while (p < end && (*p == ' ' || *p == '\t')) ++p;
      if (p >= end) return;

      const char* tokStart = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      rf.iids.emplace_back(tokStart, p);

      int colsThisRow = 0;
      while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t')) ++p;
        if (p >= end) break;
        char* next;
        double v = std::strtod(p, &next);
        if (next == p)
          throw std::runtime_error(
              filename + " line " + std::to_string(lno) +
              ": non-numeric value after IID");
        rf.vals.push_back(v);
        ++colsThisRow;
        p = next;
      }

      if (rf.iids.size() == 1) {
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

    // Parse the first data line (already in `line` from Phase 1)
    if (!line.empty())
      parseLegacy(line, lineNo);

    // Parse remaining lines
    while (std::getline(ifs, line)) {
      ++lineNo;
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      parseLegacy(line, lineNo);
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

  if (iidCol < 0)
    throw std::runtime_error(filename + ": no IID column found in header");
  if (covarCols.empty())
    throw std::runtime_error(filename + ": no covariate columns found in header");

  // ── Parse data rows ────────────────────────────────────────────────
  RawFile rf;
  rf.nCols = static_cast<int>(covarCols.size());
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

void SubjectData::loadEigenVecs(const std::string& filename) {
  if (m_hasRawEigen)
    throw std::runtime_error("loadEigenVecs already called");

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

  if (iidCol < 0)
    throw std::runtime_error(filename + ": no IID column found in header");
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


// ══════════════════════════════════════════════════════════════════════
// finalize — intersect, build bitmask, reorder
// ══════════════════════════════════════════════════════════════════════

void SubjectData::finalize() {
  if (m_finalized)
    throw std::runtime_error("SubjectData::finalize already called");
  if (!m_hasRawResid && !m_hasRawCovar && !m_hasRawEigen)
    throw std::runtime_error("SubjectData: no data files loaded before finalize");

  // ── Build IID index maps for each loaded file ──────────────────────
  auto buildMap = [](const std::vector<std::string>& iids) {
    std::unordered_map<std::string, uint32_t> m;
    m.reserve(iids.size());
    for (uint32_t i = 0; i < static_cast<uint32_t>(iids.size()); ++i)
      m.emplace(iids[i], i);
    return m;
  };

  std::unordered_map<std::string, uint32_t> residMap, covarMap, eigenMap;
  if (m_hasRawResid)  residMap  = buildMap(m_rawResid.iids);
  if (m_hasRawCovar)  covarMap  = buildMap(m_rawCovar.iids);
  if (m_hasRawEigen)  eigenMap  = buildMap(m_rawEigen.iids);

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
    if (m_hasRawEigen  && eigenMap.find(iid)  == eigenMap.end())  continue;
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

  // ── Eigenvectors ───────────────────────────────────────────────────
  if (m_hasRawEigen) {
    const int nc = m_rawEigen.nCols;
    const double* src = m_rawEigen.vals.data();
    m_PCs.resize(N, nc);
    for (Eigen::Index i = 0; i < N; ++i) {
      uint32_t r = rowIn(eigenMap, m_famIIDs[usedFamIndices[i]]);
      for (int c = 0; c < nc; ++c)
        m_PCs(i, c) = src[r * nc + c];
    }
  }

  // ── Free raw storage ───────────────────────────────────────────────
  m_rawResid  = RawFile{};
  m_rawCovar   = RawFile{};
  m_rawEigen  = RawFile{};

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
