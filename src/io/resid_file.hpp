// resid_file.hpp — Residual file I/O
//
// File format (whitespace-delimited, no header):
//   col 1: subject ID   (string, required)
//   col 2: residual     (double, required) — null-model / martingale residual
//   col 3: weight       (double, optional, default 1.0) — sampling weight
//   col 4: indicator    (double, optional, default 0.0) — 0/1 event indicator
//
// Lines starting with '#' and blank lines are ignored.
// Methods that only need (subject, residual) can use 2-column files;
// weight defaults to 1.0 and indicator to 0.0 when absent.
#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

struct ResidData {
  std::vector<std::string> subjects;
  Eigen::VectorXd          residuals;   // null-model / martingale residuals
  Eigen::VectorXd          weights;     // sampling weights  (1.0 when absent)
  Eigen::VectorXd          indicator;   // 0/1 event indicator (0.0 when absent)
};

// Parse a resid file and return a ResidData.
// Throws std::runtime_error if the file cannot be opened or a required
// column is missing / non-numeric on any data line.
ResidData loadResidFile(const std::string& filename);

// ──────────────────────────────────────────────────────────────────────
// Multi-column residual matrix loader (for SPAsqr etc.)
//
// File format (whitespace-delimited):
//   [optional header starting with '#':  #IID  col1  col2  ...]
//   subject_id  val1  val2  ...  valK
//
// Returns subject IDs and an (N × K) residual matrix.
// Throws std::runtime_error on I/O or parse errors.
// ──────────────────────────────────────────────────────────────────────

struct ResidMatrixData {
  std::vector<std::string> subjects;
  Eigen::MatrixXd          residuals;   // N × K
};

ResidMatrixData loadResidMatrix(const std::string& filename);

