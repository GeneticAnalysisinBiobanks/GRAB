// data_matrix.hpp — Design (covariate) matrix I/O, projection, and generic matrix loader
//
// File format: whitespace-delimited; first column is subject ID (skipped).
// Header line starting with '#' (e.g. "#IID COV1 COV2 ...") is skipped.
// Blank lines and '#'-prefixed lines are ignored.
#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

class DesignMatrix {
public:
  // Construct from a pre-parsed design matrix (e.g. from SubjectData::design()).
  // Pre-computes the projection matrices tX and XinvXX.
  explicit DesignMatrix(const Eigen::MatrixXd& X);

  int nRows() const { return static_cast<int>(m_X.rows()); }
  int nCols() const { return static_cast<int>(m_X.cols()); }

  // X           — (N × p) design matrix
  const Eigen::MatrixXd& X()      const { return m_X; }
  // tX          — (p × N) = X'
  const Eigen::MatrixXd& tX()     const { return m_tX; }
  // XinvXX      — (N × p) = X (X'X)^{-1}
  const Eigen::MatrixXd& XinvXX() const { return m_XinvXX; }

  // Compute covariate-adjusted genotype:
  //   adjG = G - XinvXX * (tX * G)
  // Exploits sparsity: only accumulates columns where G[j] != 0.
  // The result is written into `adjG` which must be pre-allocated (size N).
  void adjustGenotype(const double* G, const uint32_t* nzIdx, int nNz,
                      Eigen::Ref<Eigen::VectorXd> adjG) const;

private:
  Eigen::MatrixXd m_X;       // N × p
  Eigen::MatrixXd m_tX;      // p × N  (column-major → cols = subjects)
  Eigen::MatrixXd m_XinvXX;  // N × p
};

// ======================================================================
// Lightweight numeric matrix loader
// ======================================================================

// Load a headerless whitespace-delimited numeric file into an Eigen matrix.
// Lines starting with '#' and blank lines are skipped.
// Throws std::runtime_error on I/O or parse errors.
Eigen::MatrixXd loadNumericMatrix(const std::string& filename);


