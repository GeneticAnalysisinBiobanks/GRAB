// eigen_helper.hpp — Eigen utility functions for genotype and statistical computing
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <cstdint>

namespace eigen_fast {

// ============================================================
// take: out[i] = v[idx[i]]         — indexed gather (Eigen VectorXi)
// ============================================================

template<typename Derived>
inline void take(
    const Eigen::MatrixBase<Derived>& v,
    const Eigen::Ref<const Eigen::VectorXi>& idx,
    Eigen::Ref<Eigen::VectorXd> out)
{
    const int n = idx.size();
    for (int i = 0; i < n; ++i)
        out[i] = v[idx[i]];
}

// take: out[i] = v[idx[i]]         — indexed gather (uint32_t, for mtPLINK)
// Matches the std::vector<uint32_t> indexForMissing type used by PlinkCursor.

template<typename Derived>
inline void take(
    const Eigen::MatrixBase<Derived>& v,
    const std::vector<uint32_t>& idx,
    Eigen::Ref<Eigen::VectorXd> out)
{
    const int n = static_cast<int>(idx.size());
    for (int i = 0; i < n; ++i)
        out[i] = v[static_cast<Eigen::Index>(idx[i])];
}


// ============================================================
// rows: out.row(i) = A.row(idx[i]) — indexed row gather
// ============================================================

template<typename Derived>
inline void rows(
    const Eigen::MatrixBase<Derived>& A,
    const Eigen::Ref<const Eigen::VectorXi>& idx,
    Eigen::Ref<Eigen::MatrixXd> out)
{
    const int n = idx.size();
    for (int i = 0; i < n; ++i)
        out.row(i) = A.row(idx[i]);
}


// ============================================================
// cols: out.col(i) = A.col(idx[i]) — indexed column gather
// ============================================================

template<typename Derived>
inline void cols(
    const Eigen::MatrixBase<Derived>& A,
    const Eigen::Ref<const Eigen::VectorXi>& idx,
    Eigen::Ref<Eigen::MatrixXd> out)
{
    const int n = idx.size();
    for (int i = 0; i < n; ++i)
        out.col(i) = A.col(idx[i]);
}


// ============================================================
// scatter_fill: v[idx[i]] = val    — indexed scalar fill
//
// Per-marker hot path: imputes missing genotypes with 2 * altFreq.
// Uses uint32_t to match std::vector<uint32_t> indexForMissing in PlinkCursor.
// ============================================================

inline void scatter_fill(
    Eigen::Ref<Eigen::VectorXd> v,
    const std::vector<uint32_t>& idx,
    double val)
{
    for (uint32_t i : idx)
        v[static_cast<Eigen::Index>(i)] = val;
}


// ============================================================
// scatter: v[idx[i]] = src[i]      — indexed scatter (write)
// Complement of take. Uses Eigen VectorXi indices.
// ============================================================

template<typename DerivedSrc>
inline void scatter(
    Eigen::Ref<Eigen::VectorXd> v,
    const Eigen::Ref<const Eigen::VectorXi>& idx,
    const Eigen::MatrixBase<DerivedSrc>& src)
{
    const int n = idx.size();
    for (int i = 0; i < n; ++i)
        v[idx[i]] = src[i];
}


// ============================================================
// flip_geno: v[i] = 2.0 - v[i]    — strand flip (ALT <-> REF)
//
// Per-marker hot path for altFreq > 0.5 markers.
// .array() expression dispatches to AVX2 vectorized subtraction
// via Eigen's vectorization backend when compiled with -march=native.
// ============================================================

inline void flip_geno(Eigen::Ref<Eigen::VectorXd> v)
{
    v.array() = 2.0 - v.array();
}


// ============================================================
// find_gt: indices i where v[i] > threshold
//
// .eval() forces materialization into a new VectorXi before the
// local idx buffer is released, avoiding dangling-reference UB.
// ============================================================

template<typename Derived>
inline Eigen::VectorXi find_gt(
    const Eigen::MatrixBase<Derived>& v,
    double threshold)
{
    const int n = v.size();
    Eigen::VectorXi idx(n);
    int k = 0;
    for (int i = 0; i < n; ++i)
        if (v[i] > threshold)
            idx[k++] = i;
    return idx.head(k).eval();
}


// ============================================================
// find_lt: indices i where v[i] < threshold
// ============================================================

template<typename Derived>
inline Eigen::VectorXi find_lt(
    const Eigen::MatrixBase<Derived>& v,
    double threshold)
{
    const int n = v.size();
    Eigen::VectorXi idx(n);
    int k = 0;
    for (int i = 0; i < n; ++i)
        if (v[i] < threshold)
            idx[k++] = i;
    return idx.head(k).eval();
}


// ============================================================
// diag_take: out[i] = diag(A)[idx[i]]
// ============================================================

template<typename Derived>
inline void diag_take(
    const Eigen::MatrixBase<Derived>& A,
    const Eigen::Ref<const Eigen::VectorXi>& idx,
    Eigen::Ref<Eigen::VectorXd> out)
{
    const auto d = A.diagonal();
    const int n = idx.size();
    for (int i = 0; i < n; ++i)
        out[i] = d[idx[i]];
}


// ============================================================
// select_gt: values v[i] where v[i] > threshold
// ============================================================

template<typename Derived>
inline Eigen::VectorXd select_gt(
    const Eigen::MatrixBase<Derived>& v,
    double threshold)
{
    const int n = v.size();
    Eigen::VectorXd out(n);
    int k = 0;
    for (int i = 0; i < n; ++i)
        if (v[i] > threshold)
            out[k++] = v[i];
    return out.head(k).eval();
}


// ============================================================
// select_lt: values v[i] where v[i] < threshold
// ============================================================

template<typename Derived>
inline Eigen::VectorXd select_lt(
    const Eigen::MatrixBase<Derived>& v,
    double threshold)
{
    const int n = v.size();
    Eigen::VectorXd out(n);
    int k = 0;
    for (int i = 0; i < n; ++i)
        if (v[i] < threshold)
            out[k++] = v[i];
    return out.head(k).eval();
}


// ============================================================
// join_ones_left: [ones(n,1) | A]
//
// Replaces arma::join_horiz(arma::ones(N), PCs) in SPAmix/SPAmixPlus.
// Builds the design matrix with an intercept column prepended.
// ============================================================

inline Eigen::MatrixXd join_ones_left(const Eigen::MatrixXd& A)
{
    const Eigen::Index n = A.rows();
    const Eigen::Index p = A.cols();
    Eigen::MatrixXd out(n, 1 + p);
    out.col(0).setOnes();
    out.rightCols(p) = A;
    return out;
}


// ============================================================
// inv_sympd: inverse of a symmetric positive-definite matrix
//
// Uses Cholesky (LLT) decomposition instead of raw .inverse().
// Replaces arma::inv(XTX) in SPAmix's covariate projection.
// Prefer this over .inverse() for all X'X-type matrices.
// ============================================================

inline Eigen::MatrixXd inv_sympd(const Eigen::MatrixXd& A)
{
    const Eigen::Index p = A.rows();
    return A.llt().solve(Eigen::MatrixXd::Identity(p, p));
}


// ============================================================
// sqrt_diag_of_inv_sympd: sqrt( diag( A^{-1} ) )
//
// Replaces sqrt(arma::inv(XTX).diag()) in SPAmix's constructor.
// Avoids forming the full inverse matrix.
//
// Derivation:
//   A = L L'  (Cholesky)
//   A^{-1} = L'^{-1} L^{-1}
//   diag(A^{-1})[i] = || col_i(L^{-1}) ||^2
//   → sqrt gives the SE-scaling vector
// ============================================================

inline Eigen::VectorXd sqrt_diag_of_inv_sympd(const Eigen::MatrixXd& A)
{
    const Eigen::Index p = A.rows();
    Eigen::LLT<Eigen::MatrixXd> llt(A);
    // Solve L * X = I  →  X = L^{-1}
    const Eigen::MatrixXd Linv =
        llt.matrixL().solve(Eigen::MatrixXd::Identity(p, p));
    // diag(A^{-1})[i] = ||col_i(Linv)||^2  →  sqrt for SE scaling
    return Linv.colwise().squaredNorm().transpose().cwiseSqrt();
}

} // namespace eigen_fast
