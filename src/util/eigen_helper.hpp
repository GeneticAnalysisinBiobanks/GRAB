// eigen_fast_utils.hpp
#pragma once
#include <Eigen/Dense>
#include <vector>

namespace eigen_fast {

// ============================================================
// take (fast version)
// v(idx)
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


// ============================================================
// rows (fast version)
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
// cols (fast version)
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
// find (fast)
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

    return idx.head(k);
}


// ============================================================
// diag_take (fast)
// ============================================================

template<typename Derived>
inline void diag_take(
    const Eigen::MatrixBase<Derived>& A,
    const Eigen::Ref<const Eigen::VectorXi>& idx,
    Eigen::Ref<Eigen::VectorXd> out)
{
    auto d = A.diagonal();
    const int n = idx.size();

    for (int i = 0; i < n; ++i)
        out[i] = d[idx[i]];
}


// ============================================================
// logical select (fast)
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

    return out.head(k);
}

} // namespace eigen_fast
