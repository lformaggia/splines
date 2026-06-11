/**
 * @file tridiagonal.hpp
 * @brief Two-phase Thomas algorithm for tridiagonal linear systems.
 *
 * Implements the classic Thomas (LU-factorisation + back-substitution)
 * algorithm for a tridiagonal system of the form:
 *
 * @code
 *   | b[0] c[0]                    |   | x[0]   |   | d[0]   |
 *   | a[1] b[1] c[1]               |   | x[1]   |   | d[1]   |
 *   |      a[2] b[2] c[2]          | × | x[2]   | = | d[2]   |
 *   |           ...                |   | ...    |   | ...    |
 *   |               a[n-1] b[n-1]  |   | x[n-1] |   | d[n-1] |
 * @endcode
 *
 * The two-phase split matches the Fortran convention from `trif` (factor)
 * and `tris` / `tris2` (solve):
 *   1. Call `tridiag_factor(a, b, c)` — modifies `b` and `c` in place.
 *   2. Call `tridiag_solve(a, b, c, d)` — overwrites `d` with the solution.
 *
 * The solution type `T` in `tridiag_solve` may be a scalar (`Real`) or a
 * `Vec<Dim,Real>`, allowing simultaneous solution for all coordinate
 * dimensions of a spline tangent system.
 *
 * @note These functions reside in the `splines::detail` namespace; they are
 *       implementation helpers and are not part of the public API.
 */
#pragma once
#include <cassert>
#include <cmath>
#include <span>
#include <vector>
#include "vec.hpp"

namespace splines::detail {

/**
 * @brief In-place LU factorisation of a tridiagonal matrix (Thomas algorithm).
 *
 * After the call:
 * - `b[i]` holds the modified (pivot-reduced) diagonal element.
 * - `c[i]` holds `c[i] / b[i]` (pre-divided for use in back-substitution).
 * - `a` is read but not modified.
 *
 * @tparam Real  Floating-point scalar type (satisfies `RealType`).
 *
 * @param a  Sub-diagonal (`a[0]` is unused).
 * @param b  Main diagonal (modified in place).
 * @param c  Super-diagonal (`c[n-1]` is unused; modified in place).
 *
 * @pre  `a.size() == b.size() == c.size() >= 1`.
 * @pre  The matrix must be non-singular (non-zero pivots after reduction).
 */
template<RealType Real>
void tridiag_factor(std::span<Real> a, std::span<Real> b, std::span<Real> c) {
    const int n = static_cast<int>(b.size());
    assert(n >= 1);
    if (n == 1) return;

    c[0] /= b[0];
    for (int i = 1; i < n - 1; ++i) {
        b[i] -= a[i] * c[i - 1];
        c[i] /= b[i];
    }
    b[n - 1] -= a[n - 1] * c[n - 2];
}

/**
 * @brief In-place back-substitution for a factorised tridiagonal system.
 *
 * Solves `A x = d` using the LU factors produced by `tridiag_factor`.
 * On exit `d[i]` contains `x[i]`.
 *
 * The solution type `T` may be `Real` for scalar right-hand sides or
 * `Vec<Dim,Real>` to solve for all spatial dimensions simultaneously
 * (which is the typical usage when computing spline tangents).
 *
 * @tparam Real  Floating-point scalar type (satisfies `RealType`).
 * @tparam T     Element type of the right-hand-side vector: `Real` or
 *               `Vec<Dim,Real>`.  Must support `operator*`, `operator-`,
 *               and `operator*=(Real)`.
 *
 * @param a  Sub-diagonal (read only; as factorised).
 * @param b  Modified diagonal from `tridiag_factor` (read only).
 * @param c  Modified super-diagonal from `tridiag_factor` (`c[i] = c[i]/b[i]`).
 * @param d  Right-hand side on entry; solution on exit (modified in place).
 *
 * @pre  `a`, `b`, `c` must be the output of `tridiag_factor` for the same matrix.
 * @pre  `d.size() == b.size()`.
 */
template<RealType Real, typename T>
void tridiag_solve(std::span<const Real> a,
                   std::span<const Real> b,
                   std::span<const Real> c,
                   std::span<T> d) {
    const int n = static_cast<int>(d.size());
    assert(n >= 1);

    // Forward substitution (eliminate sub-diagonal)
    d[0] = d[0] * (Real(1) / b[0]);
    for (int i = 1; i < n; ++i)
        d[i] = (d[i] - a[i] * d[i - 1]) * (Real(1) / b[i]);

    // Back substitution (c[i] already stores c[i]/b[i] from factorisation)
    for (int i = n - 2; i >= 0; --i)
        d[i] = d[i] - c[i] * d[i + 1];
}

} // namespace splines::detail
