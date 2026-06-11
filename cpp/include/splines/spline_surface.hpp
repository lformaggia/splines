/**
 * @file spline_surface.hpp
 * @brief Bicubic tensor-product spline surface with chord-length parametrisation.
 *
 * `SplineSurface<Dim, Real>` interpolates an n×m lattice of knots using a
 * bicubic Hermite tensor-product surface.  The surface is C² in both
 * parametric directions.
 *
 * @par Parameter space
 * - Global parameter U ∈ [0, n−1] corresponds to columns (u-direction).
 * - Global parameter V ∈ [0, m−1] corresponds to rows    (v-direction).
 * - Within patch (pᵢ, p_j): local (u, v) = (U−pᵢ, V−p_j) ∈ [0,1]².
 *
 * @par Construction (port of Fortran `getta3`)
 * Four spline-fitting passes produce the Hermite data at every lattice node:
 * 1. **Pass 1** — For each column i (u=const), fit a v-direction spline
 *    through the knots to obtain ∂r/∂v and chord lengths in the v-direction.
 * 2. **Pass 2** — For each row j (v=const), fit a u-direction spline to
 *    obtain ∂r/∂u and chord lengths in the u-direction.
 * 3. **Pass 3** — For the boundary rows j=0 and j=m−1, fit a u-direction
 *    spline through the ∂r/∂v values to obtain ∂²r/∂u∂v at the boundary.
 * 4. **Pass 4** — For each column i, fit a v-direction spline through the
 *    ∂r/∂u values with the boundary ∂²r/∂u∂v imposed, yielding ∂²r/∂u∂v
 *    at all interior nodes.
 *
 * @par Patch evaluation (port of Fortran `evapa2` + `gpsur`)
 * Each patch is represented as a bicubic polynomial in (u, v) with
 * coefficients A[i][j] for u^i * v^j, obtained by the transformation
 * A = C * Q * Cᵀ where C is the standard 4×4 Hermite-to-monomial matrix.
 *
 * @par Fortran origins
 * Port of `getta3` (construction), `evapa2` (coefficients), `gpsur` / `evsurg`
 * (evaluation).
 */
#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <span>
#include <vector>
#include "concepts.hpp"
#include "end_conditions.hpp"
#include "spline_curve.hpp"
#include "tridiagonal.hpp"
#include "vec.hpp"

namespace splines {

/**
 * @brief Bicubic tensor-product spline surface interpolating an n×m lattice.
 *
 * @tparam Dim   Spatial dimension (typically 3, but 2 is also supported).
 * @tparam Real  Floating-point scalar type (default `double`).
 *
 * @par Example
 * @code
 * #include "splines/spline_surface.hpp"
 * using namespace splines;
 *
 * // Build a 7×5 grid on the unit sphere
 * int n = 7, m = 5;
 * std::vector<Vec3d> lattice(n * m);
 * for (int j = 0; j < m; ++j)
 *     for (int i = 0; i < n; ++i) {
 *         double u = M_PI * j / (m - 1);
 *         double v = 2*M_PI * i / (n - 1);
 *         lattice[j*n+i] = {sin(u)*cos(v), sin(u)*sin(v), cos(u)};
 *     }
 *
 * auto surf = SplineSurface<3>::interpolate(lattice, n, m);
 * auto pt = surf.eval(1.5, 2.3);   // position + all partials + normal
 * @endcode
 */
template<int Dim, RealType Real = double>
class SplineSurface {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * @brief Build a bicubic spline surface interpolating an n×m lattice.
     *
     * Lattice data must be stored row-major: `lattice[j*n + i]` is the
     * knot at column i, row j (i.e., u-index i, v-index j).
     *
     * @param lattice  Input knot positions, length n*m.
     * @param n        Number of columns (u-direction).
     * @param m        Number of rows    (v-direction).
     * @param ec_u     End condition for u-direction splines (default `NotAKnot`).
     * @param ec_v     End condition for v-direction splines (default `NotAKnot`).
     * @return         Fully constructed `SplineSurface`.
     *
     * @pre  `n >= 2 && m >= 2`.
     * @pre  `lattice.size() == n * m`.
     */
    static SplineSurface interpolate(
        std::span<const Vec<Dim, Real>> lattice,
        int                             n,
        int                             m,
        EndCondition                    ec_u = EndCondition::NotAKnot,
        EndCondition                    ec_v = EndCondition::NotAKnot);

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// @brief Number of columns in the lattice (u-direction extent).
    [[nodiscard]] int n() const noexcept { return n_; }

    /// @brief Number of rows in the lattice (v-direction extent).
    [[nodiscard]] int m() const noexcept { return m_; }

    // -----------------------------------------------------------------------
    // Evaluation
    // -----------------------------------------------------------------------

    /**
     * @brief Surface point with position, first- and second-order partials.
     *
     * All partial derivatives are with respect to the global parameters
     * U and V (not scaled by chord lengths).
     */
    struct Point {
        Vec<Dim, Real> pos;              ///< r(U, V)         — position.
        Vec<Dim, Real> ru;               ///< ∂r/∂U           — u-tangent.
        Vec<Dim, Real> rv;               ///< ∂r/∂V           — v-tangent.
        Vec<Dim, Real> ruu;              ///< ∂²r/∂U²         — u-curvature direction.
        Vec<Dim, Real> rvv;              ///< ∂²r/∂V²         — v-curvature direction.
        Vec<Dim, Real> ruv;              ///< ∂²r/∂U∂V        — twist vector.

        /**
         * @brief Unit outward normal (available only for 3-D surfaces).
         * @return  (rᵤ × r_v) / |rᵤ × r_v|, or zero if degenerate.
         */
        Vec<3, Real> normal() const requires (Dim == 3) {
            auto n = cross(ru, rv);
            Real len = norm(n);
            return (len > Real(1e-20)) ? n / len : Vec<3, Real>{};
        }
    };

    /**
     * @brief Evaluate the surface at global parameters (U, V).
     *
     * Both U and V are clamped to their valid ranges before evaluation.
     *
     * @param U  Parameter in [0, n()-1].
     * @param V  Parameter in [0, m()-1].
     * @return   `Point` with position and all first- and second-order partials.
     */
    [[nodiscard]] Point eval(Real U, Real V) const;

private:
    int n_{}, m_{};

    /// Knot positions; indexed as `coor_[j*n_ + i]`.
    std::vector<Vec<Dim, Real>>               coor_;

    /**
     * @brief Hermite tangent data at each lattice node.
     *
     * `tanret_[j*n_ + i]` is an array of three vectors:
     * - `[0]` = ∂r/∂u  (chord-scaled u-tangent)
     * - `[1]` = ∂r/∂v  (chord-scaled v-tangent)
     * - `[2]` = ∂²r/∂u∂v  (chord-scaled twist)
     */
    std::vector<std::array<Vec<Dim,Real>, 3>> tanret_;

    /**
     * @brief Chord lengths at each lattice node.
     *
     * `choret_[j*n_ + i]` is an array of two values:
     * - `[0]` = u-chord at arc (i, j)→(i+1, j)
     * - `[1]` = v-chord at arc (i, j)→(i, j+1)
     */
    std::vector<std::array<Real, 2>>          choret_;

    /**
     * @brief Evaluate a single bicubic patch at local coordinates.
     *
     * @param pi  Patch column index in [0, n_-2].
     * @param pj  Patch row    index in [0, m_-2].
     * @param u   Local u parameter in [0, 1].
     * @param v   Local v parameter in [0, 1].
     * @return    `Point` with position and all partials.
     */
    [[nodiscard]] Point evalPatch(int pi, int pj, Real u, Real v) const;
};

// ============================================================================
// Implementation helpers
// ============================================================================

namespace detail {

/**
 * @brief Hermite-to-monomial basis matrix (4×4).
 *
 * Row i of CC maps the Hermite data vector
 * `[r(0), r(1), r'(0), r'(1)]` to the polynomial coefficient of u^i.
 *
 * Equivalent to the standard Hermite basis matrix; matches the column-major
 * Fortran data: `cc / 1,0,-3,2, 0,0,3,-2, 0,1,-2,1, 0,0,-1,1 /`.
 */
constexpr double CC[4][4] = {
    { 1,  0,  0,  0},   // u^0
    { 0,  0,  1,  0},   // u^1
    {-3,  3, -2, -1},   // u^2
    { 2, -2,  1,  1}    // u^3
};

/**
 * @brief Compute 4×4 patch polynomial coefficients A = CC * Q * CCᵀ.
 *
 * Given the 4×4 Hermite corner data matrix Q (see `evalPatch` for layout),
 * returns A[i][j] such that the scalar surface component is
 * @f$ r(u,v) = \sum_{i,j} A_{ij}\, u^i v^j @f$.
 *
 * @tparam Real  Floating-point scalar type.
 * @param Q  Hermite data matrix (rows = u DOF, cols = v DOF).
 * @return   Polynomial coefficient matrix A.
 */
template<RealType Real>
std::array<std::array<Real,4>,4>
patchCoeffs(const std::array<std::array<Real,4>,4>& Q) {
    std::array<std::array<Real,4>,4> T{}, A{};
    // T = CC * Q
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 4; ++k) {
            T[i][k] = Real(0);
            for (int p = 0; p < 4; ++p)
                T[i][k] += Real(CC[i][p]) * Q[p][k];
        }
    // A = T * CC^T
    for (int i = 0; i < 4; ++i)
        for (int l = 0; l < 4; ++l) {
            A[i][l] = Real(0);
            for (int k = 0; k < 4; ++k)
                A[i][l] += T[i][k] * Real(CC[l][k]);
        }
    return A;
}

} // namespace detail

// ============================================================================
// SplineSurface::interpolate  (port of Fortran getta3)
// ============================================================================
template<int Dim, RealType Real>
SplineSurface<Dim, Real>
SplineSurface<Dim, Real>::interpolate(
    std::span<const Vec<Dim, Real>> lattice,
    int n, int m,
    EndCondition ec_u, EndCondition ec_v)
{
    assert(n >= 2 && m >= 2);
    assert(static_cast<int>(lattice.size()) == n * m);

    SplineSurface surf;
    surf.n_ = n;
    surf.m_ = m;
    surf.coor_.assign(lattice.begin(), lattice.end());
    surf.tanret_.resize(n * m);
    surf.choret_.resize(n * m);

    auto idx = [n](int i, int j) { return j * n + i; };

    ImposedTangents<Dim, Real> imposed_none{};

    // Working arrays (reused per row/column)
    std::vector<Vec<Dim, Real>> q_work(std::max(n, m));

    // -------------------------------------------------------------------
    // Pass 1: For each u=const column i, fit v-direction spline
    //         → tanret[*][1] = r_v,  choret[*][1] = v-chord
    // -------------------------------------------------------------------
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j)
            q_work[j] = surf.coor_[idx(i, j)];

        auto taus = detail::buildTangents<Dim,Real>(
            std::span{q_work.data(), size_t(m)},
            detail::computeChords<Dim,Real>(std::span{q_work.data(), size_t(m)}),
            ec_v, imposed_none);

        auto chords = detail::computeChords<Dim,Real>(
            std::span{q_work.data(), size_t(m)});

        for (int j = 0; j < m; ++j)
            surf.tanret_[idx(i, j)][1] = taus[j];
        for (int j = 0; j < m - 1; ++j)
            surf.choret_[idx(i, j)][1] = chords[j];
    }

    // -------------------------------------------------------------------
    // Pass 2: For each v=const row j, fit u-direction spline
    //         → tanret[*][0] = r_u,  choret[*][0] = u-chord
    // -------------------------------------------------------------------
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i)
            q_work[i] = surf.coor_[idx(i, j)];

        auto taus = detail::buildTangents<Dim,Real>(
            std::span{q_work.data(), size_t(n)},
            detail::computeChords<Dim,Real>(std::span{q_work.data(), size_t(n)}),
            ec_u, imposed_none);

        auto chords = detail::computeChords<Dim,Real>(
            std::span{q_work.data(), size_t(n)});

        for (int i = 0; i < n; ++i)
            surf.tanret_[idx(i, j)][0] = taus[i];
        for (int i = 0; i < n - 1; ++i)
            surf.choret_[idx(i, j)][0] = chords[i];
    }

    // -------------------------------------------------------------------
    // Pass 3: For the two boundary rows j=0 and j=m-1, fit u-direction
    //         spline through r_v values → boundary r_uv values
    // -------------------------------------------------------------------
    for (int jb : {0, m - 1}) {
        // Extract r_v values along this boundary row
        for (int i = 0; i < n; ++i)
            q_work[i] = surf.tanret_[idx(i, jb)][1];

        // Use the u-chord lengths already computed for this row
        std::vector<Real> cs_u(n - 1);
        for (int i = 0; i < n - 1; ++i)
            cs_u[i] = surf.choret_[idx(i, jb)][0];

        auto taus = detail::buildTangents<Dim,Real>(
            std::span{q_work.data(), size_t(n)},
            std::span<const Real>{cs_u},
            ec_u, imposed_none);

        for (int i = 0; i < n; ++i)
            surf.tanret_[idx(i, jb)][2] = taus[i];
    }

    // -------------------------------------------------------------------
    // Pass 4: For each u=const column i, fit v-direction spline through
    //         r_u values with imposed r_uv at j=0 and j=m-1
    //         → r_uv everywhere
    // -------------------------------------------------------------------
    for (int i = 0; i < n; ++i) {
        // q_work = r_u values along column i
        for (int j = 0; j < m; ++j)
            q_work[j] = surf.tanret_[idx(i, j)][0];

        // Use v-chord lengths from Pass 1
        std::vector<Real> cs_v(m - 1);
        for (int j = 0; j < m - 1; ++j)
            cs_v[j] = surf.choret_[idx(i, j)][1];

        // Impose r_uv at boundaries (computed in Pass 3)
        ImposedTangents<Dim, Real> imp_tang{
            surf.tanret_[idx(i, 0)][2],
            surf.tanret_[idx(i, m-1)][2]
        };

        auto taus = detail::buildTangents<Dim,Real>(
            std::span{q_work.data(), size_t(m)},
            std::span<const Real>{cs_v},
            EndCondition::Imposed, imp_tang);

        for (int j = 0; j < m; ++j)
            surf.tanret_[idx(i, j)][2] = taus[j];
    }

    return surf;
}

// ============================================================================
// SplineSurface::evalPatch
// Port of Fortran evapa2 (build coefficients) + gpsur/evsurg (evaluate).
// ============================================================================
template<int Dim, RealType Real>
typename SplineSurface<Dim, Real>::Point
SplineSurface<Dim, Real>::evalPatch(int pi, int pj, Real u, Real v) const
{
    const int n = n_;
    auto idx = [n](int i, int j) { return j * n + i; };

    // Flat indices for the four patch corners
    int i00 = idx(pi,   pj);
    int i10 = idx(pi+1, pj);
    int i01 = idx(pi,   pj+1);
    int i11 = idx(pi+1, pj+1);

    // Chord lengths used to scale the tangent/twist vectors in Q
    // (matching Fortran evapa2 variable naming: a,b,c,d)
    Real a = choret_[i00][0];   // u-chord at (pi,   pj)
    Real b = choret_[i10][1];   // v-chord at (pi+1, pj)
    Real c = choret_[i01][0];   // u-chord at (pi,   pj+1)
    Real d = choret_[i00][1];   // v-chord at (pi,   pj)

    Real ab = a*b, bc = b*c, cd = c*d, ad = a*d;

    Point pt{};

    for (int id = 0; id < Dim; ++id) {
        // ----------------------------------------------------------------
        // Build 4×4 Hermite data matrix Q for coordinate dimension `id`.
        // Layout:
        //   Q[row][col]  with  row ∈ {r(u=0), r(u=1), ru(u=0)·a, ru(u=1)·a_}
        //                      col ∈ {v=0, v=1, rv(v=0)·d, rv(v=1)·b}  (scaled)
        // ----------------------------------------------------------------
        std::array<std::array<Real,4>,4> Q{};
        Q[0][0] = coor_[i00][id];
        Q[1][0] = coor_[i10][id];
        Q[2][0] = a * tanret_[i00][0][id];
        Q[3][0] = a * tanret_[i10][0][id];

        Q[0][1] = coor_[i01][id];
        Q[1][1] = coor_[i11][id];
        Q[2][1] = c * tanret_[i01][0][id];
        Q[3][1] = c * tanret_[i11][0][id];

        Q[0][2] = d  * tanret_[i00][1][id];
        Q[1][2] = b  * tanret_[i10][1][id];
        Q[2][2] = ad * tanret_[i00][2][id];
        Q[3][2] = ab * tanret_[i10][2][id];

        Q[0][3] = d  * tanret_[i01][1][id];
        Q[1][3] = b  * tanret_[i11][1][id];
        Q[2][3] = cd * tanret_[i01][2][id];
        Q[3][3] = bc * tanret_[i11][2][id];

        // A[i][j]: coefficient of u^i * v^j
        auto A = detail::patchCoeffs<Real>(Q);

        // ----------------------------------------------------------------
        // Evaluate via nested Horner.
        // Define P_i(v) = A[i][0] + A[i][1]*v + A[i][2]*v^2 + A[i][3]*v^3
        // Then r(u,v) = P_0(v) + P_1(v)*u + P_2(v)*u^2 + P_3(v)*u^3
        // ----------------------------------------------------------------

        // dv[i]  = P_i(v)   (Horner in v for each row i)
        std::array<Real,4> dv;
        for (int i = 0; i < 4; ++i)
            dv[i] = ((A[i][3]*v + A[i][2])*v + A[i][1])*v + A[i][0];

        // dv1[i] = P_i'(v) = dP_i/dv
        std::array<Real,4> dv1;
        for (int i = 0; i < 4; ++i)
            dv1[i] = (Real(3)*A[i][3]*v + Real(2)*A[i][2])*v + A[i][1];

        // dv2[i] = P_i''(v)
        std::array<Real,4> dv2;
        for (int i = 0; i < 4; ++i)
            dv2[i] = Real(6)*A[i][3]*v + Real(2)*A[i][2];

        // Position: r = P_0 + P_1*u + P_2*u^2 + P_3*u^3
        pt.pos[id] = ((dv[3]*u + dv[2])*u + dv[1])*u + dv[0];

        // ∂r/∂u = P_1 + 2*P_2*u + 3*P_3*u^2
        pt.ru[id] = (Real(3)*dv[3]*u + Real(2)*dv[2])*u + dv[1];

        // ∂r/∂v = P_0'(v) + P_1'(v)*u + P_2'(v)*u^2 + P_3'(v)*u^3
        pt.rv[id] = ((dv1[3]*u + dv1[2])*u + dv1[1])*u + dv1[0];

        // ∂²r/∂u² = 2*P_2 + 6*P_3*u
        pt.ruu[id] = Real(6)*dv[3]*u + Real(2)*dv[2];

        // ∂²r/∂v²
        pt.rvv[id] = ((dv2[3]*u + dv2[2])*u + dv2[1])*u + dv2[0];

        // ∂²r/∂u∂v = P_1'(v) + 2*P_2'(v)*u + 3*P_3'(v)*u^2
        pt.ruv[id] = (Real(3)*dv1[3]*u + Real(2)*dv1[2])*u + dv1[1];
    }

    return pt;
}

// ============================================================================
// SplineSurface::eval
// ============================================================================
template<int Dim, RealType Real>
typename SplineSurface<Dim, Real>::Point
SplineSurface<Dim, Real>::eval(Real U, Real V) const
{
    int pi = static_cast<int>(U);
    pi = std::clamp(pi, 0, n_ - 2);
    int pj = static_cast<int>(V);
    pj = std::clamp(pj, 0, m_ - 2);

    Real u = std::clamp(U - Real(pi), Real(0), Real(1));
    Real v = std::clamp(V - Real(pj), Real(0), Real(1));

    return evalPatch(pi, pj, u, v);
}

} // namespace splines
