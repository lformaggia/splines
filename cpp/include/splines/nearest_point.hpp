/**
 * @file nearest_point.hpp
 * @brief Newton-iteration nearest-point algorithms for curves and surfaces.
 *
 * Provides `nearestOnCurve` and `nearestOnSurface`, which find the point on
 * a spline object closest to a given query point by minimising the squared
 * Euclidean distance via Newton iteration.
 *
 * @par Curve algorithm (port of Fortran `locpfs`)
 * Minimises @f$ f(t) = \tfrac{1}{2} \| \mathbf{r}(t) - \mathbf{q} \|^2 @f$
 * by iterating @f$ t \leftarrow t - f'(t) / f''(t) @f$.  If no initial
 * guess is supplied, a brute-force scan over the knots provides the seed.
 *
 * @par Surface algorithm (port of Fortran `locsug`)
 * Minimises @f$ f(U,V) = \tfrac{1}{2} \| \mathbf{r}(U,V) - \mathbf{q} \|^2 @f$
 * using a 2-D Newton step with an approximate positive-definite Hessian.
 * Steps are clamped to length 1 and the parameters are projected back into
 * the valid domain after each iteration.
 */
#pragma once
#include <algorithm>
#include <cmath>
#include <optional>
#include "spline_curve.hpp"
#include "spline_surface.hpp"

namespace splines {

// ============================================================================
// Curve nearest point
// ============================================================================

/**
 * @brief Result of a nearest-point search on a curve.
 *
 * @tparam Dim   Spatial dimension.
 * @tparam Real  Floating-point scalar type.
 */
template<int Dim, RealType Real>
struct NearestCurveResult {
    Real           param;     ///< Global curve parameter t at the nearest point.
    Vec<Dim, Real> point;     ///< Coordinates of the nearest point r(t).
    Real           distance;  ///< Euclidean distance from the query to r(t).
    bool           converged; ///< `true` if Newton iteration converged.
};

/**
 * @brief Find the point on @p curve nearest to @p query.
 *
 * Performs Newton iteration on
 * @f$ f'(t) = (\mathbf{r}(t) - \mathbf{q}) \cdot \mathbf{r}'(t) = 0 @f$.
 * The step is
 * @f[
 *   \Delta t = -\frac{(\mathbf{r}-\mathbf{q})\cdot\mathbf{r}'}
 *                    {|\mathbf{r}'|^2 + (\mathbf{r}-\mathbf{q})\cdot\mathbf{r}''}
 * @f]
 * clamped to |Δt| ≤ 1 to avoid overshooting.
 *
 * Convergence criteria (whichever is satisfied first):
 * - Distance to query ≤ 0.1% of the local chord length.
 * - Gradient magnitude < 0.1% of max(chord, distance).
 * - Parameter change < 1e-6 between consecutive iterations.
 *
 * @tparam Dim   Spatial dimension.
 * @tparam Real  Floating-point scalar type.
 *
 * @param curve  The spline curve to project onto.
 * @param query  The query point.
 * @param t0     Optional initial parameter guess.  If negative (default),
 *               a brute-force knot scan provides the seed.
 * @return       `NearestCurveResult` with the nearest parameter, point,
 *               distance, and convergence flag.
 */
template<int Dim, RealType Real>
[[nodiscard]] NearestCurveResult<Dim, Real>
nearestOnCurve(const SplineCurve<Dim, Real>& curve,
               const Vec<Dim, Real>&         query,
               Real                          t0 = Real(-1))
{
    const int n = curve.knotCount();
    const Real t_min = Real(0);
    const Real t_max = Real(n - 1);

    constexpr Real eps1 = Real(0.001);
    constexpr Real eps2 = Real(1e-6);
    constexpr Real eps3 = Real(0.001);

    // Initial guess: sample at knots and pick closest
    Real ug = t0;
    if (ug < t_min || ug > t_max) {
        Real best_dist = std::numeric_limits<Real>::max();
        for (int i = 0; i < n; ++i) {
            Real d = norm(curve.knots()[i] - query);
            if (d < best_dist) { best_dist = d; ug = Real(i); }
        }
    }

    Real ug0 = Real(-1e30);

    for (int iter = 0; iter < 200; ++iter) {
        int arc = static_cast<int>(ug);
        arc     = std::clamp(arc, 0, n - 2);
        Real u  [[maybe_unused]] = ug - Real(arc);

        Real xl = curve.chords()[arc];
        auto pt = curve.eval(ug);

        Vec<Dim, Real> v = pt.pos - query;
        Real h    = dot(v, pt.d1);
        Real dist = norm(v);

        if (dist <= eps1 * xl)
            return {ug, pt.pos, dist, true};

        Real beta = std::abs(h);
        Real xd   = std::max(xl, dist) * eps1;
        if (beta < xd)
            return {ug, pt.pos, dist, false};   // gradient zero → local min

        if (std::abs(ug0 - ug) < eps2)
            return {ug, pt.pos, dist, true};    // converged (no movement)

        ug0 = ug;
        h  /= beta;

        Real xden = dot(pt.d1, pt.d1);
        Real xpp  = dot(v, pt.d2);
        if ((xpp + xden) > eps3 * xl * xl)
            xden = xpp + xden;

        Real da = -beta / xden;
        if (std::abs(da) > Real(1)) da = std::copysign(Real(1), da);

        ug = std::clamp(ug + da * h, t_min, t_max);
    }

    auto pt = curve.eval(ug);
    return {ug, pt.pos, norm(pt.pos - query), false};
}

// ============================================================================
// Surface nearest point
// ============================================================================

/**
 * @brief Result of a nearest-point search on a surface.
 *
 * @tparam Dim   Spatial dimension.
 * @tparam Real  Floating-point scalar type.
 */
template<int Dim, RealType Real>
struct NearestSurfaceResult {
    Real           U;         ///< Global surface parameter U at the nearest point.
    Real           V;         ///< Global surface parameter V at the nearest point.
    Vec<Dim, Real> point;     ///< Coordinates of the nearest point r(U, V).
    Real           distance;  ///< Euclidean distance from the query to r(U, V).
    bool           converged; ///< `true` if Newton iteration converged.
};

/**
 * @brief Find the point on @p surf nearest to @p query.
 *
 * Performs 2-D Newton iteration on
 * @f$ \nabla f = [\mathbf{r}_u \cdot \mathbf{v},\; \mathbf{r}_v \cdot \mathbf{v}] = \mathbf{0} @f$
 * where @f$ \mathbf{v} = \mathbf{r}(U,V) - \mathbf{q} @f$.
 *
 * The Newton step uses the approximate Hessian
 * @f[
 *   H \approx \begin{bmatrix}
 *     \mathbf{r}_u \cdot \mathbf{r}_u + \mathbf{v} \cdot \mathbf{r}_{uu} &
 *     \mathbf{r}_u \cdot \mathbf{r}_v + \mathbf{v} \cdot \mathbf{r}_{uv} \\
 *     \mathbf{r}_u \cdot \mathbf{r}_v + \mathbf{v} \cdot \mathbf{r}_{uv} &
 *     \mathbf{r}_v \cdot \mathbf{r}_v + \mathbf{v} \cdot \mathbf{r}_{vv}
 *   \end{bmatrix}
 * @f]
 * with gradient-descent fallback when det(H) < 1e-30.
 * Each step is clamped to length ≤ 1 in parameter space, and (U, V) are
 * projected back into [0, n-1] × [0, m-1] after each update.
 *
 * @tparam Dim   Spatial dimension.
 * @tparam Real  Floating-point scalar type.
 *
 * @param surf  The spline surface to project onto.
 * @param query The query point.
 * @param U0    Optional initial U guess; negative (default) triggers lattice scan.
 * @param V0    Optional initial V guess; negative (default) triggers lattice scan.
 * @return      `NearestSurfaceResult` with the nearest (U, V), point, distance,
 *              and convergence flag.
 */
template<int Dim, RealType Real>
[[nodiscard]] NearestSurfaceResult<Dim, Real>
nearestOnSurface(const SplineSurface<Dim, Real>& surf,
                 const Vec<Dim, Real>&           query,
                 Real                            U0 = Real(-1),
                 Real                            V0 = Real(-1))
{
    const Real U_min = Real(0), U_max = Real(surf.n() - 1);
    const Real V_min = Real(0), V_max = Real(surf.m() - 1);

    // Brute-force initial guess over the lattice
    Real Ug = U0, Vg = V0;
    if (Ug < U_min || Ug > U_max || Vg < V_min || Vg > V_max) {
        Real best = std::numeric_limits<Real>::max();
        for (int j = 0; j < surf.m(); ++j)
            for (int i = 0; i < surf.n(); ++i) {
                auto pt = surf.eval(Real(i), Real(j));
                Real d  = norm(pt.pos - query);
                if (d < best) { best = d; Ug = Real(i); Vg = Real(j); }
            }
    }

    constexpr Real tol = Real(1e-8);
    constexpr int  MAX_ITER = 200;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        auto pt = surf.eval(Ug, Vg);
        Vec<Dim, Real> v = pt.pos - query;
        Real dist = norm(v);

        if (dist < tol)
            return {Ug, Vg, pt.pos, dist, true};

        // Gradient of 0.5*||r - q||²:  g = [ru·v, rv·v]
        Real gu = dot(pt.ru, v);
        Real gv = dot(pt.rv, v);

        if (std::sqrt(gu*gu + gv*gv) < tol * dist)
            return {Ug, Vg, pt.pos, dist, true};

        // Newton step using approximate Hessian from 1st-order terms
        // H ≈ [[ru·ru, ru·rv], [ru·rv, rv·rv]]  (positive-definite)
        Real H11 = dot(pt.ru, pt.ru) + dot(v, pt.ruu);
        Real H12 = dot(pt.ru, pt.rv) + dot(v, pt.ruv);
        Real H22 = dot(pt.rv, pt.rv) + dot(v, pt.rvv);

        Real detH = H11*H22 - H12*H12;
        Real dU, dV;
        if (std::abs(detH) > Real(1e-30)) {
            dU = -(H22*gu - H12*gv) / detH;
            dV = -(H11*gv - H12*gu) / detH;
        } else {
            // Fall back to gradient step
            dU = -gu;
            dV = -gv;
        }

        // Clamp step size to avoid overshooting
        Real step = std::sqrt(dU*dU + dV*dV);
        if (step > Real(1)) { dU /= step; dV /= step; }

        Real Ug_new = std::clamp(Ug + dU, U_min, U_max);
        Real Vg_new = std::clamp(Vg + dV, V_min, V_max);

        if (std::abs(Ug_new - Ug) < tol && std::abs(Vg_new - Vg) < tol)
            return {Ug, Vg, pt.pos, dist, true};

        Ug = Ug_new;
        Vg = Vg_new;
    }

    auto pt = surf.eval(Ug, Vg);
    return {Ug, Vg, pt.pos, norm(pt.pos - query), false};
}

} // namespace splines
