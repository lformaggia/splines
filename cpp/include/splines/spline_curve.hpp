/**
 * @file spline_curve.hpp
 * @brief Parametric cubic spline curve with chord-length parametrisation.
 *
 * `SplineCurve<Dim, Real>` interpolates an ordered sequence of @em knots in
 * `Dim`-dimensional space using a C² cubic Hermite spline.  The curve uses a
 * **chord-length parametrisation**: global parameter
 * \f$ t = \text{arc\_index} + u \f$ with \f$ u \in [0,1] \f$.
 *
 * @par Tangent convention
 * All stored tangent vectors are @em dr/dc (derivative with respect to chord
 * length @em c, not the local parameter @em u).  Within arc @em i the Hermite
 * formula scales them as `dr/du = chords[i] * tau[i]`.
 *
 * @par End conditions
 * Five boundary conditions are supported via `EndCondition` (see
 * `end_conditions.hpp`).  The default is `NotAKnot`, which produces the
 * "natural" polynomial extension favoured by most numerical texts.
 *
 * @par Fortran origins
 * Port of Fortran routines `trcoe` (tangent system assembly),
 * `evtan` / `trif` + `tris2` (solve), and `getp2` (evaluation).
 *
 * @note The `detail::buildTangents` function also serves `SplineSurface`
 *       during the four-pass surface construction.
 */
#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <optional>
#include <span>
#include <vector>
#include "concepts.hpp"
#include "end_conditions.hpp"
#include "tridiagonal.hpp"
#include "vec.hpp"

namespace splines {

/**
 * @brief Parametric cubic spline curve interpolating an ordered set of knots.
 *
 * @tparam Dim   Spatial dimension (e.g. 2 for planar, 3 for space curves).
 * @tparam Real  Floating-point scalar type (default `double`).
 *
 * @par Parameter space
 * Global parameter @f$ t \in [0,\, \text{arcCount}()] @f$.  Integer values of
 * @f$ t @f$ correspond exactly to knots; each interval @f$ [i, i+1] @f$ is one
 * cubic arc.
 *
 * @par Example
 * @code
 * #include "splines/spline_curve.hpp"
 * using namespace splines;
 *
 * std::vector<Vec3d> pts = {{0,0,0},{1,2,0},{3,1,0},{4,3,0}};
 * auto curve = SplineCurve<3>::interpolate(pts, EndCondition::NotAKnot);
 *
 * // Sample 100 points along the curve
 * for (int k = 0; k <= 100; ++k) {
 *     double t = 3.0 * k / 100.0;
 *     Vec3d pos = curve.position(t);
 * }
 * double L = curve.totalLength();
 * @endcode
 */
template<int Dim, RealType Real = double>
class SplineCurve {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /**
     * @brief Build a cubic spline interpolating @p knots.
     *
     * Computes chord lengths, assembles and solves the tridiagonal tangent
     * system with the requested end condition, and returns the fully
     * initialised curve.
     *
     * @param knots    Ordered sequence of interpolation points.  At least 2
     *                 knots are required; `Bessel` and `NotAKnot` require ≥ 3.
     * @param ec       End condition (default `EndCondition::NotAKnot`).
     * @param imposed  Tangent vectors at the endpoints, used only when
     *                 `ec == EndCondition::Imposed`.
     * @return         Fully constructed `SplineCurve`.
     *
     * @pre  `knots.size() >= 2`.
     */
    static SplineCurve interpolate(
        std::span<const Vec<Dim, Real>> knots,
        EndCondition                    ec      = EndCondition::NotAKnot,
        ImposedTangents<Dim, Real>      imposed = {});

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// @brief Number of knots (control points).
    [[nodiscard]] int knotCount() const noexcept { return static_cast<int>(knots_.size()); }

    /// @brief Number of cubic arcs (= knotCount() - 1).
    [[nodiscard]] int arcCount()  const noexcept { return knotCount() - 1; }

    /// @brief Read-only view of the knot positions.
    [[nodiscard]] std::span<const Vec<Dim, Real>> knots()    const noexcept { return knots_; }

    /**
     * @brief Read-only view of the tangent vectors at each knot.
     *
     * Tangents are stored as @em dr/dc (with respect to chord length).
     * Multiply by `chords()[i]` to obtain `dr/du` at knot @em i.
     */
    [[nodiscard]] std::span<const Vec<Dim, Real>> tangents() const noexcept { return taus_; }

    /// @brief Read-only view of the chord lengths of each arc.
    [[nodiscard]] std::span<const Real>            chords()  const noexcept { return chords_; }

    // -----------------------------------------------------------------------
    // Evaluation
    // -----------------------------------------------------------------------

    /**
     * @brief Position and derivatives at a parameter value.
     *
     * All derivative fields (`d1`, `d2`, `d3`) are with respect to the
     * local arc parameter @em u (not chord length).
     */
    struct Point {
        Vec<Dim, Real> pos;   ///< r(u)       — position.
        Vec<Dim, Real> d1;    ///< dr/du      — first derivative (velocity).
        Vec<Dim, Real> d2;    ///< d²r/du²    — second derivative (acceleration).
        Vec<Dim, Real> d3;    ///< d³r/du³    — third derivative (jerk).
        Real           speed; ///< |dr/du|    — speed (Euclidean norm of d1).
    };

    /**
     * @brief Evaluate the curve at global parameter @p t.
     *
     * @p t is clamped to [0, arcCount()] before evaluation.
     *
     * @param t  Global parameter value.
     * @return   `Point` containing position and first three derivatives.
     */
    [[nodiscard]] Point          eval(Real t)     const;

    /**
     * @brief Return only the position at parameter @p t.
     * @param t  Global parameter value.
     */
    [[nodiscard]] Vec<Dim, Real> position(Real t) const { return eval(t).pos; }

    /**
     * @brief Return the unit tangent direction at parameter @p t.
     * @param t  Global parameter value.
     * @return   Normalised first derivative dr/du.
     * @warning  Undefined if the curve has a cusp (speed = 0) at @p t.
     */
    [[nodiscard]] Vec<Dim, Real> tangent(Real t)  const { return normalized(eval(t).d1); }

    /**
     * @brief Arc length of a single arc via adaptive Romberg quadrature.
     *
     * Integrates the speed polynomial |dr/du| over arc [@p arc, @p arc+1].
     * The integrand is a degree-4 polynomial in u, so Romberg converges
     * rapidly (typically in fewer than 10 refinement levels).
     *
     * @param arc  Zero-based arc index in [0, arcCount()-1].
     * @param eps  Relative tolerance for the Romberg scheme (default 1e-6).
     * @return     Arc length of the requested segment.
     *
     * @pre  0 <= arc < arcCount()
     */
    [[nodiscard]] Real arcLength(int arc, Real eps = Real(1e-6)) const;

    /**
     * @brief Total arc length of the entire curve.
     *
     * Sums `arcLength()` over all arcs.
     *
     * @param eps  Tolerance forwarded to each arc integration (default 1e-6).
     * @return     Total arc length from the first to the last knot.
     */
    [[nodiscard]] Real totalLength(Real eps = Real(1e-6)) const;

private:
    std::vector<Vec<Dim, Real>> knots_;   ///< Knot positions.
    std::vector<Vec<Dim, Real>> taus_;    ///< Tangents dr/dc at each knot.
    std::vector<Real>           chords_;  ///< Chord length of each arc.
};

// ============================================================================
// Implementation
// ============================================================================

namespace detail {

/**
 * @brief Compute Euclidean chord lengths between consecutive knots.
 *
 * A minimum value of 1e-10 is enforced to prevent division-by-zero when
 * two knots coincide.
 *
 * @tparam Dim   Spatial dimension.
 * @tparam Real  Scalar type.
 * @param q  Ordered knot sequence (length n).
 * @return   Vector of n-1 chord lengths.
 */
template<int Dim, RealType Real>
std::vector<Real> computeChords(std::span<const Vec<Dim, Real>> q) {
    const int na = static_cast<int>(q.size()) - 1;
    std::vector<Real> cs(na);
    constexpr Real eps = Real(1e-10);
    for (int i = 0; i < na; ++i)
        cs[i] = std::max(norm(q[i + 1] - q[i]), eps);
    return cs;
}

/**
 * @brief Assemble and solve the tridiagonal tangent system.
 *
 * Constructs the tridiagonal system for the C² chord-parametrised cubic
 * spline and solves it via the Thomas algorithm.  Returns the tangent
 * vectors τᵢ = dr/dc at each knot.
 *
 * @par Interior equations (knot i, 1 ≤ i ≤ n−2)
 * @code
 *   cs[i]*tau[i-1] + 2*(cs[i]+cs[i-1])*tau[i] + cs[i-1]*tau[i+1]
 *     = 3/(cs[i]*cs[i-1]) * (cs[i-1]^2*Δright + cs[i]^2*Δleft)
 * @endcode
 *
 * @par Boundary equations
 * Determined by @p ec; see `end_conditions.hpp`.
 *
 * @note For `NotAKnot` the Fortran formula creates a zero modified diagonal
 *       (exploiting pre-IEEE `0/0 = 0`).  This implementation uses a
 *       mathematically equivalent but non-degenerate reformulation derived
 *       by eliminating τ[1] (start) or τ[n-2] (end) via the adjacent C²
 *       equation.
 *
 * @note For n = 3 with `NotAKnot`, both boundary conditions encode the same
 *       constraint (the single interior knot is always not-a-knot).  The
 *       implementation falls back to `Natural` in that case.
 *
 * @tparam Dim   Spatial dimension.
 * @tparam Real  Scalar type.
 * @param q    Knot positions.
 * @param cs   Chord lengths (length n-1).
 * @param ec   End condition selector.
 * @param imp  Imposed tangent vectors (only used when ec == Imposed).
 * @return     Vector of n tangent vectors τᵢ = dr/dc.
 */
template<int Dim, RealType Real>
std::vector<Vec<Dim, Real>>
buildTangents(std::span<const Vec<Dim, Real>> q,
              std::span<const Real>           cs,
              EndCondition                    ec,
              const ImposedTangents<Dim, Real>& imp) {
    const int n = static_cast<int>(q.size());
    assert(n >= 2);

    constexpr Real eps = Real(1e-10);

    // Not-a-knot with n==3 is degenerate when applied at both ends
    // (both conditions encode the same constraint). Fall back to Natural.
    EndCondition ec_eff = (n == 3 && ec == EndCondition::NotAKnot)
                          ? EndCondition::Natural : ec;

    std::vector<Real>           a(n, Real(0)), b(n, Real(0)), c(n, Real(0));
    std::vector<Vec<Dim, Real>> rhs(n);

    // ----------------------------------------------------------------
    // Interior rows  (knots 1 .. n-2, 0-indexed)
    // a[i] = cs[i]   (right-arc chord)
    // b[i] = 2*(cs[i]+cs[i-1])
    // c[i] = cs[i-1] (left-arc chord)
    // rhs  = 3/(cs[i]*cs[i-1]) * (cs[i-1]^2*Δ_right + cs[i]^2*Δ_left)
    // ----------------------------------------------------------------
    for (int i = 1; i < n - 1; ++i) {
        Real ci  = std::max(cs[i],     eps);   // right arc
        Real ci1 = std::max(cs[i - 1], eps);   // left  arc
        a[i] = ci;
        b[i] = Real(2) * (ci + ci1);
        c[i] = ci1;
        rhs[i] = (Real(3) / (ci * ci1)) *
                 (ci1 * ci1 * (q[i + 1] - q[i]) +
                  ci  * ci  * (q[i]     - q[i - 1]));
    }

    // ----------------------------------------------------------------
    // Boundary conditions – start
    // ----------------------------------------------------------------
    auto applyStart = [&]() {
        switch (ec_eff) {
        case EndCondition::Natural:
            a[0] = Real(0);
            b[0] = std::max(Real(2) * cs[0], eps);
            c[0] = cs[0];
            rhs[0] = Real(3) * (q[1] - q[0]);
            break;

        case EndCondition::Imposed:
            a[0] = Real(0); b[0] = Real(1); c[0] = Real(0);
            rhs[0] = imp.start.value_or(Vec<Dim, Real>{});
            break;

        case EndCondition::Bessel: {
            assert(n >= 3 && "Bessel requires >= 3 knots");
            Real beta = cs[0] + cs[1];          // b(2)/2 in Fortran
            a[0] = Real(0); b[0] = Real(1); c[0] = Real(0);
            rhs[0] = -q[0] * ((Real(2)*cs[0] + cs[1]) / (beta * cs[0]))
                   +  q[1] * (beta / (cs[0] * cs[1]))
                   -  q[2] * (cs[0] / (cs[1] * beta));
            break;
        }

        case EndCondition::Quadratic:
            a[0] = Real(0); b[0] = Real(1); c[0] = Real(1);
            rhs[0] = Real(2) * (q[1] - q[0]) / cs[0];
            break;

        case EndCondition::NotAKnot: {
            assert(n >= 3 && "NotAKnot requires >= 3 knots");
            // Fortran trcoe uses a degenerate formula relying on pre-IEEE 0/0=0.
            // Correct derivation: combine the d³r/dc³ continuity condition at
            // knot 1 with the C² equation there to eliminate tau[2], yielding a
            // proper non-degenerate first row.
            Real h0 = cs[0], h1 = cs[1];
            a[0] = Real(0);
            b[0] = h1;
            c[0] = h0 + h1;
            // numerator = Fortran rhs_start; divide by (h0+h1) to normalise
            Vec<Dim, Real> num0 = (q[2]-q[1]) * (h0*h0/h1)
                                + (q[1]-q[0]) * (h1/h0) * (Real(3)*h0 + Real(2)*h1);
            rhs[0] = num0 * (Real(1) / (h0 + h1));
            break;
        }
        }
    };

    // ----------------------------------------------------------------
    // Boundary conditions – end
    // ----------------------------------------------------------------
    auto applyEnd = [&]() {
        switch (ec_eff) {
        case EndCondition::Natural:
            a[n-1] = cs[n-2];
            b[n-1] = std::max(Real(2) * cs[n-2], eps);
            c[n-1] = Real(0);
            rhs[n-1] = Real(3) * (q[n-1] - q[n-2]);
            break;

        case EndCondition::Imposed:
            a[n-1] = Real(0); b[n-1] = Real(1); c[n-1] = Real(0);
            rhs[n-1] = imp.end.value_or(Vec<Dim, Real>{});
            break;

        case EndCondition::Bessel: {
            assert(n >= 3 && "Bessel requires >= 3 knots");
            Real beta = cs[n-2] + cs[n-3];     // b(n-1)/2 in Fortran
            a[n-1] = Real(0); b[n-1] = Real(1); c[n-1] = Real(0);
            rhs[n-1] =  q[n-1] * ((Real(2)*cs[n-2] + cs[n-3]) / (beta * cs[n-2]))
                      - q[n-2] * (beta / (cs[n-2] * cs[n-3]))
                      + q[n-3] * (cs[n-2] / (cs[n-3] * beta));
            break;
        }

        case EndCondition::Quadratic:
            a[n-1] = Real(1); b[n-1] = Real(1); c[n-1] = Real(0);
            rhs[n-1] = Real(2) * (q[n-1] - q[n-2]) / cs[n-2];
            break;

        case EndCondition::NotAKnot: {
            assert(n >= 3 && "NotAKnot requires >= 3 knots");
            // Symmetric to the start fix: eliminate tau[n-3] using the C²
            // condition at row n-2 to obtain a non-degenerate last row.
            Real hm2 = cs[n-2], hm3 = cs[n-3];
            a[n-1] = hm2 + hm3;
            b[n-1] = hm3;
            c[n-1] = Real(0);
            Vec<Dim, Real> numN = (q[n-2]-q[n-3]) * (hm2*hm2/hm3)
                                + (q[n-1]-q[n-2]) * (hm3/hm2) * (Real(3)*hm2 + Real(2)*hm3);
            rhs[n-1] = numN * (Real(1) / (hm2 + hm3));
            break;
        }
        }
    };

    applyStart();
    applyEnd();

    // ----------------------------------------------------------------
    // Factor + solve
    // ----------------------------------------------------------------
    detail::tridiag_factor<Real>(a, b, c);
    detail::tridiag_solve<Real>(
        std::span<const Real>{a},
        std::span<const Real>{b},
        std::span<const Real>{c},
        std::span<Vec<Dim, Real>>{rhs});

    return rhs;
}

} // namespace detail

// ----------------------------------------------------------------------------
// SplineCurve::interpolate
// ----------------------------------------------------------------------------
template<int Dim, RealType Real>
SplineCurve<Dim, Real>
SplineCurve<Dim, Real>::interpolate(
    std::span<const Vec<Dim, Real>> knots,
    EndCondition                    ec,
    ImposedTangents<Dim, Real>      imposed)
{
    assert(static_cast<int>(knots.size()) >= 2 && "SplineCurve requires at least 2 knots");

    SplineCurve curve;
    curve.knots_.assign(knots.begin(), knots.end());
    curve.chords_ = detail::computeChords<Dim, Real>(knots);
    curve.taus_   = detail::buildTangents<Dim, Real>(
                        knots, curve.chords_, ec, imposed);
    return curve;
}

// ----------------------------------------------------------------------------
// SplineCurve::eval
//
// Port of Fortran getp2.  Coefficients (per dimension):
//   a1 = 2*(q[is] - q[is+1]) + cs*(tau[is]+tau[is+1])   (u^3 coeff of x-q[is])
//   a2 = -3*(q[is]-q[is+1]) - cs*(2*tau[is]+tau[is+1])  (u^2 coeff)
//   a3 = cs * tau[is]                                     (u^1 coeff)
//   x(u) = ((a1*u + a2)*u + a3)*u + q[is]
// ----------------------------------------------------------------------------
template<int Dim, RealType Real>
typename SplineCurve<Dim, Real>::Point
SplineCurve<Dim, Real>::eval(Real t) const
{
    const int n = knotCount();
    // Clamp to valid range
    int arc = static_cast<int>(t);
    arc = std::clamp(arc, 0, n - 2);
    Real u = t - Real(arc);
    u = std::clamp(u, Real(0), Real(1));

    const Real cs = chords_[arc];
    const Vec<Dim, Real>& q0  = knots_[arc];
    const Vec<Dim, Real>& q1  = knots_[arc + 1];
    const Vec<Dim, Real>& t0  = taus_[arc];
    const Vec<Dim, Real>& t1  = taus_[arc + 1];

    // r12 = q[is] - q[is+1]  (Fortran convention)
    Vec<Dim, Real> r12 = q0 - q1;
    Vec<Dim, Real> a1  = Real(2)*r12  + cs*(t0 + t1);
    Vec<Dim, Real> a2  = Real(-3)*r12 - cs*(Real(2)*t0 + t1);
    Vec<Dim, Real> a3  = cs * t0;

    Point p;
    p.pos  = ((a1*u + a2)*u + a3)*u + q0;
    p.d1   = (Real(3)*a1*u + Real(2)*a2)*u + a3;
    p.d2   = Real(6)*a1*u + Real(2)*a2;
    p.d3   = Real(6)*a1;
    p.speed = norm(p.d1);
    return p;
}

// ----------------------------------------------------------------------------
// SplineCurve::arcLength
//
// Integrates |dr/du| over arc [arc, arc+1] using adaptive Romberg quadrature.
// The integrand is sqrt(C1 + C2*u + C3*u^2 + C4*u^3 + C5*u^4) on [0,1].
// Coefficients are obtained by squaring dr/du (a degree-2 polynomial in u)
// component-wise and summing (port of Fortran coeff.f).
// ----------------------------------------------------------------------------
template<int Dim, RealType Real>
Real SplineCurve<Dim, Real>::arcLength(int arc, Real eps) const
{
    assert(arc >= 0 && arc < arcCount());

    // Polynomial coefficients of |dr/du|^2 = a1 + a2*u + a3*u^2 + a4*u^3 + a5*u^4
    // (Port of Fortran coeff.f)
    const Real cs = chords_[arc];
    const Vec<Dim, Real>& q0 = knots_[arc];
    const Vec<Dim, Real>& q1 = knots_[arc + 1];
    const Vec<Dim, Real>& t0 = taus_[arc];
    const Vec<Dim, Real>& t1 = taus_[arc + 1];

    Real C1{}, C2{}, C3{}, C4{}, C5{};
    for (int id = 0; id < Dim; ++id) {
        Real r12 = q0[id] - q1[id];
        Real p   = cs * t0[id];                                 // = a3 coeff
        Real z   = cs * t1[id];
        Real s   = Real(3) * (Real(2)*r12 + p + z);            // 3*a1
        Real f   = Real(2) * (Real(-3)*r12 - Real(2)*p - z);   // 2*a2  (= 2*a2 from getp2)
        // Note: dr/du = s*u^2 + f*u + p  (checking: at u=0, dr/du = p = cs*t0 ✓)
        C1 += p*p;
        C2 += Real(2)*p*f;
        C3 += f*f + Real(2)*p*s;
        C4 += Real(2)*f*s;
        C5 += s*s;
    }

    // Romberg quadrature of sqrt(C1 + C2*u + C3*u^2 + C4*u^3 + C5*u^4) on [0,1]
    // Based on Numerical Recipes adaptive Romberg (port of qrombs/trapzds/polints).
    auto integrand = [&](Real u) -> Real {
        return std::sqrt(std::max(C1 + u*(C2 + u*(C3 + u*(C4 + u*C5))), Real(0)));
    };

    // Adaptive Romberg (up to 20 refinements, 2-point Richardson extrapolation)
    constexpr int JMAX = 20;
    constexpr int K    = 5;   // use last K points for extrapolation

    std::array<Real, JMAX + 1> s_trap{}, h_arr{};
    h_arr[0] = Real(1);

    // Trapezoid with progressive refinements
    Real trap_s{};
    int  it = 1;
    for (int j = 0; j < JMAX; ++j) {
        if (j == 0) {
            trap_s = Real(0.5) * (integrand(Real(0)) + integrand(Real(1)));
        } else {
            Real del = Real(1) / Real(it);
            Real x   = Real(0.5) * del;
            Real sum{};
            for (int k = 0; k < it; ++k, x += del)
                sum += integrand(x);
            trap_s = Real(0.5) * (trap_s + sum * del);
            it *= 2;
        }
        s_trap[j] = trap_s;
        h_arr[j]  = (j == 0) ? Real(1) : Real(0.25) * h_arr[j - 1];

        // Richardson extrapolation once we have at least K points
        if (j >= K - 1) {
            int start = j - (K - 1);
            // Neville polynomial interpolation to h=0
            std::array<Real, K> c_arr, d_arr;
            for (int i = 0; i < K; ++i)
                c_arr[i] = d_arr[i] = s_trap[start + i];

            Real y = s_trap[start];
            Real dy{};
            int ns = 0;
            Real dif = std::abs(h_arr[start]);
            for (int i = 0; i < K; ++i) {
                Real dift = std::abs(h_arr[start + i]);
                if (dift < dif) { ns = i; dif = dift; }
            }
            y = s_trap[start + ns];
            --ns;
            for (int m = 1; m < K; ++m) {
                for (int i = 0; i < K - m; ++i) {
                    Real ho  = h_arr[start + i];
                    Real hp  = h_arr[start + i + m];
                    Real w   = c_arr[i + 1] - d_arr[i];
                    Real den = ho - hp;
                    if (den == Real(0)) break;
                    den      = w / den;
                    d_arr[i] = hp * den;
                    c_arr[i] = ho * den;
                }
                dy = (2 * (ns + 1) < K - m) ? c_arr[ns + 1] : d_arr[ns--];
                y += dy;
            }
            if (std::abs(dy) <= eps * std::abs(y)) return y;
        }
        s_trap[j + 1] = s_trap[j];
        h_arr[j + 1]  = Real(0.25) * h_arr[j];
    }
    return s_trap[JMAX - 1];  // return best estimate even if not converged
}

template<int Dim, RealType Real>
Real SplineCurve<Dim, Real>::totalLength(Real eps) const {
    Real len{};
    for (int i = 0; i < arcCount(); ++i)
        len += arcLength(i, eps);
    return len;
}

} // namespace splines
