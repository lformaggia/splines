/**
 * @file arc_length.hpp
 * @brief Free-function wrappers for spline arc-length computation.
 *
 * Provides thin, free-function wrappers around `SplineCurve::arcLength` and
 * `SplineCurve::totalLength`.  The actual computation (Romberg quadrature of
 * the speed polynomial |dr/du|) lives in `spline_curve.hpp`.  Include this
 * header when you prefer the free-function call style.
 *
 * @see SplineCurve::arcLength
 * @see SplineCurve::totalLength
 */
#pragma once
#include "spline_curve.hpp"

namespace splines {

/**
 * @brief Compute the arc length of a single arc of a cubic spline.
 *
 * Integrates the speed |dr/du| over the arc from knot @p arc to knot
 * @p arc+1 using adaptive Romberg quadrature.
 *
 * @tparam Dim  Spatial dimension.
 * @tparam Real Floating-point scalar type.
 *
 * @param curve  The spline curve to measure.
 * @param arc    Zero-based arc index in [0, curve.arcCount()-1].
 * @param eps    Relative convergence tolerance for the Romberg scheme
 *               (default 1e-6).
 * @return       Arc length of the requested segment.
 *
 * @pre  0 <= arc < curve.arcCount()
 */
template<int Dim, RealType Real>
[[nodiscard]] Real arcLength(const SplineCurve<Dim, Real>& curve,
                             int  arc,
                             Real eps = Real(1e-6)) {
    return curve.arcLength(arc, eps);
}

/**
 * @brief Compute the total arc length of a cubic spline.
 *
 * Sums `arcLength()` over all arcs of the curve.
 *
 * @tparam Dim  Spatial dimension.
 * @tparam Real Floating-point scalar type.
 *
 * @param curve  The spline curve to measure.
 * @param eps    Relative convergence tolerance forwarded to each arc
 *               integration (default 1e-6).
 * @return       Total arc length from the first to the last knot.
 */
template<int Dim, RealType Real>
[[nodiscard]] Real totalLength(const SplineCurve<Dim, Real>& curve,
                               Real eps = Real(1e-6)) {
    return curve.totalLength(eps);
}

} // namespace splines
