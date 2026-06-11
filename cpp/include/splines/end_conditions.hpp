/**
 * @file end_conditions.hpp
 * @brief Boundary condition types for cubic spline interpolation.
 *
 * Defines the `EndCondition` enumeration, which selects how the first and
 * last tangent vectors of a spline are determined, and the `ImposedTangents`
 * helper struct used when tangent vectors are supplied explicitly.
 */
#pragma once
#include <optional>
#include "concepts.hpp"

namespace splines {

/**
 * @brief Selects the end condition used when constructing a cubic spline.
 *
 * For a curve with @e n knots the tridiagonal tangent system has @e n
 * unknowns.  The interior equations are fixed by C² continuity; the two
 * boundary equations are determined by the choice below.
 *
 * | Value       | Mathematical constraint |
 * |-------------|------------------------|
 * | `Natural`   | Zero second derivative at each endpoint: r''(0) = r''(L) = 0. |
 * | `Bessel`    | Tangent at each endpoint is the tangent of the parabola passing through the first (or last) three knots. |
 * | `NotAKnot`  | Third derivative is continuous at the second knot (start side) and second-to-last knot (end side), effectively forcing arcs 0–1 and (n−3)–(n−2) to be the same cubic. Requires ≥ 3 knots. |
 * | `Quadratic` | The second derivative at each endpoint equals the second derivative at its neighbour, reducing the spline to a quadratic polynomial at the ends. |
 * | `Imposed`   | Tangent vectors are supplied explicitly via `ImposedTangents`. |
 *
 * @note With exactly 3 knots `NotAKnot` degenerates (both boundary conditions
 *       encode the same constraint).  The library falls back to `Natural`
 *       automatically in that case.
 */
enum class EndCondition {
    Natural,    ///< Zero second derivative at each endpoint.
    Bessel,     ///< Tangent from parabola through first/last three knots.
    NotAKnot,   ///< Third-derivative continuity at knots 1 and n-2.
    Quadratic,  ///< Second derivative at endpoint equals that at its neighbour.
    Imposed,    ///< Tangent vectors supplied explicitly via ImposedTangents.
};

// Forward-declare Vec so ImposedTangents can reference it.
template<int Dim, RealType Real> struct Vec;

/**
 * @brief Holds optional user-supplied tangent vectors for `EndCondition::Imposed`.
 *
 * When `EndCondition::Imposed` is selected, the tridiagonal system sets
 * the first row to `tau[0] = start` and the last row to `tau[n-1] = end`.
 * Either endpoint may be left as `std::nullopt`; in that case the
 * corresponding tangent defaults to the zero vector.
 *
 * Tangent vectors are expressed as @em dr/dc (derivative with respect to
 * chord length), consistent with the internal tangent convention.
 *
 * @tparam Dim  Spatial dimension.
 * @tparam Real Floating-point scalar type (default `double`).
 */
template<int Dim, RealType Real = double>
struct ImposedTangents {
    std::optional<Vec<Dim, Real>> start; ///< Tangent at the first knot (dr/dc).
    std::optional<Vec<Dim, Real>> end;   ///< Tangent at the last knot  (dr/dc).
};

} // namespace splines
