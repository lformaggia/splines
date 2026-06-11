/**
 * @file splines.hpp
 * @brief Master include for the splines C++20 library.
 *
 * Including this single header pulls in the entire public API:
 *
 * | Header               | Contents |
 * |----------------------|----------|
 * | `concepts.hpp`       | `RealType` concept |
 * | `vec.hpp`            | `Vec<Dim,Real>` spatial vector |
 * | `end_conditions.hpp` | `EndCondition` enum, `ImposedTangents` |
 * | `tridiagonal.hpp`    | Thomas algorithm (`detail` namespace) |
 * | `spline_curve.hpp`   | `SplineCurve<Dim,Real>` |
 * | `arc_length.hpp`     | Free-function arc-length wrappers |
 * | `geometry.hpp`       | Curvature, torsion, surface normal, principal κ |
 * | `spline_surface.hpp` | `SplineSurface<Dim,Real>` |
 * | `nearest_point.hpp`  | `nearestOnCurve`, `nearestOnSurface` |
 *
 * @par Minimal example
 * @code
 * #include "splines/splines.hpp"
 * using namespace splines;
 *
 * // Interpolate a planar curve through four points
 * std::vector<Vec2d> pts = {{0,0},{1,2},{3,1},{4,3}};
 * auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);
 *
 * // Sample the curve at 100 equally-spaced parameter values
 * for (int k = 0; k <= 100; ++k) {
 *     double t = 3.0 * k / 100.0;
 *     Vec2d p = curve.position(t);
 * }
 * @endcode
 */
#pragma once

#include "splines/concepts.hpp"
#include "splines/vec.hpp"
#include "splines/end_conditions.hpp"
#include "splines/tridiagonal.hpp"
#include "splines/spline_curve.hpp"
#include "splines/arc_length.hpp"
#include "splines/geometry.hpp"
#include "splines/spline_surface.hpp"
#include "splines/nearest_point.hpp"
