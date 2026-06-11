/**
 * @file geometry.hpp
 * @brief Differential-geometry quantities for spline curves and surfaces.
 *
 * All functions operate on derivative vectors returned by
 * `SplineCurve::eval` or `SplineSurface::eval` and have no dependency
 * on the spline objects themselves — they are pure mathematical utilities.
 *
 * @par Curve quantities
 * | Function          | Formula |
 * |-------------------|---------|
 * | `curvature` (2-D) | Signed κ = (d1[0]·d2[1] − d1[1]·d2[0]) / |d1|³ |
 * | `curvature` (3-D) | Unsigned κ = |d1 × d2| / |d1|³ |
 * | `torsion`   (3-D) | τ = (d1 × d2)·d3 / |d1 × d2|² |
 *
 * @par Surface quantities (3-D only)
 * | Function                | Formula / method |
 * |-------------------------|------------------|
 * | `surfaceNormal`         | n = (rᵤ × r_v) / |rᵤ × r_v| |
 * | `principalCurvatures`   | k₁, k₂ via 1st and 2nd fundamental forms |
 *
 * @note Port of Fortran routines: `getk`, `getk2`, `curvps`, `evsrn`.
 */
#pragma once
#include <cmath>
#include "vec.hpp"
#include "concepts.hpp"

namespace splines {

// ============================================================================
// Curve geometry
// ============================================================================

/**
 * @brief Signed curvature of a 2-D parametric curve.
 *
 * \f[
 *   \kappa = \frac{d_1^{[0]} d_2^{[1]} - d_1^{[1]} d_2^{[0]}}{|\mathbf{d}_1|^3}
 * \f]
 *
 * The sign encodes the turning direction (positive = counter-clockwise).
 *
 * @tparam Real  Floating-point scalar type.
 * @param d1  First derivative dr/du (the velocity vector).
 * @param d2  Second derivative d²r/du².
 * @return    Signed curvature κ, or 0 if |d1| < 1e-20 (degenerate point).
 */
template<RealType Real>
[[nodiscard]] Real curvature(const Vec<2, Real>& d1, const Vec<2, Real>& d2) {
    Real cross2 = d1[0]*d2[1] - d1[1]*d2[0];   // z-component of d1 × d2
    Real speed3 = std::pow(norm(d1), Real(3));
    if (speed3 < Real(1e-20)) return Real(0);
    return cross2 / speed3;
}

/**
 * @brief Unsigned curvature of a 3-D parametric curve.
 *
 * \f[
 *   \kappa = \frac{|\mathbf{d}_1 \times \mathbf{d}_2|}{|\mathbf{d}_1|^3}
 * \f]
 *
 * @tparam Real  Floating-point scalar type.
 * @param d1  First derivative dr/du.
 * @param d2  Second derivative d²r/du².
 * @return    Curvature κ ≥ 0, or 0 if |d1| < 1e-20.
 */
template<RealType Real>
[[nodiscard]] Real curvature(const Vec<3, Real>& d1, const Vec<3, Real>& d2) {
    Real cross_norm = norm(cross(d1, d2));
    Real speed3     = std::pow(norm(d1), Real(3));
    if (speed3 < Real(1e-20)) return Real(0);
    return cross_norm / speed3;
}

/**
 * @brief Torsion of a 3-D parametric curve.
 *
 * \f[
 *   \tau = \frac{(\mathbf{d}_1 \times \mathbf{d}_2) \cdot \mathbf{d}_3}
 *               {|\mathbf{d}_1 \times \mathbf{d}_2|^2}
 * \f]
 *
 * For a helix r(t) = (a cos t, a sin t, b t), torsion = b / (a² + b²).
 *
 * @tparam Real  Floating-point scalar type.
 * @param d1  First derivative dr/du.
 * @param d2  Second derivative d²r/du².
 * @param d3  Third derivative d³r/du³.
 * @return    Torsion τ, or 0 if the binormal is degenerate (|d1 × d2| < 1e-20).
 */
template<RealType Real>
[[nodiscard]] Real torsion(const Vec<3, Real>& d1,
                           const Vec<3, Real>& d2,
                           const Vec<3, Real>& d3) {
    Vec<3, Real> c12 = cross(d1, d2);
    Real den = norm_sq(c12);
    if (den < Real(1e-20)) return Real(0);
    return dot(c12, d3) / den;
}

// ============================================================================
// Surface geometry (3-D only)
// ============================================================================

/**
 * @brief Unit outward normal to a 3-D parametric surface.
 *
 * \f[
 *   \hat{\mathbf{n}} = \frac{\mathbf{r}_u \times \mathbf{r}_v}
 *                           {|\mathbf{r}_u \times \mathbf{r}_v|}
 * \f]
 *
 * @tparam Real  Floating-point scalar type.
 * @param ru  Partial derivative ∂r/∂u.
 * @param rv  Partial derivative ∂r/∂v.
 * @return    Unit normal vector, or the zero vector if the cross product
 *            is degenerate (|rᵤ × r_v| < 1e-20).
 */
template<RealType Real>
[[nodiscard]] Vec<3, Real> surfaceNormal(const Vec<3, Real>& ru,
                                         const Vec<3, Real>& rv) {
    Vec<3, Real> n = cross(ru, rv);
    Real len = norm(n);
    if (len < Real(1e-20)) return Vec<3, Real>{};
    return n / len;
}

/**
 * @brief Result type for `principalCurvatures`.
 *
 * k1 ≥ k2 by convention (larger principal curvature first).
 */
struct PrincipalCurvatures {
    double k1; ///< Larger principal curvature.
    double k2; ///< Smaller principal curvature.
};

/**
 * @brief Principal curvatures of a 3-D parametric surface at a point.
 *
 * Computes k₁ and k₂ via the first and second fundamental forms of
 * differential geometry:
 *
 * - First form: E = rᵤ·rᵤ,  F = rᵤ·r_v,  G = r_v·r_v
 * - Second form: L = n̂·rᵤᵤ,  M = n̂·rᵤ_v,  N = n̂·r_vv
 * - Mean curvature: H = (EN − 2FM + GL) / (2(EG − F²))
 * - Gaussian curvature: K = (LN − M²) / (EG − F²)
 * - Principal curvatures: k₁,₂ = H ± √(H² − K)
 *
 * @tparam Real  Floating-point scalar type.
 * @param ru   ∂r/∂u
 * @param rv   ∂r/∂v
 * @param ruu  ∂²r/∂u²
 * @param ruv  ∂²r/∂u∂v
 * @param rvv  ∂²r/∂v²
 * @return     `{k1, k2}` with k1 ≥ k2, or `{0, 0}` if the surface is
 *             degenerate (EG − F² < 1e-20).
 */
template<RealType Real>
[[nodiscard]] PrincipalCurvatures
principalCurvatures(const Vec<3, Real>& ru,
                    const Vec<3, Real>& rv,
                    const Vec<3, Real>& ruu,
                    const Vec<3, Real>& ruv,
                    const Vec<3, Real>& rvv) {
    // First fundamental form coefficients
    Real E = dot(ru,  ru);
    Real F = dot(ru,  rv);
    Real G = dot(rv,  rv);
    Real det1 = E*G - F*F;
    if (std::abs(det1) < Real(1e-20))
        return {0.0, 0.0};

    // Unit normal
    Vec<3, Real> n = surfaceNormal(ru, rv);

    // Second fundamental form coefficients
    Real L = dot(n, ruu);
    Real M = dot(n, ruv);
    Real N = dot(n, rvv);

    // Mean and Gaussian curvature
    Real H = (E*N - Real(2)*F*M + G*L) / (Real(2) * det1);
    Real K = (L*N - M*M) / det1;

    Real disc = std::max(H*H - K, Real(0));
    Real sq   = std::sqrt(disc);
    return {H + sq, H - sq};
}

} // namespace splines
