#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numbers>
#include <vector>
#include "splines/spline_curve.hpp"
#include "splines/spline_surface.hpp"
#include "splines/geometry.hpp"

using namespace splines;
using Catch::Matchers::WithinAbs;
namespace num = std::numbers;

// Circle knots helper
static std::vector<Vec2d> circleKnots2D(int N, double R = 1.0) {
    std::vector<Vec2d> pts(N);
    for (int i = 0; i < N; ++i) {
        double t = 2.0*num::pi * i / N;
        pts[i] = Vec2d{R*std::cos(t), R*std::sin(t)};
    }
    return pts;
}

// Helix knots: r(t) = (a*cos(t), a*sin(t), b*t)
static std::vector<Vec3d> helixKnots(int N, double a = 1.0, double b = 0.5,
                                      double t_max = 2.0*num::pi) {
    std::vector<Vec3d> pts(N);
    for (int i = 0; i < N; ++i) {
        double t = t_max * i / (N - 1);
        pts[i] = Vec3d{a*std::cos(t), a*std::sin(t), b*t};
    }
    return pts;
}

// ============================================================================
// Curve curvature
// ============================================================================

TEST_CASE("Curvature of 2D circle spline ~ 1/R", "[geometry]") {
    double R = 2.0;
    auto pts = circleKnots2D(16, R);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    // Sample curvature at several arcs mid-points
    int pass = 0;
    for (int arc = 0; arc < curve.arcCount(); ++arc) {
        auto pt = curve.eval(double(arc) + 0.5);
        double kappa = curvature(pt.d1, pt.d2);
        if (std::abs(std::abs(kappa) - 1.0/R) < 0.15/R) ++pass;
    }
    // At least 70% of arcs should have curvature near 1/R
    REQUIRE(pass >= curve.arcCount() * 7 / 10);
}

TEST_CASE("Curvature of 3D circle spline ~ 1/R", "[geometry]") {
    double R = 1.5;
    int N = 16;
    std::vector<Vec3d> pts(N);
    for (int i = 0; i < N; ++i) {
        double t = 2.0*num::pi * i / N;
        pts[i] = Vec3d{R*std::cos(t), R*std::sin(t), 0.0};
    }
    auto curve = SplineCurve<3>::interpolate(pts, EndCondition::NotAKnot);

    int pass = 0;
    for (int arc = 0; arc < curve.arcCount(); ++arc) {
        auto pt = curve.eval(double(arc) + 0.5);
        double kappa = curvature(pt.d1, pt.d2);
        if (std::abs(kappa - 1.0/R) < 0.15/R) ++pass;
    }
    REQUIRE(pass >= curve.arcCount() * 7 / 10);
}

TEST_CASE("Curvature of straight line is zero", "[geometry]") {
    std::vector<Vec2d> pts;
    for (int i = 0; i <= 5; ++i) pts.push_back(Vec2d{double(i), 0.0});
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::Natural);

    for (int arc = 0; arc < curve.arcCount(); ++arc) {
        auto pt = curve.eval(double(arc) + 0.5);
        double kappa = curvature(pt.d1, pt.d2);
        REQUIRE_THAT(std::abs(kappa), WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("Torsion of helix ~ b/(a^2+b^2)", "[geometry]") {
    // Theoretical torsion: tau = b/(a^2+b^2)
    double a = 1.0, b = 0.5;
    double expected_tau = b / (a*a + b*b);

    auto pts = helixKnots(40, a, b);
    auto curve = SplineCurve<3>::interpolate(pts, EndCondition::Bessel);

    int pass = 0;
    for (int arc = 2; arc < curve.arcCount() - 2; ++arc) {
        auto pt = curve.eval(double(arc) + 0.5);
        double tau = torsion(pt.d1, pt.d2, pt.d3);
        if (std::abs(tau - expected_tau) < 0.1 * expected_tau) ++pass;
    }
    REQUIRE(pass >= (curve.arcCount() - 4) / 2);
}

// ============================================================================
// Surface geometry
// ============================================================================

static std::vector<Vec3d> sphereLattice(int n, int m, double R = 1.0) {
    std::vector<Vec3d> pts(n * m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) {
            double u = num::pi * j / (m - 1);
            double v = 2.0*num::pi * i / (n - 1);
            pts[j*n+i] = Vec3d{
                R*std::sin(u)*std::cos(v),
                R*std::sin(u)*std::sin(v),
                R*std::cos(u)
            };
        }
    return pts;
}

TEST_CASE("Surface normal on flat plane is (0,0,1)", "[geometry]") {
    int n = 5, m = 5;
    std::vector<Vec3d> lattice(n * m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            lattice[j*n+i] = Vec3d{double(i), double(j), 0.0};

    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    for (double U : {1.0, 2.0, 3.0})
        for (double V : {1.0, 2.0, 3.0}) {
            auto pt = surf.eval(U, V);
            auto nrm = surfaceNormal(pt.ru, pt.rv);
            REQUIRE_THAT(std::abs(nrm[2]), WithinAbs(1.0, 1e-10));
        }
}

TEST_CASE("Surface normal on sphere points radially outward", "[geometry]") {
    int n = 9, m = 7;
    double R = 2.0;
    auto lattice = sphereLattice(n, m, R);
    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    // Sample at lattice knots (exact interpolation points)
    int pass = 0;
    for (int j = 1; j < m - 1; ++j) {
        for (int i = 1; i < n - 1; ++i) {
            auto pt  = surf.eval(double(i), double(j));
            auto nrm = pt.normal();
            // Position unit vector
            auto rhat = normalized(pt.pos);
            // Normal should align with radial direction
            double align = std::abs(dot(nrm, rhat));
            if (align > 0.9) ++pass;
        }
    }
    int total = (m-2)*(n-2);
    REQUIRE(pass >= total * 7 / 10);
}

TEST_CASE("Principal curvatures of sphere ~ 1/R", "[geometry]") {
    int n = 9, m = 7;
    double R = 2.0;
    auto lattice = sphereLattice(n, m, R);
    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    // Sample at a mid-interior point
    auto pt = surf.eval(double(n/2), double(m/2));
    auto [k1, k2] = principalCurvatures(pt.ru, pt.rv, pt.ruu, pt.ruv, pt.rvv);

    // Both curvatures should be near 1/R (sphere has k1=k2=1/R)
    REQUIRE_THAT(std::abs(k1), WithinAbs(1.0/R, 0.5/R));
    REQUIRE_THAT(std::abs(k2), WithinAbs(1.0/R, 0.5/R));
}
