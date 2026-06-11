#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numbers>
#include <vector>
#include "splines/spline_curve.hpp"
#include "splines/spline_surface.hpp"
#include "splines/nearest_point.hpp"

using namespace splines;
using Catch::Matchers::WithinAbs;
namespace num = std::numbers;

static std::vector<Vec2d> circleKnots2D(int N, double R = 1.0) {
    std::vector<Vec2d> pts(N);
    for (int i = 0; i < N; ++i) {
        double t = 2.0 * num::pi * i / N;
        pts[i] = Vec2d{R * std::cos(t), R * std::sin(t)};
    }
    return pts;
}

static std::vector<Vec3d> sphereLattice(int n, int m, double R = 1.0) {
    std::vector<Vec3d> pts(n * m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) {
            double u = num::pi * j / (m - 1);
            double v = 2.0 * num::pi * i / (n - 1);
            pts[j*n+i] = Vec3d{R*std::sin(u)*std::cos(v),
                               R*std::sin(u)*std::sin(v),
                               R*std::cos(u)};
        }
    return pts;
}

TEST_CASE("Nearest on curve: query on curve gives distance ~ 0", "[nearest_point]") {
    auto pts = circleKnots2D(12);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    // Pick a knot as the query point — distance must be ~0
    Vec2d query = pts[3];
    auto res = nearestOnCurve(curve, query);

    REQUIRE_THAT(res.distance, WithinAbs(0.0, 1e-6));
}

TEST_CASE("Nearest on curve: point outside circle", "[nearest_point]") {
    // Unit circle spline; query = (2, 0) → nearest point ~ (1, 0), distance ~ 1
    auto pts = circleKnots2D(16);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    Vec2d query{2.0, 0.0};
    auto res = nearestOnCurve(curve, query);

    REQUIRE_THAT(res.distance, WithinAbs(1.0, 0.05));
    REQUIRE_THAT(res.point[0], WithinAbs(1.0, 0.05));
    REQUIRE_THAT(std::abs(res.point[1]), WithinAbs(0.0, 0.05));
}

TEST_CASE("Nearest on curve: 3D helix, query on curve", "[nearest_point]") {
    // Build a 3D helix and query with a point that is an interpolation knot
    int N = 20;
    double a = 1.0, b = 0.3;
    std::vector<Vec3d> pts(N);
    for (int i = 0; i < N; ++i) {
        double t = 2.0 * num::pi * i / (N - 1);
        pts[i] = Vec3d{a*std::cos(t), a*std::sin(t), b*t};
    }
    auto curve = SplineCurve<3>::interpolate(pts, EndCondition::Bessel);

    Vec3d query = pts[7];
    auto res = nearestOnCurve(curve, query);

    REQUIRE_THAT(res.distance, WithinAbs(0.0, 1e-5));
}

TEST_CASE("Nearest on curve: straight line, nearest to off-axis point", "[nearest_point]") {
    // Line along x-axis; query (2, 3) → nearest = (2, 0), distance = 3
    std::vector<Vec2d> pts;
    for (int i = 0; i <= 5; ++i) pts.push_back(Vec2d{double(i), 0.0});
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::Natural);

    Vec2d query{2.0, 3.0};
    auto res = nearestOnCurve(curve, query);

    REQUIRE_THAT(res.distance, WithinAbs(3.0, 1e-4));
    REQUIRE_THAT(res.point[0], WithinAbs(2.0, 1e-4));
    REQUIRE_THAT(res.point[1], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Nearest on surface: query on surface gives distance ~ 0", "[nearest_point]") {
    int n = 7, m = 5;
    double R = 1.5;
    auto lattice = sphereLattice(n, m, R);
    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    // Query at a lattice knot — should return ~0 distance
    Vec3d query = lattice[2*n + 3];
    auto res = nearestOnSurface(surf, query);

    REQUIRE_THAT(res.distance, WithinAbs(0.0, 1e-4));
}

TEST_CASE("Nearest on surface: point outside sphere", "[nearest_point]") {
    // Sphere of radius 1; query = (3, 0, 0) → nearest ~ (1, 0, 0), distance ~ 2
    int n = 9, m = 7;
    double R = 1.0;
    auto lattice = sphereLattice(n, m, R);
    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    Vec3d query{3.0, 0.0, 0.0};
    auto res = nearestOnSurface(surf, query);

    REQUIRE_THAT(res.distance, WithinAbs(2.0, 0.15));
    // Result should lie on the sphere
    REQUIRE_THAT(norm(res.point), WithinAbs(R, 0.1));
}

TEST_CASE("Nearest on surface: flat plane, nearest to elevated point", "[nearest_point]") {
    // Flat z=0 plane; query = (2, 2, 5) → nearest = (2, 2, 0), distance = 5
    int n = 6, m = 6;
    std::vector<Vec3d> lattice(n * m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            lattice[j*n+i] = Vec3d{double(i), double(j), 0.0};

    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    Vec3d query{2.0, 2.0, 5.0};
    auto res = nearestOnSurface(surf, query);

    REQUIRE_THAT(res.distance, WithinAbs(5.0, 0.01));
    REQUIRE_THAT(res.point[2], WithinAbs(0.0, 1e-8));
}
