#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numbers>
#include <vector>
#include "splines/spline_curve.hpp"

using namespace splines;
using Catch::Matchers::WithinAbs;
namespace num = std::numbers;

// Helper: sample N equally-spaced points on the unit circle
static std::vector<Vec2d> circleKnots(int N) {
    std::vector<Vec2d> pts(N);
    for (int i = 0; i < N; ++i) {
        double theta = 2.0 * num::pi * i / N;
        pts[i] = Vec2d{std::cos(theta), std::sin(theta)};
    }
    return pts;
}

// Helper: sample N+1 equally-spaced points on a line segment
static std::vector<Vec2d> lineKnots(int N) {
    std::vector<Vec2d> pts(N + 1);
    for (int i = 0; i <= N; ++i)
        pts[i] = Vec2d{double(i), 0.0};
    return pts;
}

TEST_CASE("SplineCurve interpolation property", "[spline_curve]") {
    // The spline must pass through every knot
    auto pts = circleKnots(8);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    for (int i = 0; i < curve.knotCount(); ++i) {
        auto p = curve.position(double(i));
        REQUIRE_THAT(p[0], WithinAbs(pts[i][0], 1e-10));
        REQUIRE_THAT(p[1], WithinAbs(pts[i][1], 1e-10));
    }
}

TEST_CASE("SplineCurve smoothness: C1 continuity at interior knots", "[spline_curve]") {
    auto pts = circleKnots(10);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    // Evaluate first derivative from left and right of each interior knot
    double eps = 1e-7;
    for (int i = 1; i < curve.knotCount() - 1; ++i) {
        auto left  = curve.eval(double(i) - eps);
        auto right = curve.eval(double(i) + eps);
        // Tangent directions should agree (C1)
        REQUIRE_THAT(left.d1[0], WithinAbs(right.d1[0], 1e-4));
        REQUIRE_THAT(left.d1[1], WithinAbs(right.d1[1], 1e-4));
    }
}

TEST_CASE("SplineCurve points on unit circle stay near circle", "[spline_curve]") {
    // With enough knots and NotAKnot conditions, the spline of circle points
    // should stay very close to the unit circle.
    auto pts = circleKnots(12);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    for (int k = 0; k <= 100; ++k) {
        double t = double(k) / 100.0 * double(curve.arcCount());
        auto p = curve.position(t);
        double r = norm(p);
        // With 12 knots the spline approximates the circle well
        REQUIRE_THAT(r, WithinAbs(1.0, 0.01));
    }
}

TEST_CASE("SplineCurve straight line: position correct", "[spline_curve]") {
    auto pts = lineKnots(5);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    // Any point should lie on y=0
    for (int k = 0; k <= 50; ++k) {
        double t = double(k) / 50.0 * double(curve.arcCount());
        auto p = curve.position(t);
        REQUIRE_THAT(p[1], WithinAbs(0.0, 1e-12));
        // x should increase linearly
        REQUIRE_THAT(p[0], WithinAbs(t, 1e-12));
    }
}

TEST_CASE("SplineCurve endpoint tangent modes", "[spline_curve]") {
    auto pts = circleKnots(8);

    // All supported end conditions should produce a valid curve
    for (auto ec : {EndCondition::Natural, EndCondition::Bessel,
                    EndCondition::NotAKnot, EndCondition::Quadratic}) {
        auto curve = SplineCurve<2>::interpolate(pts, ec);
        REQUIRE(curve.knotCount() == 8);
        REQUIRE(curve.arcCount()  == 7);

        // Interpolation must hold for all end conditions
        for (int i = 0; i < curve.knotCount(); ++i) {
            auto p = curve.position(double(i));
            REQUIRE_THAT(p[0], WithinAbs(pts[i][0], 1e-10));
            REQUIRE_THAT(p[1], WithinAbs(pts[i][1], 1e-10));
        }
    }
}

TEST_CASE("SplineCurve imposed tangent at endpoints", "[spline_curve]") {
    // At a circle knot at angle 0 (point (1,0)), the tangent should be (0,1)
    auto pts = circleKnots(8);
    ImposedTangents<2> imp;
    imp.start = Vec2d{0.0, 1.0 / 1.0};  // dr/dc at start knot; chord ≈ 1, so tangent ≈ (0,1)
    imp.end   = Vec2d{0.0, 1.0 / 1.0};

    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::Imposed, imp);
    REQUIRE(curve.knotCount() == 8);
    auto t_start = curve.tangents()[0];
    REQUIRE_THAT(t_start[0], WithinAbs(0.0, 1e-10));
}

TEST_CASE("SplineCurve 3D straight line", "[spline_curve]") {
    std::vector<Vec3d> pts;
    for (int i = 0; i <= 4; ++i)
        pts.push_back(Vec3d{double(i), double(i), 0.0});

    auto curve = SplineCurve<3>::interpolate(pts, EndCondition::Natural);

    for (int k = 0; k <= 40; ++k) {
        double t = double(k) / 40.0 * double(curve.arcCount());
        auto p = curve.position(t);
        REQUIRE_THAT(p[0], WithinAbs(t, 1e-12));
        REQUIRE_THAT(p[1], WithinAbs(t, 1e-12));
        REQUIRE_THAT(p[2], WithinAbs(0.0, 1e-12));
    }
}
