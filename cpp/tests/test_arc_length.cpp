#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numbers>
#include <vector>
#include "splines/spline_curve.hpp"
#include "splines/arc_length.hpp"

using namespace splines;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
namespace num = std::numbers;

// Helper: N equally-spaced points on circle plus closing knot (N arcs total)
static std::vector<Vec2d> circleKnots(int N, double R = 1.0) {
    std::vector<Vec2d> pts(N + 1);
    for (int i = 0; i < N; ++i) {
        double theta = 2.0 * num::pi * i / N;
        pts[i] = Vec2d{R * std::cos(theta), R * std::sin(theta)};
    }
    pts[N] = pts[0];  // close the curve
    return pts;
}

TEST_CASE("Arc length: straight line segment", "[arc_length]") {
    // Two knots: (0,0) → (3,4): exact length = 5
    std::vector<Vec2d> pts = {{0.0, 0.0}, {3.0, 4.0}};
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::Natural);

    double L = curve.arcLength(0, 1e-8);
    REQUIRE_THAT(L, WithinAbs(5.0, 0.01));
}

TEST_CASE("Arc length: semicircle approximation", "[arc_length]") {
    // 13 knots on upper semicircle (0..pi): exact arc length = pi
    int N = 13;
    std::vector<Vec2d> pts(N);
    for (int i = 0; i < N; ++i) {
        double theta = num::pi * i / (N - 1);
        pts[i] = Vec2d{std::cos(theta), std::sin(theta)};
    }
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::Bessel);

    double L = curve.totalLength(1e-6);
    REQUIRE_THAT(L, WithinAbs(num::pi, 0.01));   // within 1%
}

TEST_CASE("Arc length: full unit circle", "[arc_length]") {
    // 16 equally-spaced points on unit circle: total length ≈ 2*pi
    auto pts = circleKnots(16);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    double L = curve.totalLength(1e-6);
    REQUIRE_THAT(L, WithinAbs(2.0 * num::pi, 0.05));   // within ~1%
}

TEST_CASE("Arc length: parabola arc known formula", "[arc_length]") {
    // r(t) = (t, t^2/2) for t in [0, 1].
    // Exact arc length = integral_0^1 sqrt(1 + t^2) dt
    //                  = [t/2*sqrt(1+t^2) + 1/2*ln(t + sqrt(1+t^2))]_0^1
    //                  = sqrt(2)/2 + ln(1+sqrt(2))/2 ≈ 1.1478943...
    double exact = std::sqrt(2.0)/2.0 + std::log(1.0 + std::sqrt(2.0))/2.0;

    int N = 20;
    std::vector<Vec2d> pts(N + 1);
    for (int i = 0; i <= N; ++i) {
        double t = double(i) / N;
        pts[i] = Vec2d{t, t*t / 2.0};
    }

    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::Bessel);
    double L = curve.totalLength(1e-7);
    REQUIRE_THAT(L, WithinAbs(exact, 0.005));
}

TEST_CASE("Arc length: circle with radius R", "[arc_length]") {
    double R = 3.5;
    auto pts = circleKnots(20, R);
    auto curve = SplineCurve<2>::interpolate(pts, EndCondition::NotAKnot);

    double L = curve.totalLength(1e-6);
    REQUIRE_THAT(L, WithinAbs(2.0 * num::pi * R, 0.2));  // within ~1%
}
