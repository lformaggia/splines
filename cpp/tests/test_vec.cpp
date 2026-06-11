#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "splines/vec.hpp"

using namespace splines;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Vec2 construction and indexing", "[vec]") {
    Vec2d v{1.0, 2.0};
    REQUIRE(v[0] == 1.0);
    REQUIRE(v[1] == 2.0);
}

TEST_CASE("Vec3 construction and indexing", "[vec]") {
    Vec3d v{3.0, -1.5, 4.2};
    REQUIRE(v[0] == 3.0);
    REQUIRE(v[1] == -1.5);
    REQUIRE(v[2] == 4.2);
}

TEST_CASE("Vec addition", "[vec]") {
    Vec3d a{1, 2, 3}, b{4, 5, 6};
    auto c = a + b;
    REQUIRE_THAT(c[0], WithinAbs(5.0, 1e-15));
    REQUIRE_THAT(c[1], WithinAbs(7.0, 1e-15));
    REQUIRE_THAT(c[2], WithinAbs(9.0, 1e-15));
}

TEST_CASE("Vec subtraction", "[vec]") {
    Vec3d a{4, 5, 6}, b{1, 3, 5};
    auto c = a - b;
    REQUIRE_THAT(c[0], WithinAbs(3.0, 1e-15));
    REQUIRE_THAT(c[1], WithinAbs(2.0, 1e-15));
    REQUIRE_THAT(c[2], WithinAbs(1.0, 1e-15));
}

TEST_CASE("Vec scalar multiply", "[vec]") {
    Vec3d v{1, 2, 3};
    auto u = v * 2.0;
    REQUIRE_THAT(u[0], WithinAbs(2.0, 1e-15));
    REQUIRE_THAT(u[1], WithinAbs(4.0, 1e-15));
    REQUIRE_THAT(u[2], WithinAbs(6.0, 1e-15));
    // Commutative
    auto w = 3.0 * v;
    REQUIRE_THAT(w[2], WithinAbs(9.0, 1e-15));
}

TEST_CASE("Vec unary minus", "[vec]") {
    Vec3d v{1, -2, 3};
    auto u = -v;
    REQUIRE_THAT(u[0], WithinAbs(-1.0, 1e-15));
    REQUIRE_THAT(u[1], WithinAbs(2.0,  1e-15));
    REQUIRE_THAT(u[2], WithinAbs(-3.0, 1e-15));
}

TEST_CASE("Vec dot product", "[vec]") {
    Vec3d a{1, 2, 3}, b{4, 5, 6};
    double d = dot(a, b);     // 4+10+18 = 32
    REQUIRE_THAT(d, WithinAbs(32.0, 1e-14));
}

TEST_CASE("Vec norm", "[vec]") {
    Vec3d v{3, 4, 0};
    REQUIRE_THAT(norm(v), WithinAbs(5.0, 1e-14));
    REQUIRE_THAT(norm_sq(v), WithinAbs(25.0, 1e-14));
}

TEST_CASE("Vec normalized", "[vec]") {
    Vec3d v{0, 3, 4};
    auto u = normalized(v);
    REQUIRE_THAT(norm(u), WithinAbs(1.0, 1e-14));
    REQUIRE_THAT(u[0], WithinAbs(0.0,   1e-14));
    REQUIRE_THAT(u[1], WithinAbs(0.6,   1e-14));
    REQUIRE_THAT(u[2], WithinAbs(0.8,   1e-14));
}

TEST_CASE("Vec3 cross product", "[vec]") {
    Vec3d x{1, 0, 0}, y{0, 1, 0};
    auto z = cross(x, y);
    REQUIRE_THAT(z[0], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(z[1], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(z[2], WithinAbs(1.0, 1e-15));

    // Anti-commutative
    auto neg_z = cross(y, x);
    REQUIRE_THAT(neg_z[2], WithinAbs(-1.0, 1e-15));
}

TEST_CASE("Vec cross of parallel vectors is zero", "[vec]") {
    Vec3d a{1, 2, 3}, b{2, 4, 6};
    auto c = cross(a, b);
    REQUIRE_THAT(norm(c), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Vec equality", "[vec]") {
    Vec2d a{1.0, 2.0}, b{1.0, 2.0}, c{1.0, 3.0};
    REQUIRE(a == b);
    REQUIRE(a != c);
}

TEST_CASE("Vec compound assignment", "[vec]") {
    Vec3d v{1, 2, 3};
    v += Vec3d{1, 1, 1};
    REQUIRE_THAT(v[0], WithinAbs(2.0, 1e-15));
    v *= 2.0;
    REQUIRE_THAT(v[1], WithinAbs(6.0, 1e-15));
}
