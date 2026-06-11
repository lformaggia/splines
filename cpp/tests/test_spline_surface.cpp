#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numbers>
#include <vector>
#include "splines/spline_surface.hpp"

using namespace splines;
using Catch::Matchers::WithinAbs;
namespace num = std::numbers;

// Helper: build an n×m flat lattice (z=0 bilinear plane)
static std::vector<Vec3d> flatLattice(int n, int m) {
    std::vector<Vec3d> pts(n * m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            pts[j*n + i] = Vec3d{double(i), double(j), 0.0};
    return pts;
}

// Helper: build an n×m lattice on a sphere of radius R
static std::vector<Vec3d> sphereLattice(int n, int m, double R = 1.0) {
    std::vector<Vec3d> pts(n * m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) {
            double u = num::pi * j / (m - 1);          // colatitude [0, pi]
            double v = 2.0 * num::pi * i / (n - 1);    // longitude  [0, 2pi]
            pts[j*n + i] = Vec3d{
                R * std::sin(u) * std::cos(v),
                R * std::sin(u) * std::sin(v),
                R * std::cos(u)
            };
        }
    return pts;
}

TEST_CASE("SplineSurface interpolation property on flat patch", "[spline_surface]") {
    int n = 5, m = 4;
    auto lattice = flatLattice(n, m);
    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    // eval at integer parameters must reproduce lattice knots
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) {
            auto pt = surf.eval(double(i), double(j));
            REQUIRE_THAT(pt.pos[0], WithinAbs(double(i), 1e-9));
            REQUIRE_THAT(pt.pos[1], WithinAbs(double(j), 1e-9));
            REQUIRE_THAT(pt.pos[2], WithinAbs(0.0,       1e-9));
        }
}

TEST_CASE("SplineSurface flat surface: normal is (0,0,1)", "[spline_surface]") {
    auto lattice = flatLattice(5, 5);
    auto surf = SplineSurface<3>::interpolate(lattice, 5, 5);

    // Check normal at several interior points
    for (double U : {0.5, 1.3, 2.7})
        for (double V : {0.5, 1.5, 2.9}) {
            auto pt = surf.eval(U, V);
            auto n  = pt.normal();
            REQUIRE_THAT(std::abs(n[2]), WithinAbs(1.0, 0.01));
        }
}

TEST_CASE("SplineSurface flat surface: position interpolates linearly", "[spline_surface]") {
    // On a flat z=0 plane, eval(U,V) should give (U,V,0)
    auto lattice = flatLattice(6, 6);
    auto surf = SplineSurface<3>::interpolate(lattice, 6, 6);

    for (double U : {0.0, 1.0, 2.5, 4.0, 5.0})
        for (double V : {0.0, 1.0, 2.5, 4.0, 5.0}) {
            auto pt = surf.eval(U, V);
            REQUIRE_THAT(pt.pos[0], WithinAbs(U,   1e-9));
            REQUIRE_THAT(pt.pos[1], WithinAbs(V,   1e-9));
            REQUIRE_THAT(pt.pos[2], WithinAbs(0.0, 1e-9));
        }
}

TEST_CASE("SplineSurface sphere: interpolation at knot positions", "[spline_surface]") {
    int n = 7, m = 5;
    double R = 2.0;
    auto lattice = sphereLattice(n, m, R);
    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) {
            auto pt = surf.eval(double(i), double(j));
            const auto& ref = lattice[j*n + i];
            REQUIRE_THAT(pt.pos[0], WithinAbs(ref[0], 1e-5));
            REQUIRE_THAT(pt.pos[1], WithinAbs(ref[1], 1e-5));
            REQUIRE_THAT(pt.pos[2], WithinAbs(ref[2], 1e-5));
        }
}

TEST_CASE("SplineSurface sphere: interior points stay near sphere", "[spline_surface]") {
    // With a coarse lattice the spline still approximates the sphere
    int n = 9, m = 7;
    double R = 1.0;
    auto lattice = sphereLattice(n, m, R);
    auto surf = SplineSurface<3>::interpolate(lattice, n, m);

    for (int ji = 0; ji < 4; ++ji)
        for (int ii = 0; ii < 4; ++ii) {
            double U = double(ii) * (n - 1) / 3.0;
            double V = double(ji) * (m - 1) / 3.0;
            auto pt = surf.eval(U, V);
            double r = norm(pt.pos);
            REQUIRE_THAT(r, WithinAbs(R, 0.1));
        }
}

TEST_CASE("SplineSurface 2D variant (planar curve surface)", "[spline_surface]") {
    // A 2D surface is just a 2D lattice
    int n = 4, m = 3;
    std::vector<Vec2d> lattice(n * m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i)
            lattice[j*n+i] = Vec2d{double(i)*0.5, double(j)*0.5};

    auto surf = SplineSurface<2>::interpolate(lattice, n, m);
    // Interpolation property
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) {
            auto pt = surf.eval(double(i), double(j));
            REQUIRE_THAT(pt.pos[0], WithinAbs(double(i)*0.5, 1e-9));
            REQUIRE_THAT(pt.pos[1], WithinAbs(double(j)*0.5, 1e-9));
        }
}
