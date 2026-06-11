#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>
#include "splines/tridiagonal.hpp"
#include "splines/vec.hpp"

using namespace splines;
using namespace splines::detail;
using Catch::Matchers::WithinAbs;

// Solve a known 4×4 tridiagonal system and check the solution.
//
// Equation:  [2 1 0 0] [x0]   [3]
//            [1 2 1 0] [x1] = [6]
//            [0 1 2 1] [x2]   [6]
//            [0 0 1 2] [x3]   [5]
//
// Solution (verified by substitution): x = [1, 1, 2, 1.5]  -- let's compute:
// From last row: x2 + 2*x3 = 5.
// From 3rd:  x1 + 2*x2 + x3 = 6.
// From 2nd:  x0 + 2*x1 + x2 = 6.
// From 1st:  2*x0 + x1 = 3.
// Trial x0=1, x1=1: 2+1=3 ✓.  x0+2x1+x2=1+2+x2=6 → x2=3. x1+2x2+x3=1+6+x3=6 → x3=-1.
// x2+2x3=3-2=1 ≠ 5. Let me recompute.
//
// Actually solve properly.  Let α,β,γ,δ be unknowns:
//   2α + β = 3               → β = 3 - 2α
//   α + 2β + γ = 6           → α + 2(3-2α) + γ = 6 → γ = 3α
//   β + 2γ + δ = 6           → (3-2α) + 6α + δ = 6 → δ = 3 - 4α
//   γ + 2δ = 5               → 3α + 2(3-4α) = 5 → -5α = -1 → α = 1/5 = 0.2
//
//   α = 0.2, β = 2.6, γ = 0.6, δ = 2.2

TEST_CASE("Tridiagonal factorize and solve (scalar)", "[tridiagonal]") {
    std::vector<double> a = {0.0, 1.0, 1.0, 1.0};
    std::vector<double> b = {2.0, 2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0, 1.0, 0.0};
    std::vector<double> rhs_d = {3.0, 6.0, 6.0, 5.0};

    // Wrap in Vec<1> to use the template (trivial 1-D)
    using V = Vec<1, double>;
    std::vector<V> d = {{rhs_d[0]}, {rhs_d[1]}, {rhs_d[2]}, {rhs_d[3]}};

    tridiag_factor<double>(a, b, c);
    tridiag_solve<double>(std::span<const double>{a},
                          std::span<const double>{b},
                          std::span<const double>{c},
                          std::span{d});

    REQUIRE_THAT(d[0][0], WithinAbs(0.2, 1e-12));
    REQUIRE_THAT(d[1][0], WithinAbs(2.6, 1e-12));
    REQUIRE_THAT(d[2][0], WithinAbs(0.6, 1e-12));
    REQUIRE_THAT(d[3][0], WithinAbs(2.2, 1e-12));
}

TEST_CASE("Tridiagonal vector RHS (Vec<3>)", "[tridiagonal]") {
    // Same system but 3-D RHS: each component is independent.
    std::vector<double> a = {0.0, 1.0, 1.0, 1.0};
    std::vector<double> b = {2.0, 2.0, 2.0, 2.0};
    std::vector<double> c = {1.0, 1.0, 1.0, 0.0};

    std::vector<Vec3d> d = {
        {3.0,  6.0,  3.0},
        {6.0, 12.0,  6.0},
        {6.0, 12.0,  6.0},
        {5.0, 10.0,  5.0}
    };

    tridiag_factor<double>(a, b, c);
    tridiag_solve<double>(std::span<const double>{a},
                          std::span<const double>{b},
                          std::span<const double>{c},
                          std::span{d});

    // Component 0: same as above
    REQUIRE_THAT(d[0][0], WithinAbs(0.2, 1e-12));
    REQUIRE_THAT(d[3][0], WithinAbs(2.2, 1e-12));

    // Component 1: RHS is 2x component 0 → solution is 2x
    REQUIRE_THAT(d[0][1], WithinAbs(0.4, 1e-12));
    REQUIRE_THAT(d[3][1], WithinAbs(4.4, 1e-12));

    // Component 2: same as component 0
    REQUIRE_THAT(d[0][2], WithinAbs(0.2, 1e-12));
}

TEST_CASE("Tridiagonal 2×2 system", "[tridiagonal]") {
    // [3 1] [x]   [5]     x=(5-1y)/3,  y=(7-1x)/2
    // [1 2] [y] = [7]     y = (7 - (5-y)/3)/2 → 6y = 21 - 5 + y → 5y=16 → y=3.2, x=0.6
    std::vector<double> a = {0.0, 1.0};
    std::vector<double> b = {3.0, 2.0};
    std::vector<double> c = {1.0, 0.0};
    using V = Vec<1, double>;
    std::vector<V> d = {{5.0}, {7.0}};

    tridiag_factor<double>(a, b, c);
    tridiag_solve<double>(std::span<const double>{a},
                          std::span<const double>{b},
                          std::span<const double>{c},
                          std::span{d});

    REQUIRE_THAT(d[0][0], WithinAbs(0.6, 1e-12));
    REQUIRE_THAT(d[1][0], WithinAbs(3.2, 1e-12));
}

TEST_CASE("Tridiagonal 1×1 trivial system", "[tridiagonal]") {
    std::vector<double> a = {0.0};
    std::vector<double> b = {5.0};
    std::vector<double> c = {0.0};
    using V = Vec<1, double>;
    std::vector<V> d = {{10.0}};

    tridiag_factor<double>(a, b, c);
    tridiag_solve<double>(std::span<const double>{a},
                          std::span<const double>{b},
                          std::span<const double>{c},
                          std::span{d});

    REQUIRE_THAT(d[0][0], WithinAbs(2.0, 1e-12));
}
