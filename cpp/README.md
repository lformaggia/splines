# splines — C++20 Cubic Spline Library

A header-only C++20 library for parametric cubic spline curves and bicubic
tensor-product surfaces, including arc-length computation, differential
geometry, and nearest-point projection.  Port of the Fortran `libsplin`
library.

---

## Contents

```
cpp/
├── CMakeLists.txt
├── include/
│   └── splines/
│       ├── splines.hpp          # Master include (pulls in everything)
│       ├── concepts.hpp         # RealType concept
│       ├── vec.hpp              # Vec<Dim,Real> spatial vector
│       ├── end_conditions.hpp   # EndCondition enum + ImposedTangents
│       ├── tridiagonal.hpp      # Thomas algorithm (detail namespace)
│       ├── spline_curve.hpp     # SplineCurve<Dim,Real>
│       ├── arc_length.hpp       # Free-function arc-length wrappers
│       ├── geometry.hpp         # Curvature, torsion, surface normal, κ₁/κ₂
│       ├── spline_surface.hpp   # SplineSurface<Dim,Real>
│       └── nearest_point.hpp   # nearestOnCurve, nearestOnSurface
└── tests/
    ├── CMakeLists.txt
    ├── test_vec.cpp
    ├── test_tridiagonal.cpp
    ├── test_spline_curve.cpp
    ├── test_arc_length.cpp
    ├── test_spline_surface.cpp
    ├── test_geometry.cpp
    └── test_nearest_point.cpp
```

---

## Mathematical Background

### Parametrisation

Both curves and surfaces use **chord-length parametrisation**.

**Curves** — For a sequence of n knots, the global parameter t ∈ [0, n−1].
Integer values of t correspond exactly to knots; each interval [i, i+1] is
one cubic Hermite arc.  Internally the library stores tangent vectors as
`tau = dr/dc` (derivative with respect to chord length c).  Within arc i the
Hermite formula uses `dr/du = chords[i] * tau[i]`.

**Surfaces** — For an n×m lattice, global parameters are U ∈ [0, n−1] and
V ∈ [0, m−1].  Within patch (pi, pj) the local coordinates are
(u, v) = (U−pi, V−pj) ∈ [0,1]².

### Spline Construction

The tangent vectors at the knots satisfy a tridiagonal linear system.
The interior equations enforce C² continuity; the two boundary equations
are selected by the `EndCondition` enum:

| End condition | Description |
|---------------|-------------|
| `Natural`     | Zero second derivative at each endpoint (r'' = 0). |
| `Bessel`      | Tangent equals the derivative of the interpolating parabola through the first (or last) three knots. |
| `NotAKnot`    | Third derivative continuous at the second and second-to-last knot, making the first/last pair of arcs a single cubic. Requires ≥ 3 knots. |
| `Quadratic`   | Second derivative at each endpoint equals that at the adjacent interior knot. |
| `Imposed`     | User-supplied tangent vectors via `ImposedTangents`. |

The tridiagonal system is solved with the Thomas (LU + back-substitution)
algorithm (`tridiagonal.hpp`).

### Surface Construction (`SplineSurface`)

Construction uses four spline-fitting passes (port of Fortran `getta3`):

1. **Pass 1** — v-direction splines through each column → ∂r/∂v, v-chords.
2. **Pass 2** — u-direction splines through each row → ∂r/∂u, u-chords.
3. **Pass 3** — u-direction splines through the ∂r/∂v values on the two
   boundary rows → boundary ∂²r/∂u∂v.
4. **Pass 4** — v-direction splines through the ∂r/∂u values with imposed
   boundary ∂²r/∂u∂v → interior ∂²r/∂u∂v.

Each patch is then a bicubic polynomial with coefficients `A = CC * Q * CCᵀ`
where CC is the 4×4 Hermite-to-monomial matrix.

### Arc Length

Arc length is computed by adaptive Romberg quadrature of the speed
`|dr/du|` over the arc.  The integrand is the square root of a degree-4
polynomial in the local parameter u (port of Fortran `coeff.f` + `qrombs`).

### Nearest-Point Algorithms

| Function          | Algorithm |
|-------------------|-----------|
| `nearestOnCurve`  | Newton iteration on `(r(t)−q)·r'(t) = 0`. |
| `nearestOnSurface`| 2-D Newton on `[rᵤ·v, r_v·v] = 0` with approximate Hessian. |

Both functions fall back to a brute-force lattice/knot scan when no initial
guess is provided.

---

## Requirements

| Requirement | Version |
|-------------|---------|
| C++ compiler | C++20 support required (GCC ≥ 11, Clang ≥ 13, MSVC ≥ 19.29) |
| CMake | ≥ 3.20 |
| Catch2 | v3.5.2 (fetched automatically by CMake) |

No other external dependencies.  The library itself is header-only.

---

## Compilation

```bash
cd cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

For a debug build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)
```

The `splines` CMake target is an `INTERFACE` library; link against it in
your own project with:

```cmake
add_subdirectory(path/to/splines/cpp)
target_link_libraries(my_target PRIVATE splines)
```

---

## Running the Tests

```bash
cd cpp/build
ctest --output-on-failure
```

Or run individual test executables directly:

```bash
./tests/test_spline_curve
./tests/test_arc_length
./tests/test_spline_surface
./tests/test_geometry
./tests/test_nearest_point
```

### Test summary (49 tests, 7 binaries)

| Binary                  | What is tested |
|-------------------------|----------------|
| `test_vec`              | `Vec` arithmetic, dot, cross, norm, normalised |
| `test_tridiagonal`      | Thomas algorithm for scalar and vector RHS |
| `test_spline_curve`     | Interpolation property, C¹ continuity, circle, line, all end conditions |
| `test_arc_length`       | Straight line, semicircle, full circle, parabola (known exact values) |
| `test_spline_surface`   | Flat patch, sphere — interpolation and interior geometry |
| `test_geometry`         | Curvature (2-D/3-D), torsion (helix), surface normal, principal κ |
| `test_nearest_point`    | On-curve/on-surface projections, off-curve distance (circle, line, sphere, plane) |

---

## API Quick Reference

### `Vec<Dim, Real>` (`vec.hpp`)

```cpp
Vec<3, double> a{1.0, 2.0, 3.0};
Vec<3, double> b{4.0, 5.0, 6.0};

double d = dot(a, b);          // inner product
double n = norm(a);            // Euclidean norm
Vec3d  u = normalized(a);      // unit vector
Vec3d  c = cross(a, b);        // cross product (3-D only)
```

Convenience aliases: `Vec2d`, `Vec3d`, `Vec2f`, `Vec3f`.

---

### `SplineCurve<Dim, Real>` (`spline_curve.hpp`)

```cpp
#include "splines/spline_curve.hpp"
using namespace splines;

std::vector<Vec3d> pts = {{0,0,0},{1,2,0},{3,1,0},{4,3,1}};
auto curve = SplineCurve<3>::interpolate(pts, EndCondition::NotAKnot);

int    n  = curve.knotCount();    // number of knots
int    na = curve.arcCount();     // number of arcs = n-1
double L  = curve.totalLength();  // total arc length

auto pt = curve.eval(1.5);  // evaluate at t=1.5
Vec3d pos    = pt.pos;      // position
Vec3d d1     = pt.d1;       // dr/du
Vec3d d2     = pt.d2;       // d²r/du²
double speed = pt.speed;    // |dr/du|

Vec3d tang = curve.tangent(1.5);  // unit tangent

// Imposed end conditions
ImposedTangents<3> imp{Vec3d{1,0,0}, Vec3d{0,1,0}};
auto c2 = SplineCurve<3>::interpolate(pts, EndCondition::Imposed, imp);
```

---

### `SplineSurface<Dim, Real>` (`spline_surface.hpp`)

```cpp
#include "splines/spline_surface.hpp"
using namespace splines;

int n = 7, m = 5;
std::vector<Vec3d> lattice(n * m);
// ... fill lattice[j*n + i] = point at column i, row j ...

auto surf = SplineSurface<3>::interpolate(lattice, n, m);

auto pt = surf.eval(1.5, 2.3);
Vec3d pos  = pt.pos;   // position r(U,V)
Vec3d ru   = pt.ru;    // ∂r/∂U
Vec3d rv   = pt.rv;    // ∂r/∂V
Vec3d ruu  = pt.ruu;   // ∂²r/∂U²
Vec3d rvv  = pt.rvv;   // ∂²r/∂V²
Vec3d ruv  = pt.ruv;   // ∂²r/∂U∂V
Vec3d norm = pt.normal(); // unit outward normal (3-D only)
```

---

### `geometry.hpp`

```cpp
#include "splines/geometry.hpp"
using namespace splines;

auto pt = curve.eval(t);

// Curve geometry
double kappa = curvature(pt.d1, pt.d2);           // κ (2-D signed, 3-D unsigned)
double tau   = torsion(pt.d1, pt.d2, pt.d3);      // τ (3-D)

// Surface geometry
auto spt = surf.eval(U, V);
Vec3d n    = surfaceNormal(spt.ru, spt.rv);
auto [k1, k2] = principalCurvatures(spt.ru, spt.rv,
                                    spt.ruu, spt.ruv, spt.rvv);
```

---

### `nearest_point.hpp`

```cpp
#include "splines/nearest_point.hpp"
using namespace splines;

Vec3d query{3.0, 0.0, 0.0};

// Nearest point on curve
auto cr = nearestOnCurve(curve, query);
double t    = cr.param;
Vec3d  p    = cr.point;
double dist = cr.distance;
bool   ok   = cr.converged;

// Nearest point on surface  (provide initial guess to skip lattice scan)
auto sr = nearestOnSurface(surf, query, 2.0, 1.5);
double U = sr.U, V = sr.V;
```

---

## Generating HTML Documentation

The code is documented with [Doxygen](https://www.doxygen.nl/).  To generate
HTML documentation:

```bash
cd cpp
doxygen Doxyfile   # if a Doxyfile is present
# or generate a default one first:
doxygen -g Doxyfile
# edit Doxyfile: set INPUT = include/splines, RECURSIVE = YES,
#   EXTRACT_ALL = YES, GENERATE_LATEX = NO
doxygen Doxyfile
# open html/index.html
```

---

## Notes on the Fortran Port

The Fortran `libsplin` library used pre-IEEE-754 floating-point semantics
where `0/0 = 0`.  Several routines relied on this to simplify the not-a-knot
boundary condition assembly.  This C++20 port replaces those degenerate
formulae with mathematically equivalent but numerically sound alternatives:

- **`NotAKnot` start row** — the Fortran formula creates a zero modified
  diagonal after LU factorisation.  The C++ version derives a non-degenerate
  first row by eliminating τ[1] via the C² equation at knot 1.
- **`NotAKnot` end row** — symmetric fix at the last knot.
- **`NotAKnot` with n = 3** — both boundary conditions encode the same
  constraint; the library falls back to `Natural` automatically.
