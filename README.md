# splines

This repository contains two related implementations of a spline-geometry toolkit:

- a legacy Fortran library centered on cubic parametric curves and bicubic tensor-product surfaces,
- a modern C++20 header-only port of the same mathematical core,
- tests, debugging programs, and notes produced while validating and repairing the original Fortran code.

The codebase is aimed at scientific and geometric computing tasks such as interpolation of point data, curve and surface evaluation, curvature and derivative analysis, arc-length computation, and nearest-point / point-location queries.

## Repository Overview

At the top level, the repository is organized around the original Fortran implementation and the newer C++ port:

```text
splines/
├── cpp/              C++20 header-only spline library and unit tests
├── docs/             Reserved documentation area (currently empty)
├── fortran/          Legacy Fortran spline library, test programs, and fix notes
├── tools/            Reserved tooling area (currently empty)
├── build-aux/        Autotools helper files for legacy build support
├── autom4te.cache/   Autotools cache artifacts
└── README.md         Repository-level guide
```

The `docs/` and `tools/` directories exist but currently do not contain tracked files. Most of the active project content is in `fortran/` and `cpp/`.

## Main Components

### `fortran/`

This directory contains the original spline library, written primarily in fixed-form Fortran 77 style with a few newer Fortran test/debug drivers.

Important contents:

- `Makefile.in`, `configure.in`, `myrules.make`: legacy Autoconf and make-based build system.
- `*.f`: core spline, geometry, conversion, evaluation, and point-location routines.
- `test_splines.f90`: main Fortran regression/debug test program.
- `debug_*.f90`: focused diagnostic programs for tangent and arc-length issues.
- `README.md`: detailed library-level description of the Fortran routines and data representations.
- `INDEX.txt`, `TEST_REPORT.txt`, `TESTING_SUMMARY.txt`: testing guides and reports.
- `fixes/`: investigation notes and bug documentation.
- `new/`: alternate or revised implementations of a few routines, including parametrization helpers and spline builders.

The Fortran build produces a static library:

- `libsplin.a`

Based on `myrules.make`, the library is assembled from routines covering:

- spline construction and parametrization,
- tridiagonal linear-system solves,
- curve evaluation and curvature,
- bicubic surface interpolation and evaluation,
- explicit polynomial patch and arc representations,
- arc-length calculations,
- nearest-point / closest-point searches,
- surface edge, iso-curve, and coordinate conversion utilities.

Representative Fortran source files include:

- `psplin.f`, `pspl.f`, `pspli2.f`: spline construction routines,
- `evtan.f`, `trcoe.f`, `tridi.f`: tangent and linear-system infrastructure,
- `getp2.f`, `gettan.f`, `getk.f`, `slen.f`: curve evaluation and differential geometry,
- `getta3.f`, `surfb.f`, `evapa2.f`, `evapa3.f`: surface interpolation/evaluation,
- `locpfs.f`, `locsur.f`, `locsu2.f`, `locsug.f`, `locps3.f`, `locps4.f`: point-location routines,
- `evps1d.f`, `evsurg.f`, `patc1d.f`, `gpsur.f`: explicit polynomial curve/surface forms.

### `cpp/`

This directory contains a modern C++20 port of the spline functionality. The C++ library is header-only and is built with CMake.

Important contents:

- `CMakeLists.txt`: top-level CMake configuration for the C++ port.
- `include/splines/`: public headers implementing the library.
- `tests/`: unit tests for vectors, tridiagonal solves, curves, surfaces, geometry, arc length, and nearest-point queries.
- `README.md`: detailed C++-specific overview, mathematical notes, API summary, and build instructions.

The public headers include:

- `splines.hpp`: umbrella include,
- `vec.hpp`: fixed-size vector type and operations,
- `end_conditions.hpp`: spline endpoint-condition definitions,
- `tridiagonal.hpp`: Thomas solver implementation,
- `spline_curve.hpp`: cubic spline curve API,
- `spline_surface.hpp`: bicubic spline surface API,
- `arc_length.hpp`, `geometry.hpp`, `nearest_point.hpp`: analysis and query utilities.

The CMake target exported by this subproject is an `INTERFACE` library named `splines`.

## What The Software Does

Across both implementations, the repository focuses on geometric splines rather than general-purpose numerical splines. The core functionality includes:

- interpolating curves through knot points using cubic parametric splines,
- generating tangent vectors under several endpoint conditions,
- evaluating positions, first and higher derivatives, tangents, and curvature,
- estimating or integrating arc length,
- interpolating tensor-product bicubic surfaces over rectangular point lattices,
- evaluating surface points, partial derivatives, normals, and curvature-related quantities,
- converting between interpolation-oriented and explicit polynomial representations,
- finding closest points on curves and surfaces.

The Fortran implementation also includes application-oriented helpers for packed storage layouts and geometry database I/O conventions inherited from the original system.

## Build Paths

### Building the Fortran library

From `fortran/`:

```sh
autoreconf -fi
./configure
make
```

This generates `libsplin.a` using the object list defined in `fortran/myrules.make`.

A test executable can also be built via the makefile target:

```sh
make test_splines
```

### Building the C++ library and tests

From `cpp/`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Run the C++ test suite with:

```sh
cd build
ctest --output-on-failure
```

## Testing and Project Status

The repository contains both a mature C++ test suite and an in-progress validation effort for the legacy Fortran code.

### C++ status

The `cpp/tests/` directory contains dedicated test binaries covering:

- vector algebra,
- tridiagonal solves,
- spline curves,
- spline surfaces,
- arc length,
- geometry utilities,
- nearest-point algorithms.

The C++ side appears to be structured as the cleaned-up, modernized implementation of the spline algorithms.

### Fortran status

The `fortran/` directory includes explicit test reports and fix investigations. Those files document that the legacy library builds, but some numerical and correctness issues remain under investigation, including:

- a documented indexing bug in `trcoe.f` for quadratic boundary conditions,
- arc-length discrepancies,
- NaN/Infinity failures in some tangent-related paths,
- incomplete or failing surface-related cases in the current test programs.

For the current state of that work, the most relevant files are:

- `fortran/TEST_REPORT.txt`
- `fortran/TESTING_SUMMARY.txt`
- `fortran/fixes/FIXES_APPLIED.txt`
- `fortran/fixes/fix_trcoe_quadratic_bug.txt`
- `fortran/fixes/fix_arc_length_investigation.txt`
- `fortran/fixes/fix_nan_infinity_investigation.txt`

## Recommended Entry Points

If you are new to the repository, the most useful starting points are:

1. `fortran/README.md` for the detailed description of the original library and routine catalog.
2. `cpp/README.md` for the modern API, mathematical background, and build/test instructions.
3. `fortran/INDEX.txt` for the testing artifacts and debugging notes around the legacy implementation.

## Summary

In practical terms, this repository is both:

- a preservation and repair effort for a legacy spline-geometry Fortran library, and
- a modern C++20 reimplementation of the same core computational ideas.

If you want the original algorithms and packed routine set, work in `fortran/`. If you want a cleaner interface and contemporary build/test flow, start in `cpp/`.
