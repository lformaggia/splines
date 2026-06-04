# splines

This workspace contains a legacy Fortran spline geometry library centered on:

- parametric cubic spline curves,
- tensor-product bicubic spline surfaces,
- conversion between spline, polynomial, Ferguson, and Bezier forms,
- geometric evaluation routines,
- closest-point / point-location routines,
- simple geometry database I/O used by the original application stack.

The build currently produces a static library:

- `libsplin.a`

The source style is fixed-form Fortran 77 with extensive inline documentation from the original codebase.

## Workspace Layout

- `*.f`: primary library routines.
- `new/*.f`: alternate or newer variants of a few spline/surface builders and parametrization helpers.
- `configure.in`, `configure`, `Makefile.in`, `myrules.make`: Autoconf and make-based build system.
- `libsplin.a`: static archive produced by `make`.

## Build

The library now builds with a normal modern toolchain such as `gfortran`.

```sh
autoreconf -fi
./configure
make
```

This generates `libsplin.a` from the object list defined in `myrules.make`.

## What The Library Does

At a high level, the code supports two related geometry pipelines.

### 1. Curve pipeline

Input data:

- knot points `q`
- tangent vectors `t`
- chord lengths `cs`

Core operations:

- build cubic parametric splines through points,
- evaluate points, tangents, curvature, and arc length,
- convert spline arcs into explicit polynomial coefficient form,
- find the closest point on a curve to a given spatial point.

### 2. Surface pipeline

Input data:

- lattice of surface knot points,
- tangent/cross-derivative arrays,
- chord-length scaling data,
- or explicit polynomial patch coefficients.

Core operations:

- interpolate bicubic tensor-product surfaces,
- convert patches between internal forms,
- evaluate points, derivatives, normals, and curvature,
- extract edge and iso-parametric curves,
- find the closest point on a surface to a given spatial point.

## Internal Geometry Representations

Understanding the routine list is much easier if you keep the main data structures in mind.

### Cubic spline curve form

Common arrays:

- `q(ndimn,*)`: knot coordinates.
- `t(ndimn,*)`: tangent vectors at knots.
- `cs(*)`: chord lengths or parametrization increments.

This form is used by routines such as `psplin`, `pspl`, `getp2`, `gettan`, `getk`, `slen`, and `locpfs`.

### Polynomial curve form

Common arrays:

- `xln1d(3,*)`: packed polynomial coefficients for curve arcs.
- `isn(*)`: pointers into `xln1d`.

This is the explicit arc-by-arc polynomial representation used by `evps1d`, `patc1d`, `xmin1d`, `xlenst`, `fundis`, and database routines.

### Surface interpolation form

Common arrays:

- `coor(ndimn,*)`: grid of knot coordinates.
- `tanret(ndimn,3,n,m)`: surface derivatives at knots.
- `choret(2,n,m)`: chord-length scaling data.

This form is produced by `getta3` and used by `surfb`, `evapa2`, `locsur`, `locsu2`, and `locsug`.

### Polynomial patch surface form

Common arrays:

- `aps(3,*)`: packed patch coefficients and metadata.
- `keypa(*)`: pointers to patches inside `aps`.

This is the explicit polynomial patch representation used by `evsurg`, `gpsur`, `gders1`, `locps3`, `locps4`, `distps`, `finded`, `getiso`, and the database readers/writers.

### Global versus local parametric coordinates

Several routines use a mixed parameter convention:

- the integer part identifies the arc or patch index,
- the fractional part is the local coordinate inside that arc or patch.

Helpers such as `loc2gl` and `gtlcor` convert between local and global parameter descriptions.

## Build System Notes

- `configure.in` detects the Fortran compiler and fills the legacy make variables expected by `Makefile.in`.
- `Makefile.in` builds object files and archives them into `libsplin.a`.
- `myrules.make` is the authoritative list of objects that belong to the library.
- The makefiles still contain some legacy install-time paths from the original environment, but they do not block library compilation.

## Routine Catalog

The list below covers the primary routines in the workspace. When a file has weak or missing structured documentation, the description is inferred from the implementation and nearby comments.

### Curve construction and parametrization

#### `psplin`

Builds a cubic parametric spline through knot points. It computes chord lengths, evaluates or inserts endpoint tangents, solves the tridiagonal system for knot tangents, and can also compute arc lengths.

#### `pspl`

Variant of `psplin` that supports an arbitrary knot sequence rather than relying only on the default parametrization path.

#### `pspli2`

Spline interpolation helper for scalar or component-wise data stored in `q`. It reuses the cubic spline machinery but is documented as the variant for interpolation of variable values rather than direct geometric coordinates.

#### `cholen`

Computes the chord length of each spline arc. Used as the first parametrization step before tangent estimation.

#### `evtan`

Computes knot tangents for a parametric cubic spline by solving the spline system. This is one of the central curve interpolation routines.

#### `trcoe`

Builds the tridiagonal coefficient arrays used by `psplin` / `pspli2` before the actual solve.

#### `trif`

First routine inside `tridi.f`. Performs the decomposition phase for a tridiagonal linear system.

#### `tris`

Second routine inside `tridi.f`. Performs back-substitution for one right-hand side after `trif`.

#### `tris2`

Third routine inside `tridi.f`. Multi-right-hand-side version of `tris`, useful for vector-valued geometry systems.

#### `gtimp`

Copies endpoint geometric information into the imposed-tangent work array used by spline construction.

#### `putimp`

Small helper that copies tangent values into the `timp` storage used for imposed end conditions.

#### `movt`

Surface analogue of `gtimp` used by `getta3` to place imposed derivative values into working storage.

#### `new/uniform`

Computes a uniform parametrization for spline interpolation in the `new/` branch.

#### `new/centripet`

Computes centripetal parametrization increments for spline interpolation.

#### `new/correctt`

Adjusts imposed end tangents so their effect is less sensitive to the chosen parametrization.

### Curve evaluation and differential geometry

#### `getp2`

Evaluates a spline curve at a parameter value and returns position plus several derivatives.

#### `gettan`

Computes the first derivative vector at a point on a cubic spline curve.

#### `gett2`

Computes the unit tangent at a knot by normalizing the derivative.

#### `getpoi`

Evaluates point position, first derivative, and an estimated curvature on a cubic spline curve.

#### `getk`

Computes curvature on a cubic spline.

#### `getk2`

Computes curvature, torsion, normal, and binormal from derivative data at a curve point.

#### `curva`

Estimates curvature at spline knots based on cubic spline data.

#### `coeff`

Builds the quartic polynomial coefficients for `|x_u|^2`, used in arc-length integration.

#### `slen`

Computes spline arc lengths using a Romberg-style integration process.

#### `lengt`

Numerical integration support library used by `slen`. It contains adaptive Romberg-style integration helpers derived from Numerical Recipes-era algorithms.

#### `xlenst`

Returns the length of a polynomial curve segment between parameter values `u1` and `u2`.

#### `patc1d`

Extracts the polynomial coefficient array for one cubic spline arc.

#### `evps1d`

Evaluates a polynomial curve in explicit coefficient form and returns position and derivatives in 3D.

### Curve nearest-point and search routines

#### `fundis`

Objective function for curve closest-point search. It measures the distance between a spatial point and a point on a polynomial curve.

#### `xmin1d`

One-dimensional minimizer used to locate the nearest point on a polynomial curve. It minimizes `fundis`.

#### `brasf`

Brackets spline arcs that may contain a minimum of the squared distance to a query point. It is a cheap preprocessing step before a more refined closest-point solve.

#### `locpfs`

Newton-style closest-point search on a Ferguson/cubic spline curve. It iterates in global curve parameter space and reports whether it found a point on the curve, a local minimum, or hit a boundary/extremum.

### Surface construction and patch conversion

#### `getta3`

Major surface interpolation routine. It builds a bicubic tensor-product spline surface through a lattice of points, supports several boundary-condition types, and computes knot derivatives and scaling data.

#### `surfb`

Builds a natural bicubic surface through a lattice of points and stores it directly in the general packed surface representation.

#### `evapa2`

Converts a bicubic spline patch from interpolation-form data (`tanret`, `choret`, `coor`) into explicit polynomial patch coefficients.

#### `evapa3`

Variation of `evapa2` that fills the coefficient matrix `A` directly for a patch. Useful when code wants a local patch coefficient tensor instead of the packed database-like format.

#### `evapat`

Builds polynomial patch coefficients starting from the Ferguson tensor-product representation.

#### `evaqfg`

Builds the Ferguson `Q` matrix from the internal bicubic patch representation.

#### `evppar`

Older conversion routine kept for compatibility. The header explicitly marks it as obsolete and recommends `evapa2` or `evapa3`.

#### `fertob`

Converts a Ferguson tensor-product patch into Bezier control points.

#### `buiaps`

Converts a bicubic spline surface patch representation into the library’s general polynomial surface representation.

#### `buicur`

Converts a cubic spline curve into the general polynomial curve representation.

### Surface evaluation and differential geometry

#### `gpsur`

Evaluates the position vector on a bicubic surface patch.

#### `evsurg`

Evaluates a polynomial surface patch and returns position and derivatives.

#### `gders1`

Computes first and second derivatives of a surface patch at a point.

#### `evsrn`

Computes the surface normal from derivative information.

#### `curvps`

Evaluates directional curvature on a surface using the first and second fundamental forms.

#### `xmaxpa`

Computes patch-corner knot coordinates and bounding-box-like min/max extents for a surface.

### Surface curve extraction and parameter utilities

#### `finded`

Extracts the polynomial curve corresponding to one surface edge.

#### `getiso`

Extracts the polynomial curve corresponding to an iso-parametric line on a surface.

#### `gtlcor`

Converts between global and local parametric coordinates for surfaces and curves.

#### `loc2gl`

Converts between global curve coordinates and local arc coordinates.

### Surface nearest-point and optimization routines

#### `distps`

Distance objective for line-restricted surface search. It evaluates `||R(x)-xp||` where the surface parameters are constrained to move along a line in the parametric plane.

#### `limiuv`

Computes allowable search bounds along a parametric search line so the search remains inside a rectangular parameter box.

#### `brauv`

Bracketing routine used by `xminli`. It finds an interval along a search direction that contains a local minimum of the surface-point distance function.

#### `xminli`

One-dimensional line minimizer on a surface. It minimizes distance along a line in the parametric plane and is used inside higher-level 2D searches.

#### `locsur`

Closest-point routine on an interpolated bicubic surface using Newton-Raphson plus steepest descent.

#### `locsu2`

Closest-point routine on an interpolated bicubic surface using a more direct Newton solve based on the local metric/Hessian approximation.

#### `locsug`

Closest-point routine on an interpolated bicubic surface using a Newton-like scheme with a conjugate-gradient-style search direction update.

#### `locps3`

Closest-point routine on a packed polynomial surface using a Powell-like direction-set method over the 2D parameter plane.

#### `locps4`

Higher-level closest-point routine on a packed polynomial surface. It is documented as a modified Powell method with repeated line searches based on `xminli`.

### Geometry database I/O and indexing

These routines assume a legacy direct-access record-based file layout used by the original software environment.

#### `readca`

Top-level reader for a CATIA neutral file. It parses geometry and stores it into the internal database structures.

#### `reseca`

Reads one curve entity from a CATIA neutral file and converts it into the internal curve database format.

#### `resuca`

Reads one surface entity from a CATIA neutral file and converts it into the internal surface database format.

#### `buicvf`

Writes a curve definition into the database file structure.

#### `buisgf`

Writes a surface definition into the database file structure.

#### `loas1d`

Loads a polynomial curve from the database into `isn` and `xln1d`.

#### `loas2d`

Loads a polynomial surface from the database into `keypa` and `aps`.

#### `builtg`

Scans the database key files and builds local-to-global numbering maps for curves and surfaces.

#### `ilocli`

Looks up the local index of a curve or surface from a local-to-global pointer array.

### Utility, compatibility, and low-level support

#### `testre`

Small helper that copies selected cross-derivative values from `tanret` into a temporary array. It appears to be a debugging or compatibility utility used during surface tangent assembly.

#### `cholen`, `coeff`, `trcoe`, `trif`, `tris`, `tris2`, `putimp`, `movt`

These are the main low-level assembly helpers behind the higher-level interpolation routines.

### `new/` variants

The `new/` directory contains alternate implementations or experiments around the main spline builders:

- `new/psplin`: alternate `psplin`.
- `new/pspli2`: alternate `pspli2`.
- `new/evtan`: alternate tangent solver.
- `new/surfb`: alternate surface builder.
- `new/uniform`: uniform parametrization helper.
- `new/centripet`: centripetal parametrization helper.
- `new/correctt`: imposed-tangent correction helper.

These should be treated as related but separate variants, not automatically interchangeable drop-ins without checking call signatures and output conventions.

## Object List In The Static Library

The current archive contents are defined by `myrules.make` and include:

- curve interpolation and evaluation: `psplin`, `pspli2`, `pspl`, `getp2`, `gettan`, `gett2`, `getpoi`, `getk`, `getk2`, `curva`, `slen`, `coeff`, `xlenst`, `lengt`, `patc1d`, `evps1d`
- surface interpolation and evaluation: `getta3`, `surfb`, `evapa2`, `evapa3`, `evapat`, `evaqfg`, `evppar`, `gpsur`, `evsurg`, `evsrn`, `gders1`, `fertob`, `curvps`
- nearest-point and search routines: `brasf`, `locpfs`, `xmin1d`, `fundis`, `locsur`, `locsu2`, `locsug`, `locps3`, `locps4`, `xminli`, `distps`, `brauv`, `limiuv`
- database and extraction helpers: `readca`, `reseca`, `resuca`, `buicvf`, `buisgf`, `loas1d`, `loas2d`, `builtg`, `ilocli`, `finded`, `getiso`, `xmaxpa`, `gtlcor`, `loc2gl`

## Practical Entry Points

If you need a starting point for using or modernizing this code, these are the most important routines to read first.

### For curve interpolation

- `psplin`
- `evtan`
- `getp2`
- `slen`
- `patc1d`

### For surface interpolation

- `getta3`
- `surfb`
- `evapa2`
- `gpsur`
- `gders1`

### For closest-point queries

- curves: `brasf`, `xmin1d`, `locpfs`
- surfaces from interpolation form: `locsur`, `locsu2`, `locsug`
- surfaces from polynomial patch form: `locps3`, `locps4`, `xminli`

### For database import/export

- `readca`
- `reseca`
- `resuca`
- `buicvf`
- `buisgf`
- `loas1d`
- `loas2d`

## Caveats

- This is legacy fixed-form Fortran and still emits several old-style loop warnings with modern compilers.
- Some routines are documented much better than others; a few small helpers had to be described from code inspection.
- The database routines assume direct-access files and the original record layout.
- The `new/` directory contains variants, not a clearly separated replacement API.
- The codebase predates modules and explicit interfaces, so callers must match argument conventions exactly.

## Suggested Modernization Path

If this library is going to be maintained further, the highest-value cleanup steps are:

1. add explicit interfaces or wrap the code in modern Fortran modules,
2. document the packed `aps`, `keypa`, `xln1d`, and `isn` layouts more formally,
3. separate database I/O from pure geometry kernels,
4. decide whether the `new/` routines replace or complement the top-level ones,
5. add regression tests for curve evaluation, surface evaluation, and closest-point searches.

