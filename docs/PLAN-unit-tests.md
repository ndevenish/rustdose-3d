# Plan: Port Java unit tests to Rust

Port the ~45 Java unit tests that have no Rust equivalent. Place them as
`#[cfg(test)] mod tests` inline in the source files (not in `tests/`), since
they exercise internal/concrete types rather than the public pipeline API.

## 1. DDM tests — `src/ddm.rs` (7 tests, trivial)

Port `DDMTests.java`:
- `DdmSimple::calc_decay()` always returns 1.0
- `DdmLinear::calc_decay()` in [0, 1] for positive dose
- `DdmLeal::calc_decay()` in [0, 1] for positive dose
- `DdmLinear` monotonically decreases with dose
- `DdmLeal` monotonically decreases with dose
- `DdmLinear` returns 1.0 at zero dose
- `DdmLeal` returns 1.0 at zero dose

## 2. Histogram tests — `src/output/exposure_summary.rs` (3 tests, easy)

Port `HistogramTest.java`. Rust histogram is built into `ExposureSummary`:
- `dose_hist_break()` returns correct bin boundaries
- Binning via `summary_observation()` places values in correct bins
- `dose_hist_normalised()` sums to 1.0

## 3. Container tests — `src/container.rs` (2 tests, easy)

Port `ContainerTests.java`:
- Transparent container: zero attenuation at any energy
- Elemental/mixture consistency for alanine (if both are implemented)

## 4. Wedge tests — `src/wedge.rs` (2 tests, easy)

Port `WedgeTest.java`:
- Default values for optional parameters (angular resolution, offsets)
- Required parameters (start/end angle, exposure time)

## 5. Experiment observer tests — `src/experiment.rs` (2 tests, easy)

Port `ExperimentTest.java`:
- Observer pattern with null crystal/beam/wedge doesn't panic
- Multiple observers receive events in correct order

## 6. Factory tests — `src/beam/mod.rs`, `src/crystal/mod.rs` (8 tests, easy)

Port `BeamFactoryTest.java` and `CrystalFactoryTest.java`:
- `create_beam()` returns Ok for TopHat, Gaussian configs
- `create_beam()` returns Err for invalid/empty type
- `create_crystal()` returns Ok for Cuboid, Spherical configs
- `create_crystal()` returns Err for invalid/empty type

Skip `ClassFactoryTest.java` (Java reflection-based, not applicable to Rust).

## 7. Beam intensity tests — `src/beam/tophat.rs`, `src/beam/gaussian.rs`, `src/beam/experimental.rs` (7 tests, moderate)

Port `BeamExperimentalTest.java`:
- Horizontal-only intensity profile
- Vertical-only intensity profile
- Combined horizontal + vertical profile
- Bilinear interpolation with known values

Add basic TopHat/Gaussian checks:
- TopHat returns uniform intensity inside aperture, zero outside
- Gaussian peak at centre, falls off with distance
- Circular collimation clips correctly

## 8. Crystal geometry tests — `src/crystal/cuboid.rs`, `src/crystal/polyhedron.rs` (6 tests, moderate)

Port `CrystalCuboidTest.java`:
- `find_depth()` at known voxel positions matches expected values
- 180° rotation symmetry: `find_depth(angle)` ≈ `find_depth(angle + π)`
- 360° rotation invariance for on-axis voxels

Port `CrystalPolyhedronTests.java`:
- Polyhedron depth matches cuboid depth for box-shaped OBJ
- Concave geometry depth (needs horseshoe OBJ fixture or inline mesh)

## 9. CoefCalc gap-fill — `src/coefcalc/compute.rs` (3 tests, moderate)

Port scenarios from `CoefCalcTests.java` not covered by existing
`tests/coefcalc_coefficients.rs`:
- Pure water H:O ratio check
- Heavy atom counting with monomer multiplication
- Cryo surrounding (water-based, oil-based) if container attenuation is wired

## Implementation order

1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 (easiest/fastest first)

## Notes

- Existing `tests/*.rs` files are integration tests (full pipeline vs golden
  files). Those stay as-is.
- Crystal geometry tests (§8) require constructing crystals from config, which
  needs a valid CoefCalc and unit cell. Use the insulin test fixture config
  as a template for minimal construction.
- BeamExperimental tests (§7) need small inline intensity grids (2×2, 3×3).
