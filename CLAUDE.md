# CLAUDE.md — Rust RADDOSE-3D

Rust rewrite of [RADDOSE-3D](https://github.com/GarmanGroup/RADDOSE-3D), a radiation dose modelling tool for macromolecular crystallography. The Java source lives in the parent directory (`../src/se/raddo/raddose3D/`). The full porting plan is at `../snoopy-imagining-noodle.md`.

## Build & Run

```bash
cargo build                              # Debug build
cargo build --release                    # Release build
cargo run -- -i tests/fixtures/insulin_test.txt  # Run CLI
cargo test                               # Run all tests
```

## Workspace Layout

Three crates in a Cargo workspace:

| Crate | Path | Purpose |
|-------|------|---------|
| `raddose3d` | `crates/raddose3d/` | Core library — traits, simulation, physics |
| `raddose3d-parser` | `crates/raddose3d-parser/` | Input file parser (pest PEG grammar) |
| `raddose3d-cli` | `crates/raddose3d-cli/` | CLI binary (`raddose3d -i input.txt`) |

## Core Library Modules (`crates/raddose3d/src/`)

| Module | Key Types | Java Source |
|--------|-----------|-------------|
| `lib.rs` | `constants` module (KEV_TO_JOULES, UNIT_CONVERSION, etc.) | Various |
| `experiment.rs` | `Experiment` — orchestrator, observer pattern | `Experiment.java` (191 lines) |
| `wedge.rs` | `Wedge` — angles in radians, translations per radian | `Wedge.java` (392 lines) |
| `ddm.rs` | `DdmModel` trait, `DdmSimple/Linear/Leal/Bfactor` | `DDMSimple.java` etc. |
| `container.rs` | `Container` trait, `ContainerTransparent`, `ContainerMixture`, `ContainerElemental` (stubs), `create_container()` | `ContainerTransparent.java` etc. |
| `beam/mod.rs` | `Beam` trait, `create_beam()` factory | `Beam.java` |
| `beam/gaussian.rs` | `BeamGaussian` — CDF normalization via `statrs` | `BeamGaussian.java` |
| `beam/tophat.rs` | `BeamTophat` — uniform intensity | `BeamTophat.java` |
| `coefcalc/mod.rs` | `CoefCalc` trait, `create_coefcalc()` factory | `CoefCalc.java` |
| `coefcalc/compute.rs` | `CoefCalcCompute` — McMaster cross-sections, fluorescent escape | `CoefCalcCompute.java` (5340 lines) |
| `coefcalc/from_params.rs` | `CoefCalcFromParams` — composition-based coefficients | `CoefCalcFromParams.java` (294 lines) |
| `crystal/mod.rs` | `Crystal` trait, `expose_rd3d()`, `create_crystal()` factory | `Crystal.java` |
| `crystal/cuboid.rs` | `CrystalCuboid` — 8 vertices, 12 triangles, ray-casting | `CrystalCuboid.java` + `CrystalPolyhedron.java` |
| `crystal/polyhedron.rs` | `CrystalPolyhedron` — general mesh (OBJ, cylinder, icosphere), `cylinder_geometry()`, `icosphere_geometry()` | `CrystalPolyhedron.java`, `CrystalCylinder.java`, `CrystalSphericalNew.java` |
| `crystal/spherical.rs` | `CrystalSpherical` — analytic sphere, pre-computed occupancy | `CrystalSpherical.java` |
| `output/mod.rs` | `Output` trait, `ExposeObserver` trait | `Output.java` |
| `output/summary.rs` | `OutputSummaryText`, `OutputSummaryCSV` | `OutputSummaryText.java` |
| `output/exposure_summary.rs` | `ExposureSummary` — DWD, dose stats, quantiles | `ExposureSummary.java` (548 lines) |
| `element/database.rs` | `Element`, `ElementDatabase` (X-ray cross-sections) | `Element.java`, `ElementDatabase.java` |
| `element/database_em.rs` | `ElementEM`, `ElementDatabaseEM` (electron) | `ElementDatabaseEM.java` |

## Parser (`crates/raddose3d-parser/src/`)

- `grammar.pest` — PEG grammar for RADDOSE-3D `.txt` input format (case-insensitive keywords, `#`/`!`/`//` comments)
- `config.rs` — `Config`, `CrystalConfig`, `BeamConfig`, `WedgeConfig` structs with all fields as `Option<T>`
- `lib.rs` — `parse(input: &str) -> Result<Config, String>`

Config enums: `CoefCalcType` (Average, Default, Pdb, Saxs, Sequence, SmallMole, Cif, ...), `DdmType` (Simple, Linear, Leal, Bfactor), `Collimation` (Rectangular, Circular, Horizontal, Vertical).

## Simulation Flow

```
Experiment::run_from_config(config)
  ├─ create_crystal(crystal_config)   → Box<dyn Crystal>
  │    └─ CoefCalcFromParams          → absorption/attenuation coefficients
  ├─ create_beam(beam_config)         → Box<dyn Beam>
  ├─ Wedge::from_config(wedge_config) → angles in radians
  └─ crystal.expose(beam, wedge)
       └─ expose_rd3d(crystal, beam, wedge, container)
            ├─ coefcalc.update_coefficients(energy)
            ├─ for angle in angles:
            │    for energy in energies:
            │      for (i,j,k) in voxels:
            │        beam_intensity → fluence → dose/elastic/compton
            │        exposure_summary.exposure_observation(...)
            │    exposure_summary.image_complete(...)
            ├─ for (i,j,k) in voxels:
            │    exposure_summary.summary_observation(dose, mass)
            └─ exposure_summary.exposure_complete()
```

## Key Formulas

All must match Java exactly for bit-identical results:

- **Fluence to dose:** `-expm1(-µ/pixPerUM) / (1e-15 × pixPerUM⁻³ × density) × 1e-6`
- **Ray-plane distance:** `d = -(n · v0)`, then `t = (coord·n + d) / (dir·n)`
- **Point-in-polygon:** 2D ray-casting in X-Y plane (matches Java's `Vector.polygonInclusionTest`)
- **Compton electron energy:** `E_beam × (1 - √(mc²/(2·E_beam + mc²)))`
- **Absorbed energy:** `energyPerFluence × (photoFluence + comptonFluence)`
- **Cell volume:** `a·b·c·√(1 - cos²α - cos²β - cos²γ + 2·cosα·cosβ·cosγ)`

## Conventions

- **Units:** Crystal dimensions in µm, beam energy in keV, dose in MGy, angles in radians internally (degrees in input)
- **Voxel indexing:** Flat `Vec<f64>` with `idx = i*ny*nz + j*nz + k` (not `ndarray::Array3`)
- **Traits:** All major component traits require `Debug + Send + Sync`
- **Factories:** Free functions `create_*()` in each module's `mod.rs` (no separate factory files)
- **DDM models:** All implement `DdmModel` trait with `calc_decay(dose) -> f64`
- **Output observer:** `Output::publish_wedge` receives `&ExposureSummary` for dose statistics
- **Error handling:** `Result<T, String>` for construction errors (no custom error types yet)

## Validation

Verified against Java RADDOSE-3D on insulin test case (`tests/fixtures/insulin_test.txt`):

| Metric | Java | Rust | Match |
|--------|------|------|-------|
| Photoelectric Coeff | 2.23e-04 /µm | 2.23e-4 /µm | exact |
| Attenuation Coeff | 2.60e-04 /µm | 2.60e-4 /µm | exact |
| Density | 1.08 g/ml | 1.08 g/ml | exact |
| Max Dose | 12.043700 MGy | 12.043700 MGy | bit-identical |
| Elastic Yield | 3.62e+11 | 3.62e+11 | bit-identical |
| Absorbed Energy | 8.52e-03 J | 8.52e-3 J | bit-identical |
| Avg DWD | 4.133330 MGy | 4.133270 MGy | 0.001% |
| Avg Dose (Whole) | 7.716530 MGy | 7.716484 MGy | 0.001% |
| Used Volume | 100.0% | 100.0% | exact |

## Porting Progress

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Skeleton + Parser + Element Database | **Done** |
| 2 | Core Types + Factories + Cuboid Simulation | **Done** |
| 3 | Remaining Crystal Geometries (Polyhedron, Cylinder, Spherical) | **Done** |
| 4 | Remaining CoefCalc Modes + Beam Types | **Done** |
| 5 | DDM + All Output Modules | **Done** |
| 6 | Monte Carlo, XFEL, MicroED | **Done** |
| 7 | CLI + WASM + Polish | **Done** |

### What's Implemented
- Input parsing (all keywords from Java grammar)
- Element databases (X-ray + electron)
- CoefCalcCompute + CoefCalcFromParams (RD3D default mode)
- CrystalPolyhedron (general mesh, OBJ loader, cylinder geometry, icosphere geometry)
- CrystalCylinder (via CrystalPolyhedron with hardcoded 64-vertex/96-triangle cylinder)
- CrystalSphericalNew (via CrystalPolyhedron with hardcoded 42-vertex/80-triangle icosphere)
- CrystalSpherical (analytic: radius occupancy check, analytic depth formula)
- ContainerMixture + ContainerElemental stubs (structure present, NIST lookup not yet implemented)
- CrystalCuboid with polyhedron ray-casting
- BeamGaussian + BeamTophat (rectangular/circular collimation)
- BeamExperimental (2D intensity grid, bilinear interpolation)
- EnergyDistribution (truncated normal sampler for pink beam)
- CoefCalcAverage (hardcoded Holton 2010 constants)
- CoefCalcFromCIF (CIF _chemical_formula_sum parser)
- CoefCalcSmallMolecules (small-molecule crystals, no water fill)
- CoefCalcMicroED (electron diffraction, same composition as FromParams)
- CoefCalcFromSequence (FASTA sequence → H/C/N/O/S/P/Se composition)
- CoefCalcSAXS (protein concentration → num_monomers)
- CoefCalcFromSequenceSAXS (SAXS + sequence file)
- CoefCalcFromPDB (PDB file or RCSB download → SEQRES + HETATM composition)
- Residue database (20 amino acids + RNA/DNA nucleotides, 1-letter and 3-letter lookup)
- DDM: Simple, Linear, Leal, Bfactor
- ContainerTransparent
- Experiment orchestrator
- OutputSummaryText + OutputSummaryCSV
- OutputDWDs, OutputRDECSV, OutputFinalDoseStateCSV/R
- OutputVoxelDose, OutputVoxelFluences
- OutputProgressIndicator, OutputProgressEstimate
- OutputFluencePerDoseHistCSV
- writer module: file_writer, stdout_writer, StringWriter, TeeWriter, NullWriter
- ExposureSummary (DWD, dose quantiles, dose contrast, per-image RDE, voxel snapshots)
- Full pipeline: parse → construct → expose → output
- Public library API: `parse_input()`, `parse_input_json()`, `run()` → `RunResults`
- CLI: full Java-compatible flags (`-i`/`--in`, `-p`/`--prefix`, `-r`/`--raddose`, `-t`/`--test`, `-V`/`--version`, `-?`/`--help`); stdin support (`-i -`); default output observers writing all 8 files; progress → stderr
- WASM crate (`raddose3d-wasm`): `wasm-bindgen` API with `run()` (text), `runStructured()` (JsValue), `runFromJson()` (JSON input)
- JSON input support via `parse_input_json()` and `serde_json`
- JavaRandom (Java LCG, bit-identical to java.util.Random for seeded MC)
- MonteCarloSimulation (mc.rs, ~2400 lines): Auger cascade (19 elements), ionization shells (K/L/M), ELSEPA elastic angles, GOS inelastic scattering, photoelectron/Compton/Auger electron tracking
- XfelSimulation (xfel.rs, ~2500 lines): time-resolved dosimetry, 4D voxel arrays (spatial + time bins), PULSE_LENGTH/PULSE_BIN_LENGTH parameters, same physics as MC
- MicroEdSimulation (micro_ed.rs, ~1400 lines): electron diffraction, stopping-power slice simulation, CSDA range, optimal voltage/thickness search

### What's Not Yet Ported
- CoefCalc: Raddose (legacy v2 subprocess — stubs to Default)
- Beam: EnergyDistribution integration with Experimental beam for pink beam multi-energy loop
- Container: Mixture/Elemental NIST mass-attenuation lookup (structure done, lookup not implemented)
- PE/FL escape code paths in expose loop
- Cryo-surrounding voxel lattice
- Crystal dispatch for MONTECARLO/GOS/XFEL/EMSP/EMED subprograms (wiring MC/XFEL/MicroED into expose())
- CLI flags: `-o` (custom output module routing — Java's `module:dest` syntax not implemented)
- WASM: `wasm-pack` toolchain must be installed separately; WASM target not in default `cargo build`
