# CLAUDE.md — Rust RADDOSE-3D

Rust rewrite of [RADDOSE-3D](https://github.com/GarmanGroup/RADDOSE-3D), a radiation dose modelling tool for macromolecular crystallography. The Java source lives in the parent directory (`../src/se/raddo/raddose3D/`). The full porting plan is at `../snoopy-imagining-noodle.md`.

## Build & Run

```bash
cargo build                              # Debug build
cargo build --release                    # Release build
cargo run -- -i tests/fixtures/insulin_test.txt  # Run CLI
cargo test                               # Run all tests
```

WASM is a separate build step and intentionally not part of `cargo build` — it requires `wasm-pack` (which bundles `wasm-opt`) and the `wasm32-unknown-unknown` rustup target:

```bash
rustup target add wasm32-unknown-unknown
curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
cd crates/raddose3d-wasm && wasm-pack build --target web
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

### TODO

#### Leal DDM golden test

The Leal DDM (`DdmLeal` in `ddm.rs`) was rewritten to implement the correct Leal et al. (2012) model with BEST_DATA + gamma/b0/beta parameters. There is currently no golden validation test covering this model — the existing unit tests only check generic properties (range, monotonicity, zero-dose = 1) rather than bit-identical Java output. A proper golden test needs a fixture where:
- DDM Leal is used with non-trivial DecayParam values
- The doses are low enough that `calcDecay` is not numerically zero (i.e. avoid extreme doses where floating-point underflow makes both Java and Rust DWD = 0)
- The Java golden Summary.txt output is captured and committed

#### Progress bar refactor + WASM progress callbacks

The progress bar is hardcoded in `expose_rd3d()` (`crystal/mod.rs`) via `print!()` to stdout. This is the wrong layer — progress is a UI concern and belongs outside the core simulation. Two specific problems:

1. **`OutputProgressIndicator` / `ExposeObserver` are not wired up.** The `ExposeObserver` trait (`output/mod.rs`) and `OutputProgressIndicator` (`output/progress.rs`) exist and are correct, but `expose_rd3d` never calls into any `ExposeObserver`. The trait implementations are dead code for their primary purpose; actual progress comes from the hardcoded `print!()`.

2. **Core library has a hidden stdout side-effect.** Any embedder (including WASM) gets this silently. In WASM the `print!()` calls compile but produce no output.

**What needs doing:**
- Remove the hardcoded `print!()` progress block from `expose_rd3d` and either (a) wire `ExposeObserver` into the simulation loop properly, or (b) accept an optional `Box<dyn FnMut(usize, usize)>` callback (images_complete, images_total).
- Add `runWithProgress(input: &str, on_progress: js_sys::Function)` to `raddose3d-wasm`, calling the JS function with `(imagesComplete, imagesTotal)` after each angle step.
- Note: since WASM runs on the main JS thread and blocks it, a live browser progress bar requires running the simulation in a Web Worker and relaying progress via `postMessage`. The callback is still useful for Node.js consumers and as groundwork for the worker approach.

### What's Not Yet Ported
- CoefCalc: Raddose (legacy v2 subprocess — stubs to Default)
- Beam: EnergyDistribution integration with Experimental beam for pink beam multi-energy loop
- Container: Mixture/Elemental NIST mass-attenuation lookup (structure done, lookup not implemented)
- Crystal dispatch for MONTECARLO/GOS/XFEL/EMSP/EMED subprograms (wiring MC/XFEL/MicroED into expose())
- CLI flags: `-o` (custom output module routing — Java's `module:dest` syntax not implemented)

### PE/FL Escape: Validation Status
PE/FL escape and cryo surrounding are implemented (`crystal/escape.rs`, integrated into `expose_rd3d()`). Validated on LiFePO₄ test case (1µm sphere, SMALLMOLE, 1.487 keV). Two Java bugs are intentionally reproduced for compatibility — see `docs/java-bugs-analysis.md`.

**Bug compatibility (implemented):**
- Bug 1: Hardcoded `muabsIndex = 4` in FL escape distribution (FL Escape now matches Java exactly: 8.23e0 J)
- Bug 2: `energyPerFluence` scope leak from cryo block (Absorbed Energy now matches Java exactly: 4.05e-9 J)

**Remaining variance:**
- PE Escape: matches within ~2% (random sampling noise at `PE_ANGLE_RESOLUTION=1`)
- DWD: ~0.5-1% systematic offset — from single-track count difference in `find_voxels_reached_by_pe` (float accumulation at phi loop boundary) and Java recomputing PE tracks per-angle with rotation while Rust computes once
- Max Dose: high variance in both Java and Rust (±10%) — single random track per voxel

**Notes:**
- `PE_ANGLE_RESOLUTION` (escape.rs line 21) controls track samples per voxel. Java hardcodes 1. Increasing reduces variance but costs O(n²) time.
- Gumbel PDF negative beta: At high density + low PE energy (e.g. ρ=3.6 g/mL, E_pe<0.82 keV), the Gumbel scale parameter goes negative. Java relies on sign cancellation. Rust matches this. See `gumbel_pdf()` in `escape.rs`.
- `Type Spherical` dispatch: maps to `CrystalSphericalNew` (icosphere mesh) matching Java.
