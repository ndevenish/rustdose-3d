# Plan: Port RADDOSE-3D from Java to Rust

## Context

RADDOSE-3D is a ~40k LOC Java scientific application for radiation dose modelling in crystallography (MX, SAXS, XFEL, MicroED, small-molecule). The goal is a full-feature-parity **Rust** rewrite, providing:

- **Native CLI binary** — standalone, no JVM dependency
- **Library crate** — `use raddose3d::Experiment;` for Rust consumers
- **npm package via `wasm-pack`** — for web frontend / Node.js integration

Rust was chosen over TypeScript because: (1) IEEE 754 `f64` matches Java `double` semantics exactly, enabling bit-identical validation against Java output, (2) seeded RNG can replicate Java's `java.util.Random` for Monte Carlo reproducibility, (3) no WASM boundary complexity — all physics code is native, (4) single language to maintain.

---

## Project Structure

raddose3d/ should remain as it's own git repository, with a separate history.

```
raddose3d/
├── Cargo.toml                    # Workspace root
├── crates/
│   ├── raddose3d/                # Core library crate
│   │   ├── Cargo.toml
│   │   └── src/
│   │       ├── lib.rs            # Public API
│   │       ├── experiment.rs     # Orchestrator (observer pattern)
│   │       ├── crystal/
│   │       │   ├── mod.rs        # Crystal trait + factory
│   │       │   ├── cuboid.rs
│   │       │   ├── cylinder.rs
│   │       │   ├── spherical.rs
│   │       │   ├── polyhedron.rs # Base for cylinder/spherical/polyhedron
│   │       │   └── container.rs  # Transparent, Mixture, Elemental
│   │       ├── beam/
│   │       │   ├── mod.rs        # Beam trait + factory
│   │       │   ├── gaussian.rs
│   │       │   ├── tophat.rs
│   │       │   └── experimental.rs
│   │       ├── coefcalc/
│   │       │   ├── mod.rs        # CoefCalc trait
│   │       │   ├── compute.rs    # Core cross-section math (~5k lines Java)
│   │       │   ├── from_params.rs
│   │       │   ├── from_pdb.rs
│   │       │   ├── from_cif.rs
│   │       │   ├── from_sequence.rs
│   │       │   ├── saxs.rs
│   │       │   ├── small_molecules.rs
│   │       │   ├── micro_ed.rs
│   │       │   └── average.rs
│   │       ├── ddm/              # Dose decay models
│   │       │   ├── mod.rs        # DDM trait
│   │       │   ├── simple.rs
│   │       │   ├── linear.rs
│   │       │   ├── leal.rs
│   │       │   └── bfactor.rs
│   │       ├── wedge.rs
│   │       ├── output/
│   │       │   ├── mod.rs        # Output trait (observer)
│   │       │   ├── summary_text.rs
│   │       │   ├── summary_csv.rs
│   │       │   ├── dwds.rs
│   │       │   ├── voxel_dose.rs
│   │       │   ├── dose_state_csv.rs
│   │       │   ├── dose_state_r.rs
│   │       │   ├── rde_csv.rs
│   │       │   └── progress.rs
│   │       ├── element/
│   │       │   ├── mod.rs
│   │       │   ├── database.rs        # X-ray element DB
│   │       │   └── database_em.rs     # Electron element DB
│   │       ├── simulation/
│   │       │   ├── mc.rs              # Monte Carlo
│   │       │   ├── xfel.rs
│   │       │   └── micro_ed.rs
│   │       ├── energy_distribution.rs
│   │       └── constants/             # Embedded data files
│   │           ├── mod.rs             # include_str! / lazy_static parsing
│   │           ├── MuCalcConstants.txt
│   │           ├── ElCalcConstants.csv
│   │           ├── fullelsepa.csv
│   │           ├── Intensity_values.csv
│   │           └── ...
│   │
│   ├── raddose3d-parser/         # Input file parser (separate crate)
│   │   ├── Cargo.toml            # depends on: pest or nom
│   │   └── src/
│   │       ├── lib.rs            # parse(input: &str) → Config
│   │       ├── grammar.pest      # PEG grammar (pest) for RADDOSE-3D .txt format
│   │       └── config.rs         # CrystalConfig, BeamConfig, WedgeConfig types
│   │
│   ├── raddose3d-cli/            # CLI binary
│   │   ├── Cargo.toml            # depends on: clap, raddose3d, raddose3d-parser
│   │   └── src/
│   │       └── main.rs           # Argument parsing, wires up core + parser
│   │
│   └── raddose3d-wasm/           # WASM bindings (optional, for web/npm)
│       ├── Cargo.toml            # depends on: wasm-bindgen, raddose3d, raddose3d-parser
│       └── src/
│           └── lib.rs            # #[wasm_bindgen] run(input: &str) → JsValue
│
├── constants/                    # Source data files (copied into crate at build)
│   └── ... (from Java project)
│
└── tests/
    ├── fixtures/                 # Input files from Java examples/
    ├── golden/                   # Expected outputs captured from Java jar
    ├── parser.rs
    ├── crystal.rs
    ├── coefcalc.rs
    ├── integration.rs            # End-to-end: input → output diff vs golden
    └── scripts/
        └── generate_golden.sh    # Runs Java jar on all fixtures, captures output
```

---

## Key Rust Crate Dependencies

| Crate | Purpose | Replaces |
|-------|---------|----------|
| `ndarray` | N-dimensional arrays for voxel grids (`Array3<f64>`) | Java `double[][][]` |
| `nalgebra` | 3×3 rotation matrices, vector math | Inline Java matrix code |
| `statrs` | Normal distribution CDF/PDF, statistical functions | Apache Commons Math |
| `pest` | PEG parser generator (compile-time grammar) | ANTLR 3 |
| `clap` | CLI argument parsing | Manual Java arg parsing |
| `wasm-bindgen` | WASM bindings for npm/browser target | N/A |
| `serde` + `serde_json` | JSON input/output support | N/A (new feature) |
| `rayon` | Optional parallelism for voxel loops | N/A (Java is single-threaded) |

**No large framework dependencies.** Total added dependency weight is modest.

### Java RNG Compatibility

To reproduce Java's `java.util.Random` for Monte Carlo validation:
```rust
struct JavaRandom { seed: i64 }
impl JavaRandom {
    fn next_bits(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(0x5DEECE66D) + 0xB) & 0xFFFFFFFFFFFF;
        (self.seed >> (48 - bits)) as i32
    }
    fn next_double(&mut self) -> f64 {
        let hi = self.next_bits(26) as i64;
        let lo = self.next_bits(27) as i64;
        (((hi << 27) + lo) as f64) / ((1i64 << 53) as f64)
    }
}
```
This enables bit-identical MC output when given the same seed.

---

## Implementation Phases

Note that commits should be made frequently, at least every numbered
item inside a phase, but more granular than that for large tasks.

### Phase 1: Skeleton + Parser + Element Database

**Goal:** Parse input files and load element data. No simulation yet.

1. **Initialize Cargo workspace** with crate structure above
2. **Write pest grammar** for RADDOSE-3D input format
   - Source reference: `lib/antlrworks-parsergenerator/Inputfile.g`
   - Case-insensitive keywords, FLOAT/STRING/ELEMENT tokens, `#` `!` `//` comments
   - Output: `Config { crystals: Vec<CrystalConfig>, beams: Vec<BeamConfig>, wedges: Vec<WedgeConfig> }`
   - Test against all files in `examples/` and `input.txt`
3. **Port Element + ElementDatabase**
   - Source: `src/se/raddo/raddose3D/Element.java`, `ElementDatabase.java`
   - Embed `constants/MuCalcConstants.txt` via `include_str!`
   - Module-level `lazy_static` or `OnceLock` initialization
4. **Port ElementDatabaseEM**
   - Source: `ElementDatabaseEM.java`, `ElementEM.java`
   - Embed `constants/fullelsepa.csv`, `constants/ElCalcConstants.csv`
5. **Port constants readers** — CSV/TXT parsing with embedded data

**Verification:** Parse all example input files. Query cross-sections for several elements/energies, compare against Java output values.

---

### Phase 2: Core Types + Factories + Cuboid Simulation

**Goal:** Working dose calculation for a cuboid crystal.

1. **Define core traits**
   - `Crystal` trait (expose, find_depth, get_cryst_coord, etc.)
   - `Beam` trait (beam_intensity, etc.)
   - `CoefCalc` trait (calculate_coefficients)
   - `DDM` trait (calc_decay)
   - `Output` trait (observer: publish_crystal, publish_beam, publish_wedge)
   - `Container` trait
   - `Wedge` struct
2. **Port factories** — map config type strings → concrete types
3. **Port CoefCalcCompute + CoefCalcFromParams** (default RD3D mode)
   - Source: `CoefCalcCompute.java` (5,340 lines) — port incrementally:
     a. `calculate_coefficients_all()` — main entry
     b. Element cross-section lookups
     c. Polynomial coefficient expansion
     d. Fluorescence/Auger yield calculations
4. **Port CrystalCuboid**
   - Source: `CrystalCuboid.java`
   - Voxel grid via `ndarray::Array3<f64>`
   - Occupancy, `find_depth()`, `get_cryst_coord()`
5. **Port BeamGaussian + BeamTophat**
   - Source: `BeamGaussian.java`, `BeamTophat.java`
   - Use `statrs::distribution::Normal` for CDF
6. **Port Experiment orchestrator**
   - Source: `Experiment.java` (191 lines)
7. **Port OutputSummaryText + OutputSummaryCSV** (minimum for validation)
8. **Wire up full pipeline:** parse → construct → expose → output

**Verification:** Run insulin example from README. Compare Average DWD and Max Dose against Java output — target **bit-identical** results (not ±0.1%).

---

### Phase 3: Remaining Crystal Geometries

**Goal:** Port all crystal shapes.

1. **Port CrystalPolyhedron** (3,087 lines — most complex)
   - Vertex/face mesh geometry
   - Ray-triangle intersection for voxel occupancy and depth
   - PE/FL escape factor arrays
   - Cryo-surrounding voxel lattice
   - Use `nalgebra` for rotation matrices and vector operations
2. **Port CrystalCylinder** (extends Polyhedron)
   - 32-point polygonal base approximation
3. **Port CrystalSpherical / CrystalSphericalNew**
   - SphericalNew: icosahedron mesh (42 vertices, 80 faces)
4. **Port ImportWireframeObj** for .obj polyhedra
5. **Port Container implementations**

**Verification:** Port Java tests from `tests/CrystalPolyhedronTests.java`. Compare voxel coordinates and depth values.

---

### Phase 4: Remaining CoefCalc Modes + Beam Types

**Goal:** Full coefficient calculation and beam support.

1. **Port remaining CoefCalc implementations:**
   - `CoefCalcFromPDB` (679 lines) — PDB file parsing
   - `CoefCalcFromCIF` — CIF parsing
   - `CoefCalcFromSequence` — FASTA sequence → composition
   - `CoefCalcSAXS`, `CoefCalcFromSequenceSAXS`
   - `CoefCalcSmallMolecules`
   - `CoefCalcAverage`
   - `CoefCalcRaddose` (legacy V2 formula)
   - `CoefCalcMicroED`
2. **Port BeamExperimental** — 2D intensity grid, bilinear interpolation
3. **Port EnergyDistribution** — pink beam Gaussian energy sampling

**Verification:** Compare computed absorption coefficients against Java for test compositions across all CoefCalc modes.

---

### Phase 5: DDM + All Output Modules

**Goal:** All dose decay models and output formats.

1. **Port DDM implementations:** Simple, Linear, Leal (embedded BEST dataset), Bfactor
2. **Port all Output implementations:**
   - SummaryText, SummaryCSV, DWDs, VoxelDose, VoxelFluences
   - FinalDoseStateCSV, FinalDoseStateR, RDECSV
   - ProgressIndicator, ProgressEstimate
   - FluencePerDoseHistCSV
3. **Port Writer abstraction:** stdout, file, string buffer, tee

**Verification:** Full pipeline with all output types enabled. `diff` output files against Java golden files — target byte-identical.

---

### Phase 6: Monte Carlo, XFEL, MicroED

**Goal:** Port the three advanced simulation engines (~13.5k LOC).

1. **Port MC** (~4,500 lines)
   - Auger cascade (19 elements), ionization shells (K, L1-3, M1-5)
   - Energy straggling (100 bins), angular distribution (50 bins)
   - Implement `JavaRandom` for reproducible seeded output
   - Consider `rayon` parallelization (improvement over single-threaded Java)
2. **Port XFEL** (~4,500 lines)
   - Per-pulse dosimetry, cascade chains (shared code with MC)
3. **Port MicroED** (~4,500 lines)
   - Electron beam mode, elastic/inelastic scattering

**Verification:** With `JavaRandom` and matching seed, MC/XFEL/MicroED output should be bit-identical to Java. If any randomness source differs, fall back to statistical comparison (mean/std within 1%).

---

### Phase 7: CLI + WASM + Polish

**Goal:** Production-ready CLI binary and optional WASM package.

1. **CLI (`raddose3d-cli`):**
   - Match Java arguments: `-i`, `-o`, `-p`, `-r`, `-V`, `-?`, `-t`
   - `clap` derive-based arg parsing
   - STDIN support (`-i -`)
   - Progress output to stderr
   - Cross-compile: Linux, macOS, Windows binaries via `cross`
2. **WASM package (`raddose3d-wasm`):**
   - `wasm-bindgen` API: `run(input: &str) → String` (summary) + structured `JsValue` results
   - `wasm-pack build --target web` for browser
   - `wasm-pack build --target nodejs` for Node.js
   - Publish to npm as `@raddose3d/wasm`
3. **JSON input support** — `serde` deserialization as alternative to `.txt` format
4. **Library API:**
   ```rust
   use raddose3d::{Experiment, Config, parse_input};

   // High-level
   let config = parse_input(input_text)?;
   let results = raddose3d::run(config)?;
   println!("{}", results.average_dwd);

   // Low-level
   let mut exp = Experiment::new();
   exp.set_crystal(crystal);
   exp.set_beam(beam);
   exp.expose_wedge(wedge);
   ```

**Verification:** `raddose3d -i input.txt` produces byte-identical output to `java -jar raddose3d.jar -i input.txt` for all example files. WASM output matches native binary output.

---

## Testing Strategy

### Golden File Regression
```bash
# One-time: generate golden outputs from Java
for f in tests/fixtures/*.txt; do
  java -jar raddose3d.jar -i "$f" -p "tests/golden/$(basename $f .txt)-"
done

# CI: compare Rust output against golden
cargo test --test integration
```

### Per-Function Cross-Validation
For high-risk numerical code (CoefCalcCompute, Crystal.expose), write Java test harnesses that print intermediate values, then write Rust tests that assert identical values:
```rust
#[test]
fn coefcalc_insulin() {
    let coefs = CoefCalcFromParams::new(/* insulin params */).calculate();
    // Values extracted from Java debug output:
    assert_eq!(coefs.mu_abs, 3.4056789012345678);  // exact f64
    assert_eq!(coefs.mu_att, 4.1234567890123456);
}
```

### Continuous Integration
- `cargo test` — unit + integration tests
- `cargo clippy` — lint
- `cargo fmt --check` — formatting
- Golden file diff against Java jar (requires Java in CI)
- WASM build + test via `wasm-pack test --node`

---

## Risk Areas

1. **CoefCalcCompute (5,340 lines)** — Densest physics code. Mitigation: port function-by-function with per-function value comparison against Java.

2. **CrystalPolyhedron (3,087 lines)** — Complex 3D geometry. Mitigation: port Java tests first as Rust tests, implement until they pass.

3. **Floating-point operation ordering** — Java and Rust may compile `a + b + c` differently. Mitigation: use explicit parenthesization matching the Java source. Verify with exact `f64` bit comparisons early (Phase 2) to catch ordering issues before they compound.

4. **Monte Carlo reproducibility** — Requires matching Java's RNG exactly. Mitigation: implement `JavaRandom` wrapper, verify sequence output matches `java.util.Random` for 10k values before using in simulation.

---

## Estimated Scope

| Phase | Rust LOC (est.) | Java LOC covered | Key deliverable |
|-------|-----------------|------------------|-----------------|
| 1 | ~1,200 | ~2,500 | Parser + element DB |
| 2 | ~3,500 | ~8,000 | Working cuboid simulation |
| 3 | ~2,500 | ~5,000 | All crystal shapes |
| 4 | ~3,000 | ~9,000 | All CoefCalc + beams |
| 5 | ~1,500 | ~4,000 | DDM + outputs |
| 6 | ~4,500 | ~13,500 | MC/XFEL/MicroED |
| 7 | ~1,500 | — | CLI + WASM + polish |
| **Total** | **~17,700** | **~42,000** | Full parity |

Rust is more compact than Java (~42% ratio) due to pattern matching, iterators, trait dispatch, and no getter/setter boilerplate. The physics math will be similar density.

---

## Bonus: Improvements Over Java

Things that come free or cheap with the Rust architecture:

- **Parallelism via `rayon`** — the voxel loop and MC simulations can trivially parallelize with `par_iter()`. Java version is single-threaded.
- **Memory safety** — no array bounds errors, null pointer exceptions, or heap space crashes that Java users hit with large crystals.
- **Static binary** — no JVM installation, no `java -jar`, no classpath issues. Just `raddose3d -i input.txt`.
- **WASM for browser** — compute dose directly in the browser, replacing the PHP server-side architecture entirely.
- **`serde` JSON I/O** — structured input/output for programmatic use, free with derive macros.
