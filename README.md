# Rustdose-3D

> **⚠️ Experimental / vibe-coded port — not for production use**
>
> This is an AI-assisted Rust rewrite of [RADDOSE-3D](https://github.com/GarmanGroup/RADDOSE-3D), the radiation dose modelling tool from the Garman Lab at the University of Oxford. It was written almost entirely by Claude (Opus 4 / Sonnet 4) with light human steering — think of it as a high-confidence draft rather than a reviewed implementation. Core outputs match Java on standard test cases to <0.1%, but several code paths are unvalidated, a few Java bugs are deliberately reproduced for compatibility, and a handful of features are silently stubbed. **If your results matter, cross-check against the canonical Java jar.**

Rustdose-3D models the time- and space-resolved radiation dose distribution in crystal volumes during macromolecular crystallography (MX), SAXS, electron diffraction (MicroED), and XFEL experiments. It is a complete reimplementation of [RADDOSE-3D v5](https://github.com/GarmanGroup/RADDOSE-3D) in Rust, accepting the same input file format and producing the same output files.

---

## Why this exists

- **Single binary, zero runtime deps** — no JVM, no `constants/` directory to ship alongside the jar, no heap size tuning.
- **Embeddable as a library** — `raddose3d` is a plain Rust crate; call it from anything.
- **WASM target** — runs in a browser or Node.js via `raddose3d-wasm` without modification.
- **Faster iteration** — Rust's type system caught several latent bugs during porting.

---

## Usage

The CLI accepts the same input files as the Java version:

```
./raddose3d -i path/to/input.txt
```

#### Minimal example (insulin MX experiment)

```
Crystal
Type Cuboid
Dimensions 100 100 100  # µm
PixelsPerMicron 0.5
AbsCoefCalc RD3D
UnitCell 78.02 78.02 78.02
NumMonomers 24
NumResidues 51
ProteinHeavyAtoms Zn 0.333 S 6
SolventHeavyConc P 425
SolventFraction 0.64

Beam
Type Gaussian
Flux 2e12
FWHM 20 70
Energy 12.1
Collimation Rectangular 100 100

Wedge 0 90
ExposureTime 50
```

Save as `MyInput.txt`, then:

```bash
./raddose3d -i MyInput.txt -p output_
```

Output files (`output_Summary.txt`, `output_Summary.csv`, `output_DWDs.csv`, etc.) are written to the current directory. Passing `-p -` writes the summary to stdout.

#### CLI flags

| Flag | Description |
|------|-------------|
| `-i <file>` / `--in <file>` | Input file (`-` for stdin) |
| `-p <prefix>` / `--prefix <prefix>` | Output file prefix (default: empty) |
| `-V` / `--version` | Print version |
| `-?` / `--help` | Print help |
| `-r <path>` | Path to legacy RADDOSE v2 binary (unused; accepted for compatibility) |

The `-o` flag (Java's custom output routing syntax `module:dest`) is **not implemented**.

---

## Build

Requires Rust 1.75+ and Cargo.

```bash
cargo build --release          # produces target/release/raddose3d
cargo test --release           # run the test suite
```

#### WASM

Requires `wasm-pack` and the `wasm32-unknown-unknown` rustup target:

```bash
rustup target add wasm32-unknown-unknown
curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
cd crates/raddose3d-wasm && wasm-pack build --target web
```

---

## Workspace layout

| Crate | Purpose |
|-------|---------|
| `crates/raddose3d` | Core library — physics, simulation, traits, output |
| `crates/raddose3d-parser` | PEG grammar (`pest`) for the `.txt` input format |
| `crates/raddose3d-cli` | CLI binary |
| `crates/raddose3d-wasm` | WASM bindings via `wasm-bindgen` |

---

## Validation status

Key results are validated against Java output on standard fixtures (bit-identical or within stated tolerance):

| Test case | Max Dose | Avg DWD | Notes |
|-----------|----------|---------|-------|
| Insulin (MX, Cuboid, RD3D) | exact | <0.001% | Golden |
| Polyhedron crystal (OBJ mesh) | <0.3% | <0.1% | Boundary-voxel FP variance |
| SAXS (CoefCalcSAXS) | <0.1% | <0.1% | Golden |
| PDB-derived composition | <0.1% | <0.1% | Golden |
| Multi-beam / multi-wedge | exact | <0.1% | Golden (regression for beam-ordering bug) |
| PE/FL escape | ~2% | ~1% | Random-sampling noise; see below |
| Monte Carlo | N/A | ~5% | Stochastic; seed-matched to Java LCG |
| XFEL | <0.1% | <0.1% | Time-resolved 4D dose |
| MicroED | <0.1% | <0.1% | Electron CSDA stopping power |

---

## Known divergences and incomplete implementation

This section documents where Rustdose-3D deliberately or unavoidably differs from Java RADDOSE-3D, and what is not yet ported.

### Deliberately reproduced Java bugs

These bugs exist in Java and are **intentionally replicated** so that outputs match:

**`findDepth` deduplication is a no-op** (`cuboid.rs`, `polyhedron.rs`)
Java's deduplication loop uses `==` on boxed `Double` objects, which is reference equality — always `false` for distinct objects. When a ray hits a shared triangle edge both adjacent triangles register the intersection, the even count trips the sanity-check fallback, and the voxel gets `depth=0`. Rust replicates this by not deduplicating either. Affects voxels on rotated crystal face boundaries; contributes ~0.3% Max Dose variance on coarse grids (`PixelsPerMicron < 0.5`).

**Gaussian beam circular-collimation normalization** (`beam/gaussian.rs`)
Java's `bivariateGaussianVolume` uses a polar-coordinate numerical integration (100 trapezoidal angular steps, analytical radial integral) rather than the more accurate Cartesian approach. Rust matches this algorithm. A Cartesian implementation produces ~0.047% normalization difference which is amplified to ~2.75% in DWD by steep DDM models (Leal, Bfactor).

**PE/FL escape: two silent bugs** (`crystal/escape.rs`)
- Bug 1: hardcoded `muabsIndex = 4` in the fluorescent escape distribution (should look up the correct index for each element).
- Bug 2: `energyPerFluence` from the cryo block leaks into the non-cryo path, affecting absorbed energy accounting.

Both are reproduced so that absorbed energy and FL escape outputs match Java exactly.

### Numerical divergences (not bugs, just FP path differences)

**PE escape DWD (~0.5–1% offset)**
Java recomputes photoelectron tracks per-angle with goniometer rotation applied; Rust computes them once. Additionally, float accumulation at the phi loop boundary causes a one-track count difference in `find_voxels_reached_by_pe`. These combine for a systematic ~0.5–1% DWD offset when PE escape is enabled.

**PE/FL escape Max Dose (±10%)**
Only one random PE track is sampled per voxel (`PE_ANGLE_RESOLUTION=1`, matching Java). This produces high per-run variance in both implementations. Both Java and Rust are in the same ballpark but individual runs will differ.

**Coarse-grid boundary voxels (~0.3% Max Dose)**
The exact set of voxels that receive `depth=0` due to the broken `findDepth` deduplication depends on whether near-zero intersection distances (~1e-14) round positive or negative, which differs between Java and Rust floating-point paths. Cannot be eliminated without matching every FP operation exactly.

### Incomplete / stubbed features

**Container: Mixture and Elemental** (`container.rs`)
`ContainerMixture` and `ContainerElemental` parse correctly but do not perform NIST mass-attenuation lookup for the surrounding material. They fall through to transparent behaviour. `ContainerTransparent` works fully.

**Pink beam (EnergyDistribution + Experimental beam)**
`EnergyDistribution` (truncated normal energy sampler) is implemented and parses, but the integration with `BeamExperimental` for the per-energy multi-loop is not wired up. Pink beam inputs will run with the central energy only.

**Legacy RADDOSE v2 subprocess (`CoefCalc RDV2` / `RDV3`)**
Java can shell out to the old RADDOSE v2 binary to compute absorption coefficients. This is silently stubbed to `Default` / `RD3D` behaviour. The `-r` CLI flag is accepted but ignored.

**Custom output routing (`-o module:dest`)**
Java's `-o` flag allows routing specific output modules to specific files. Not implemented; the eight default output files are always written to `{prefix}{Name}.{ext}`.

**Progress bar architecture**
The progress indicator (`[ 0%....100% ]`) is printed directly inside the simulation loop via `print!()` rather than being routed through the `OutputProgressIndicator` observer. This is a layering violation — the observer exists and is correct but is not called. Embedders and WASM consumers get the print side-effect silently (no output in WASM).

**Leal DDM: no golden validation**
`DdmLeal` implements the Leal et al. (2012) model correctly but no golden test compares it against Java output at non-trivial doses. Unit tests only check monotonicity and boundary properties.

---

## Interpreting the output

See the [RADDOSE-3D user guide](https://github.com/GarmanGroup/RADDOSE-3D/blob/master/doc/user-guide.pdf) — the output format is identical. The most informative quantity is the **Average Diffraction Weighted Dose (DWD)**, not Max Dose. Max Dose is the value comparable to legacy RADDOSE v2 output.

---

## Input format notes (same as Java)

- Comments: `#`, `!`, or `//`
- Keywords are case-insensitive
- Multiple `Beam` / `Wedge` pairs are supported; each `Beam` applies to the `Wedge` blocks that follow it until the next `Beam`
- `AngularResolution` defaults to `1.0°`; set it smaller for wedges under ~20°
- `PixelsPerMicron` defaults to `0.5`; increase for small crystals (<20 µm), decrease for large samples to manage memory
- Flux is the **post-collimation flux at the sample**, not the source flux

---

## Relationship to Java RADDOSE-3D

This port shares no code with the Java implementation. It was produced by reading the Java source and reimplementing the same physics and algorithms in Rust, with the goal of bit-identical or near-identical outputs. Where Java behaviour is clearly a bug (e.g. the `findDepth` deduplication) the bug is reproduced anyway to preserve output compatibility — the intent is a drop-in replacement, not a corrected simulation.

For the canonical, peer-reviewed implementation please use the [Java version](https://github.com/GarmanGroup/RADDOSE-3D). For questions about RADDOSE-3D science, the [original papers](#citations) and the Garman Lab are the authoritative sources.

---

## Citations

If you use results from this tool in published work, please cite the original RADDOSE-3D papers — this port is a reimplementation, not new science:

- Zeldin, O. B., Gerstel, M., & Garman, E. F. (2013). RADDOSE-3D: time- and space-resolved modelling of dose in macromolecular crystallography. *J. Appl. Cryst.*, 46, 1225–1230. https://doi.org/10.1107/S0021889813011461

- Dickerson, J. L., McCubbin P. T. N, Brooks-Bartlett J. C., Garman E. F. (2024). Doses for X-ray and electron diffraction: New features in RADDOSE-3D including intensity decay models. *Protein Science*, 33(7), e5005. https://doi.org/10.1002/pro.5005

- Bury, C. S., Brooks‐Bartlett, J. C., Walsh, S. P. and Garman, E. F. (2018). Estimate your dose: RADDOSE‐3D. *Protein Science*, 27, 217–228. https://doi.org/10.1002/pro.3302

Additional citations for SAXS, XFEL, Monte Carlo, small molecules, and fluorescent escape can be found in the [Java README](https://github.com/GarmanGroup/RADDOSE-3D/blob/master/README.md).
