# Plan: Drastically Improve the RADDOSE-3D Differential Fuzzer

## Context

The fuzzer in `fuzz/` is a differential testing harness comparing Java vs Rust RADDOSE-3D implementations. It has two generation strategies (grammar-based 50%, mutation-based 50%) but both are conservative: narrow parameter ranges, heavy bias toward common code paths, no boundary/edge-value injection, no structural mutation, and the mutation path has no seeds to work from. The result is good coverage of the "happy path" but poor exploration of edge cases, rare code paths, and boundary conditions where bugs actually hide.

## Changes

### 1. Boundary Value Injection in `generate.py`

Add a `_boundary_perturb(cfg)` post-processing pass that runs with ~30% probability on grammar-generated configs and replaces sampled values with edge-case values:

- **Near-zero dimensions**: 0.1, 0.5, 1.0 um (tiny crystals)
- **Very large dimensions**: 1000, 5000 um (triggers voxel cap)
- **Extreme energy**: 1.0 keV, 50.0 keV (edge of tables)
- **Extreme flux**: 1e6, 1e16
- **PixelsPerMicron extremes**: 0.001 (very coarse), 10.0 (very fine, hits voxel cap)
- **Zero-angle wedge**: start=0, end=0 (single image)
- **Full rotation**: 0-360 with fine resolution
- **Solvent fraction**: 0.01, 0.99 (near-empty/near-full)
- **Unit cell**: 1.0 A (tiny), 500.0 A (huge)
- **num_residues**: 1, 2000
- **num_monomers**: 1, 200
- **Heavy atom counts**: 0.001 (trace), 1000 (massive)
- **DDM parameters**: near-zero (0.001), large (100.0)
- **MC sim_electrons**: minimum (100), very high (capped by budget)
- **Exposure time**: 1e-8 (femtosecond), 1000.0 (long)
- **FWHM**: 0.1 (pencil beam), 10000 (flood)

Implementation: a dict mapping config field names to lists of boundary values. `_boundary_perturb` picks 1-3 fields at random and replaces their values. Re-check cost after perturbation; if over budget, retry with different fields.

**File**: `generate.py`, new method `GrammarGenerator._boundary_perturb()` (~60 lines), called from `generate()`.

### 2. Structural Mutation in `generate.py`

Add a new generation strategy `"structural"` (and include it in `"both"`) that takes a grammar-generated config and applies structural changes:

- **Add extra wedge** to a segment (tests multi-wedge handling)
- **Add extra segment** (new beam + wedge block)
- **Remove optional fields** (heavy atoms, solvent conc, DDM, RNA/DNA)
- **Duplicate a segment** (identical beam, different angles)
- **Swap crystal type** while keeping other params (tests dimension reinterpretation)
- **Add container block** (currently not fuzzed at all) — `Container MIXTURE Ethanol 0.5` etc.
- **Multi-wedge on restricted modes**: add a second wedge to XFEL/MONTECARLO configs (tests error handling or multi-wedge support)

Implementation: `_structural_mutate(cfg)` takes a Config, applies 1-2 random structural transforms, returns modified Config. Called as a third strategy option.

**File**: `generate.py`, new method `GrammarGenerator._structural_mutate()` (~80 lines). Update `_generate()` in `fuzz.py` to include structural as a strategy option.

### 3. Widen Parameter Ranges in `generate.py`

Current ranges are conservative. Widen them:

| Parameter | Current | New |
|-----------|---------|-----|
| Cuboid dims | [5, 300] | [0.5, 2000] |
| Cylinder height | [200, 3000] | [10, 10000] |
| Cylinder diameter | [100, 2000] | [5, 5000] |
| Spherical dims | [5, 300] | [0.5, 2000] |
| Polyhedron dims | hardcoded 30.0 | [5, 500] |
| Energy (non-EMSP) | [5.0, 25.0] | [1.0, 100.0] |
| Flux | [1e10, 1e14] | [1e6, 1e16] |
| EMSP energy | {100,200,300,400} | [20, 1000] |
| EMSP flux | [1e5, 1e8] | [1e3, 1e10] |
| Unit cell | [20, 200] | [3, 800] |
| num_residues | [20, 500] | [1, 5000] |
| num_monomers | [1, 48] | [1, 200] |
| Solvent fraction | [0.3, 0.85] | [0.01, 0.99] |
| PixelsPerMicron | fixed sets | continuous [0.001, 10.0] |
| Angular resolution | fixed set | continuous [0.1, 90.0] |
| Protein conc (SAXS) | [0.5, 50.0] | [0.01, 500.0] |
| DDM gamma | [0.1, 2.0] | [0.001, 100.0] |
| DDM beta | [0.01, 2.0] | [0.001, 50.0] |
| DDM b0 | [0.1, 10.0] | [0.001, 500.0] |

The cost model will naturally reject configs that are too expensive, so wider ranges are safe. Configs that were previously impossible to generate (e.g., sub-micron crystals, high-energy beams) can now be explored.

**File**: `generate.py` lines 403-467 and scattered throughout `_sample()` and `_sample_beam()`.

### 4. Flatten Categorical Weights

Current weights strongly bias common paths. Flatten to increase coverage of rare paths:

| Choice | Current | New |
|--------|---------|-----|
| Subprogram | 4:1:1:1 (57% std) | 2:1:1:1 (40% std) |
| Crystal type | 3:1:1:1 (50% Cuboid) | 2:1:1:1 (33% Cuboid) |
| CoefCalc | 3:1:1 (60% RD3D) | 2:1:1 (40% RD3D) |
| DDM | 2:1:1:1 (40% Simple) | 1:1:1:1 (25% each) |
| Beam type | 50/50 | 50/50 (unchanged) |
| Segments | 3:2:1 (50% single) | 2:2:1 (40% single) |

Also remove the forced constraints that limit exploration:
- Cylinder currently forces SAXSseq — change to SAXSseq (2x), RD3D, CIF
- Polyhedron currently forces RD3D — change to RD3D (2x), SMALLMOLE, CIF

**File**: `generate.py` lines 383-401, 462-467, 485. Mirror in `strategies.py`.

### 5. Bootstrap Seed Corpus

Create `corpus/seeds/` with 10-15 representative input files covering every major code path. These enable the mutation strategy which currently has nothing to work with:

- `standard_cuboid.txt` — basic Cuboid + RD3D
- `standard_spherical.txt` — Spherical + RD3D
- `standard_cylinder.txt` — Cylinder + SAXSseq
- `standard_polyhedron.txt` — Polyhedron + RD3D
- `montecarlo_basic.txt` — MC with PE escape
- `montecarlo_fl_escape.txt` — MC with FL escape
- `xfel_basic.txt` — XFEL mode
- `microed_basic.txt` — EMSP/MicroED
- `smallmole.txt` — SMALLMOLE coefcalc
- `cif_input.txt` — CIF coefcalc
- `multi_wedge.txt` — multiple wedge/beam segments
- `heavy_atoms.txt` — many heavy protein atoms + solvent
- `ddm_leal.txt` — Leal dose decay model
- `ddm_bfactor.txt` — Bfactor dose decay model
- `large_crystal.txt` — large crystal near voxel cap

Source these from existing test fixtures in `raddose3d/tests/fixtures/` and `examples/`, plus hand-craft a few to fill gaps.

**Files**: New directory `fuzz/corpus/seeds/`, ~15 `.txt` files.

### 6. Text-Level Mutation Improvements in `generate.py`

Enhance `mutate_text()` beyond numeric perturbation:

- **Keyword swap** (10% chance): replace a crystal type keyword with another (e.g., "Cuboid" -> "Spherical"), or a beam type, or DDM type
- **Line deletion** (5% chance per line): remove a random non-block-header line
- **Line duplication** (5% chance): duplicate a random line (tests duplicate keyword handling)
- **Comment injection** (3% chance): insert `# fuzz` comment lines (tests parser robustness)
- **Case mutation** (5% chance per keyword): randomize case of a keyword
- **Whitespace mutation** (10% chance): add extra spaces, tabs, or blank lines
- **Block reordering** (5% chance): swap Crystal/Beam/Wedge block order

**File**: `generate.py`, extend `mutate_text()` function (~60 additional lines).

### 7. Container Block Fuzzing

The grammar generator never generates Container blocks, missing an entire code path. Add container generation:

- 30% chance of adding a Container block after Crystal
- Container types: `MIXTURE` with material name, `ELEMENTAL` with element/count pairs, or `NONE`
- Materials: common ones from the codebase (water, ethanol, paraffin oil, etc.)
- Container thickness: [1, 1000] um
- Container density: [0.5, 5.0] g/cm3

**File**: `generate.py`, new `_sample_container()` method, called from `_sample()`. Add container field to Config dataclass and rendering.

### 8. Coverage Tracking and Guided Generation

Add a simple coverage-feedback mechanism to steer generation toward under-explored paths:

- Track a **coverage map**: (subprogram, crystal_type, coefcalc, ddm, has_container, n_segments) tuples seen so far
- After every N iterations (e.g., 50), compute which tuples have 0 or low coverage
- Bias the next batch of grammar generations toward under-covered tuples by temporarily overriding the categorical weights
- Log coverage stats alongside existing statistics

This is not full coverage-guided fuzzing (no code instrumentation), but a lightweight feedback loop that ensures all reachable combinations get tested.

**File**: New class `CoverageTracker` in `generate.py` (~50 lines). Integrate into `fuzz.py` main loop.

### 9. Negative/Invalid Value Testing

Add a `"chaos"` generation mode (~10% of iterations when strategy is `"both"`) that deliberately generates potentially-invalid inputs to test error handling:

- Negative dimensions, negative energy, negative flux
- Zero PixelsPerMicron
- Wedge end < start
- Missing required fields (no energy, no flux)
- Unknown keywords (`CRYSTALTYPE FooBar`)
- Empty blocks (`Crystal\nEnd`)
- Extremely long element lists (20+ elements)

Expected behavior: both Java and Rust should either produce matching output or both crash/error gracefully. The interesting case is when one crashes and the other doesn't.

**File**: New function `generate_chaos()` in `generate.py` (~40 lines). Wire into `fuzz.py` strategy selection.

### 10. Sync `strategies.py` with All Changes

Mirror all changes from `generate.py` into the Hypothesis strategies so that `test_differential.py` benefits from the same improvements. This includes:
- Wider ranges
- Flatter weights
- Container blocks
- Boundary value injection (as a Hypothesis `map()` post-processor)

**File**: `strategies.py`, update `raddose_config()` composite strategy to match.

## File Summary

| File | Changes |
|------|---------|
| `fuzz/generate.py` | Boundary injection, structural mutation, wider ranges, flatter weights, container blocks, enhanced text mutation, coverage tracker, chaos mode |
| `fuzz/strategies.py` | Mirror all generate.py changes for Hypothesis |
| `fuzz/fuzz.py` | New strategy options (structural, chaos), coverage tracking integration, updated strategy selection |
| `fuzz/compare.py` | No changes needed |
| `fuzz/harness.py` | No changes needed |
| `fuzz/test_differential.py` | No changes needed (benefits automatically from strategies.py) |
| `fuzz/corpus/seeds/*.txt` | ~15 new seed files |

## Implementation Order

1. **Widen ranges + flatten weights** (items 3, 4) — smallest diff, biggest immediate coverage gain
2. **Bootstrap seed corpus** (item 5) — enables mutation path
3. **Boundary value injection** (item 1) — targets edge cases
4. **Enhanced text mutation** (item 6) — more mutation diversity
5. **Container block fuzzing** (item 7) — new code path coverage
6. **Structural mutation** (item 2) — novel config shapes
7. **Negative/chaos mode** (item 9) — robustness testing
8. **Coverage tracking** (item 8) — guided exploration
9. **Sync strategies.py** (item 10) — keep Hypothesis tests in sync (do incrementally alongside each step)

## Verification

After each step:
1. `cd /workspace/fuzz && python -c "from generate import GrammarGenerator; g = GrammarGenerator(); print(g.generate())"` — verify generation still works
2. `python fuzz.py --iterations 20 --strategy both --timeout 30` — smoke test with short run
3. `python -m pytest test_differential.py --hypothesis-seed=0 -x -k test_no_crashes` — verify Hypothesis integration
4. Inspect generated inputs manually to confirm wider diversity
5. After all steps: full run with `--iterations 500` and compare category distribution against baseline (expect more MINOR_DIFF/MAJOR_DIFF discoveries, broader triple coverage)
