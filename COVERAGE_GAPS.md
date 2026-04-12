# Coverage Gap Report

Generated from `coverage_local4.json` (67% line coverage overall across workspace crates).

## Summary

| Category | Approx lines | % of total | Action |
|----------|-------------|-----------|--------|
| Structurally unreachable from CLI | ~375 | 8% | None — by design |
| Cryo/surrounding escape cluster | ~1,100 | 23% | Fix seed wiring |
| M-shell heavy atom branches | ~80 | 2% | Add seed |
| GOS inner-shell calculations | ~200 | 4% | Add seed |
| MicroED cylinder/polyhedron geometry | ~125 | 3% | Add seed |
| XFEL edge cases | ~160 | 3% | Add seed |
| Small scattered gaps | ~600 | 13% | Partially addressable |
| **Total uncovered** | **~4,702** | **33%** | |

---

## 1. Structurally unreachable from CLI (~375 lines, 8%)

These paths cannot be reached by any seed. They are covered by unit tests or require a different entry point.

| File | Lines | Reason |
|------|-------|--------|
| `raddose3d/src/output/progress.rs` | 143 | `OutputProgressIndicator` never instantiated in CLI output chain. Progress is printed via hardcoded `print!()` in `expose_rd3d`; this class is dead from CLI. |
| `raddose3d/src/output/fluence_per_dose_hist.rs` | 137 | `OutputFluencePerDoseHistCSV` never instantiated from CLI default output set. |
| `raddose3d/src/lib.rs` | 61 | WASM / external library API (`parse_input()`, `run()`, `parse_input_json()`). Only reachable via WASM build or as a library crate. |
| `raddose3d/src/simulation/java_random.rs` | 34 | Only used by unit tests in `tests/phase6_validation.rs`. `mc.rs` uses `rand::random()` at runtime, not `JavaRandom`. |

---

## 2. Cryo/surrounding escape cluster (~1,100 lines, 23%)

This is the single largest gap. All of it is gated on `CalcSurrounding TRUE` + `CalculatePEEscape TRUE` being active simultaneously, *and* on `is_cryo()` returning `true` from the coefcalc. The `pe_escape_surrounding.txt` seed exists but may not be triggering `is_cryo()` correctly for `CoefCalcFromParams` — needs investigation.

### Affected files

**`raddose3d/src/crystal/escape.rs`** — 661 uncovered lines (86% of file)

Key uncovered regions:
- Lines 133–148: `get_pe_energy_to_subtract` — PE energy loss accounting from fe_factors
- Lines 157–213: `get_auger_energy` — Auger cascade energy from K/L/M shells
- Lines 222–234: `get_fl_energy_release` — FL energy partition across shells
- Lines 490–512: `set_cryo_ppm` — cryo grid pixels-per-micron calculation
- Lines 514+: `setup_cryo_escape` — full cryo surrounding PE setup

**`raddose3d/src/crystal/mod.rs`** — 325 uncovered lines (54% of file)

Key uncovered regions:
- Lines 175–218: `setup_pe_escape` + cryo branch (`update_cryo_coefficients`, `cryo_fluorescent_escape_factors`)
- Lines 272–285: `expose_cryo_angle` call per rotation angle
- Lines 319–323: PE escape summary print
- Lines 559–570: `apply_escape_fl` call per voxel
- Lines 580–583: `apply_escape_pe` call per voxel
- Lines 637–699: `apply_escape_pe` and `apply_escape_fl` function bodies

**`raddose3d/src/coefcalc/compute.rs`** — ~100 uncovered lines in cryo path

Key uncovered regions:
- Lines 538–539: cryo occurrence accumulation
- Lines 606–616, 641–653: `update_cryo_coefficients` — cryo photo/coherent/compton coefficients
- Lines 689–690+: `cryo_fluorescent_escape_factors` — cryo FL factors per element
- Lines 728–730: L-edge fluorescent for cryo elements

**`raddose3d/src/coefcalc/mod.rs`** — 15 uncovered lines

Key uncovered regions:
- Lines 123–125: `photo_electric_probs_element_surrounding` default stub
- Lines 205–208: `gos_outer_shell_probs` default stub
- Lines 259–262: `get_relative_shell_probs_for_element` default stub
- Lines 273–276: `relative_shell_probs` default stub

**`raddose3d/src/simulation/mc.rs`** — ~100 uncovered lines in surrounding path

Key uncovered regions:
- Lines 370–374: lazy init of `normals_surrounding`
- Lines 451–460: surrounding vertex expansion (offset by `surrounding_thickness`)
- Lines 506–513: surrounding vertex rotation per angle
- Lines 674–681: `get_intersection_point` for surrounding geometry

### Diagnosis

The seed `corpus/seeds/pe_escape_surrounding.txt` uses `CalcSurrounding TRUE` + `CalculatePEEscape TRUE` + `Subprogram MONTECARLO`. Check whether `CoefCalcFromParams::is_cryo()` returns `true` when `CalcSurrounding TRUE` is set — if not, the cryo branch is silently skipped.

---

## 3. M-shell heavy atom branches (~80 lines, 2%)

**`raddose3d/src/coefcalc/compute.rs`** lines 770–820

Gated on: element present in crystal/solution with atomic number ≥ 73 (W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi) AND beam energy above the element's M2/M3/M4/M5 edges (typically 2–4 keV, so any standard beam energy suffices).

The existing `heavy_atoms.txt` seed uses Zn (Z=30) and S/Se — all below Z=73. A seed with `SolventHeavyConc Pt 100` or `ProteinHeavyAtoms Pb 1` would cover these branches.

**Suggested fix:** Add a seed with `ProteinHeavyAtoms Pt 2` or `SolventHeavyConc Pb 50` at standard energy (e.g. 12.1 keV).

---

## 4. GOS inner-shell calculations (~200 lines, 4%)

**`raddose3d/src/simulation/mc.rs`** lines 1142–1171, 1678–1696

Key uncovered regions:
- `get_relative_shell_probs`: inner/outer shell CDF construction for GOS — requires L/M-shell ionisation (Z > ~12 with edges below beam energy)
- Auger cascade inner-shell branch: `gos_electron_dose_v_resolved` accumulation

**`raddose3d/src/simulation/xfel.rs`** lines 1484–1510, 1532–1533, 1569–1570, 1684–1698, 1736–1738

Key uncovered regions:
- Inner-shell Auger lifetime + re-emission: lines 1484–1510
- `ionisations_old` counter (non-surrounding path): line 1532–1533
- `dose_simple` accumulation: lines 1569–1570
- GOS electron dose inner-shell accumulation: lines 1684–1698
- `inner_shell_lambda` branch in GOS inelastic: lines 1736–1738

Both `montecarlo_gos.txt` and `xfel_basic.txt` seeds exist but use light elements (C/N/O/H/S). GOS inner-shell logic fires when `gos_inner_lambda > 0` which requires heavier elements with inner-shell cross-sections. **Suggested fix:** Add GOS seed with `ProteinHeavyAtoms Fe 4` or `SolventHeavyConc Fe 1000` — Fe (Z=26) has well-tabulated inner shells.

---

## 5. MicroED cylinder/polyhedron geometry (~125 lines, 3%)

**`raddose3d/src/simulation/micro_ed.rs`**

Key uncovered regions:
- Lines 53–56, 64–67: `cross_product` / `normalised_cross` helpers — only called during mesh normal calculation for non-Cuboid crystals
- Lines 307–308: cylinder Z-dimension adjustment (`z_dimension -= 1e-6`)
- Lines 317–319: cylinder surface area formula (uses π × semi-axes)
- Lines 531–532: CSDA stopping power integration loop (`distance += energy_step / stopping_power`)
- Lines 561–562: `number_slices == 0` fallback (set to 1)
- Lines 578–580: low-energy termination branch (`avg_energy < 0.05 keV`)
- Lines 856–867: `calculate_normals` for mesh triangles
- Lines 900–906: `polygon_inclusion_test` ray-casting in mesh occupancy check
- Lines 943–945: lazy `calculate_crystal_occupancy` (non-Cuboid voxel occupancy)

All mesh geometry paths (856–945) are only reached when the crystal is not a Cuboid. The current `microed_basic.txt` and `microed_emed.txt` seeds both use `Type Cuboid`. The cylinder-specific branches (307–319) need `Type Cylinder`.

**Suggested fix:** Add `Type Cylinder` + `Subprogram EMSP` seed, and `Type Polyhedron` + `Subprogram EMSP` seed.

---

## 6. XFEL edge cases (~160 lines, 3%)

**`raddose3d/src/simulation/xfel.rs`**

Key uncovered regions:
- Lines 461–471: `get_intersection_distance` — surrounding geometry intersection for XFEL. Needs `CalcSurrounding TRUE` + `Subprogram XFEL`.
- Lines 637–638: `high_energy_angles` lazy load — ELSEPA high-energy table for heavy elements; needs high-Z element + XFEL.
- Lines 994–997: non-pulsed XFEL path — `PulseEnergy` not provided, flux × total_sec used instead. The `xfel_basic.txt` seed specifies `PulseEnergy 2e-6`, bypassing this branch.

**Suggested fix:** Add XFEL seed without `PulseEnergy` keyword. Add XFEL seed with `CalcSurrounding TRUE`.

---

## 7. Small scattered gaps (~600 lines, 13%)

These are individually small uncovered regions spread across many files. Most are error paths or rarely-triggered branches:

| File | Approx lines | Nature |
|------|-------------|--------|
| `crystal/polyhedron.rs` | 114 | OBJ error paths, edge cases in ray-casting near degenerate triangles |
| `element/database.rs` | 73 | Element lookup misses, unused cross-section table branches |
| `coefcalc/from_pdb.rs` | 65 | RCSB download path (network), SEQRES parsing edge cases |
| `raddose3d-parser/src/lib.rs` | 64 | Unparsed keyword variants, multi-beam parsing paths |
| `coefcalc/micro_ed.rs` | 58 | MicroED-specific coefcalc branches (stopping power for heavy elements) |
| `coefcalc/small_molecules.rs` | 56 | Small-molecule atom count edge cases |
| `coefcalc/from_sequence.rs` | 42 | Sequence file parse edge cases (non-standard FASTA) |
| `raddose3d-cli/src/main.rs` | 40 | Error-exit paths, `-r` flag (RDFortran subprocess), stdin `-i -` path |
| `writer.rs` | 38 | `TeeWriter` and `NullWriter` paths (not used in default output chain) |
| `element/database_em.rs` | 32 | EM element lookup misses |
| `coefcalc/from_sequence_saxs.rs` | 32 | SAXS+sequence edge cases |
| `coefcalc/saxs.rs` | 31 | SAXS solvent/heavy-atom branches |
| `output/exposure_summary.rs` | 29 | Edge statistics (zero-dose quantile, empty voxel set) |
| `coefcalc/from_cif.rs` | 28 | CIF formula parse edge cases |
| `experiment.rs` | 28 | Multi-crystal error path (lines 122–124), second crystal block |
| `beam/mod.rs` | 27 | Remaining experimental beam branches |
| `coefcalc/from_params.rs` | 26 | Non-standard unit cell angles (non-orthogonal) |
| `beam/tophat.rs` | 26 | Circular collimation + tophat edge cases |
| `container.rs` | 33 | ContainerMixture path (NIST lookup stub) |
| `beam/gaussian.rs` | 19 | Circular collimation edge cases |
| `residue.rs` | 23 | Unknown residue fallback, RNA/DNA lookup |

Many of these require either: invalid/edge-case inputs (negative cell volume), network access (RCSB download), or CLI flags not exercised by the fuzzer (`-r`, stdin).

---

## Actionable seed additions

Priority order based on lines gained per seed:

| Priority | Seed to add | Target gap | Est. lines |
|----------|------------|-----------|-----------|
| 1 | Investigate `is_cryo()` wiring in `pe_escape_surrounding.txt` | Cryo/surrounding cluster | ~1,100 |
| 2 | `Type Cylinder` + `Subprogram EMSP` microED | MicroED mesh geometry | ~80 |
| 3 | `ProteinHeavyAtoms Pt 2` at standard energy | M-shell branches | ~80 |
| 4 | `Subprogram GOS` with heavy atoms (`SolventHeavyConc Fe 1000`) | GOS inner-shell | ~100 |
| 5 | XFEL seed without `PulseEnergy` | XFEL non-pulsed path | ~10 |
| 6 | XFEL seed with `CalcSurrounding TRUE` | XFEL surrounding | ~50 |
| 7 | Multi-crystal input (two Crystal blocks) | `experiment.rs` error path | ~5 |
| 8 | Non-orthogonal unit cell angles | `from_params.rs` cell volume | ~5 |
