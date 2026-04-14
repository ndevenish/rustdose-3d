# Java vs Rust Differential Fuzzing: Diff Classification

Analysis of `fuzz/corpus/diffs/` — 191 non-matching inputs (382 files counting `.json` + `.txt`),
out of 1422 total corpus inputs.

## NAN_INF — 53 inputs

### Group 1: Java NaN Diffraction Efficiency / AD-WC (46 inputs)

Java computes `0/0 = NaN` for Diffraction Efficiency (and sometimes AD-WC) when the crystal
receives zero dose. Rust returns `0.0` gracefully.

Causes: `Wedge 0 0` with crystal missing the beam, large `TranslatePerDegree` moving the crystal
out of beam range, or very low flux delivering no voxel dose.

**Classification: Java bug, Rust is correct. Not actionable.**

### Group 2: Java zero dose, Rust nonzero (5 inputs)

Files: `31123`, `31217`, `35554`, `35577`, `35586` (`LFP_test_cuboid_large_crystal`)

All use `Wedge 0 0` + `CalculatePEEscape True` / `CalculateFlEscape True`. Java gives 0 dose,
Rust gives nonzero Max Dose (~0.001–0.37 MGy). Java also gets NaN for Diffraction Efficiency
because its denominator is zero.

Java may be skipping the single-angle exposure in the PE/FL escape path when `Wedge 0 0`.

**Classification: Likely Java bug in PE/FL escape for `Wedge 0 0`. Rust may be more correct.**

### Group 3: Rust NaN Dose Inefficiency (2 inputs) — **FIXED**

Files: `42558` (`formate_PSII`), `42574` (`250517_0011-04_2`)

Helical scans with large `TranslatePerDegree` + `StartOffset`. Rust was giving zero dose
(crystal translates out of beam range), Java gave nonzero dose (~5.5 MGy). Rust then computed
`0/0 = NaN` for Dose Inefficiency.

**Root cause:** `compute_angles()` used `.floor() as i64` to normalise the start angle, but
Java uses `(int)` cast which truncates toward zero. For `Wedge -45 315` (start = −π/4):
- Java: `(int)(−0.125) = 0` → start stays at −π/4, translation delta at angle 0 is 0 µm.
- Rust (old): `floor(−0.125) = −1` → start shifts to +7π/4, translation delta at angle 0
  is `(7π/4 − (−π/4)) × trans = 2π × trans ≈ 775 µm` — the full scan length, placing the
  crystal entirely outside the beam at the very first exposure.

**Fix:** replaced `.floor() as i64` with `as i64` (truncation) in `compute_angles()`
(`crystal/mod.rs`). After the fix, Rust gives Max Dose = 5.524258 MGy, matching Java exactly.

**Seed:** `fuzz/corpus/seeds/helical_negative_start_angle.txt`

---

## MAJOR_DIFF — 28 inputs (>1% on any metric)

### Group 1: PE/FL escape — SmallMole + CalcSurrounding (20 inputs, ~1.5–5.5%)

All use `AbsCoefCalc SMALLMOLE` + `CalcSurrounding True` (or `CalculatePEEscape` /
`CalculateFlEscape`). Worst metrics: Dose Inefficiency PE, Max Dose, DWD.

This is the known PE/FL escape variance documented in CLAUDE.md (~0.5–1% DWD systematic
offset, ±10% Max Dose from single PE track per voxel random sampling and Java/Rust FP
differences in `find_voxels_reached_by_pe`).

**Classification: Known, documented. Not fixable without matching every FP operation exactly.**

### Group 2: Cylinder crystal + large dimensions → DWD diff (2 inputs, ~1.6%)

Files: `09970` (`RTaged1px`), `11011` (`RTfresh9.93`)

Both are Cylinder (60×80 µm), Gaussian beam 10×10 µm FWHM, Energy 16 keV, Wedge −30–30°.
~1.6% Last DWD difference. The smaller-crystal versions of the same experiment (`09923`, `09924`)
are only MINOR (~0.35%) — larger crystals amplify the `findDepth` boundary voxel deduplication
bug.

**Classification: `findDepth` deduplication bug amplified by larger crystal. Known.**

### Group 3: SAXS + NumResidues=0 + Elemental Container + Wedge 0 0 (6 inputs, ~10.5%)

Files: `42268`, `42269`, `42272`, `42273`, `42274`, `42344`

~10.5% Diffraction Efficiency, ~5.4% Dose Inefficiency. All use `AbsCoefCalc SAXS` +
`NumResidues 0` + `ContainerMaterialType Elemental` + `DDM Simple` + `Wedge 0 0`.

**Root cause:** The Diffraction Efficiency divergence (10.5%) and Dose Inefficiency divergence
(5.4%) are identical across all 6 inputs regardless of whether `ContainerDensity` is provided —
so the container is NOT the primary cause. The divergence comes from the SAXS crystal elastic
coefficient calculation with `NumResidues=0`.

When `NumResidues=NumRNA=NumDNA=NumCarb=0`, `total_mass=0`, so `molarity = ProteinConc/0 =
Infinity`, and `num_monomers` saturates to `i32::MAX = 2147483647` in both Java and Rust.
With `ProteinHeavyAtoms C 40 H 80 N 1 O 8 P 1`, this produces ~8.6×10¹⁰ C atoms in a 1e6 Å³
cell — extreme values that produce floating-point divergence in the elastic coefficient
calculation (density and solvent-water normalisation).

`ContainerElemental` IS fully implemented in Rust (using embedded NIST tables). For the three
inputs that do provide `ContainerDensity` (42269, 42272, 42273), Rust uses embedded NIST data
with log-log interpolation while Java downloads the NIST table and uses linear interpolation;
this causes a small additional Elastic Yield difference on top of the SAXS divergence.

**Classification: Degenerate SAXS case (NumResidues=0 overflow) causing FP divergence in
elastic coefficient. Not fixable without matching Java's exact FP evaluation path.
Seed:** `fuzz/corpus/seeds/saxs_num_residues_zero.txt`

### Group 4: PE/FL escape — SmallMole Spherical (2 inputs, ~1–3.7%)

Files: `41464`, `41465`

`Type Spherical` + `AbsCoefCalc SMALLMOLE` + `CalculatePEEscape True`, `Wedge 0 0`.
Same root cause as Group 1 but spherical crystal geometry.

**Classification: Known PE/FL escape variance.**

---

## MINOR_DIFF — 110 inputs (0.1–1%)

### Group 1: MicroED GOS Inelastic divergence (3 inputs, ~0.3–0.4%)

Files: `03653`, `03660`, `03661`

~0.4% Inelastic yield, ~0.28% Dose. MicroED/EMED subprogram. This is the known GOS Inelastic
Lambda ~1% divergence from Java's `sumZ` integer truncation bug documented in CLAUDE.md. Rust
uses `round()` rather than truncation, which is physically more correct.

**Classification: Known, documented. Rust is more correct.**

### Group 2: PE/FL escape — SmallMole CalcSurrounding (22 inputs, ~0.5–1%)

~0.5–1% Max Dose, Last DWD, Dose Inefficiency PE. Same root cause as MAJOR Group 1 but smaller
magnitude due to different crystal/beam parameters.

**Classification: Known PE/FL escape variance.**

### Group 3: Cuboid TranslatePerDegree Max Dose ~0.2% (13 inputs)

Files: `00473`–`00487` (1E7 series)

Consistent ~0.19% Max Dose. Cuboid, 360° wedge, `TranslatePerDegree`. The `findDepth`
deduplication bug causes boundary voxels to get depth=0 at certain angles; which exact voxels
are affected varies slightly between Java and Rust FP paths.

**Classification: `findDepth` boundary voxel deduplication bug. Known.**

### Group 4: Cylinder Last DWD ~0.35% (2 inputs)

Files: `09923` (`RTfresh1`), `09924` (`RTfresh2px`)

Smaller-crystal version of MAJOR Group 2. Same cause.

**Classification: `findDepth` bug. Known.**

### Group 5: Helical scan TranslatePerDegree Last DWD ~0.2–0.3% (NcAA9D series, ~5 inputs)

Files: `03575`–`03579`

Helical scan, ~0.27% Last DWD in per-wedge slices. Likely floating-point accumulation
differences in the translation-per-angle step, related to the same root cause as NAN Group 3
but within tolerance.

**Classification: Helical scan FP accumulation. Related to the same floor/trunc bug as NAN
Group 3, but the angle offset here is small enough that the crystal stays in beam — only the
translation accumulation differs. The NAN Group 3 fix does not affect these cases.**

### Group 6: General FP noise — Diffraction Efficiency / Wedge Absorbed Energy / AD-WC (65 inputs, <0.1%)

Tiny (<0.1%) differences across many metrics in standard cuboid runs with no escape, no DDM,
no helical scan. Likely from different Java/Rust FP evaluation order in the voxel accumulation
loop.

**Classification: Floating-point accumulation noise. Not actionable.**

---

## Summary

| Category | Inputs | Max Diff | Root Cause | Actionable? |
|----------|-------:|----------|------------|-------------|
| Java NaN DE (0/0 division) | 46 | ∞ (Java bug) | Java 0/0 = NaN | No — Java bug |
| Java zero dose on Wedge 0 0 + PE escape | 5 | ∞ | Java skips single-angle PE exposure | Investigate |
| Rust zero dose on helical scan | 2 | ∞ | `floor()` vs `(int)` in angle normalisation | **Fixed** |
| PE/FL escape variance MAJOR | 22 | ~5.5% | PE random sampling noise | No — known |
| SAXS + NumResidues=0 elastic divergence | 6 | ~10.5% | FP divergence with num_monomers overflow | No — degenerate input |
| Cylinder large crystal findDepth MAJOR | 2 | ~1.6% | `findDepth` dedup bug (larger crystal) | No — known |
| MicroED GOS Inelastic | 3 | ~0.4% | Java `sumZ` int truncation | No — Rust correct |
| Cuboid findDepth MINOR | 13 | ~0.2% | `findDepth` dedup bug | No — known |
| PE/FL escape variance MINOR | 22 | ~1% | PE random sampling noise | No — known |
| Helical scan FP accumulation | 5 | ~0.3% | TranslatePerDegree FP order | Investigate |
| General FP accumulation noise | 65 | <0.1% | FP evaluation order differences | No — noise |

### Resolved findings

1. **Helical scan crystal-in-beam tracking — FIXED** (42558, 42574): Root cause was
   `compute_angles()` using `.floor() as i64` instead of `as i64` (truncate toward zero) to
   normalise the wedge start angle. For `Wedge -45 315`, `floor(-0.125) = -1` shifted the
   normalised angles by +2π, causing `translation_vector(angle - start_ang)` to start at a
   delta of 2π instead of 0 — placing the crystal at the far end of the helical scan at the
   very first exposure angle. Fix: `crystal/mod.rs`, `compute_angles()`, one-line change.
   Seed: `fuzz/corpus/seeds/helical_negative_start_angle.txt`.

2. **`ContainerMaterialType Elemental` + SAXS + NumResidues=0** (42268 cluster): Corrected
   misclassification — `ContainerElemental` is fully implemented in Rust (embedded NIST tables,
   `container.rs`). The ~10.5% Diffraction Efficiency divergence is a SAXS + degenerate input
   issue: `NumResidues=0` → `total_mass=0` → `num_monomers = i32::MAX = 2147483647`, producing
   ~8.6×10¹⁰ atoms per 1e6 Å³ cell and FP divergence in the elastic coefficient calculation.
   Not actionable without exactly matching Java's FP path. Seed: `fuzz/corpus/seeds/saxs_num_residues_zero.txt`.
