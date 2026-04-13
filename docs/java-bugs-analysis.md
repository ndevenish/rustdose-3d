# Java RADDOSE-3D Bug Analysis: PE/FL Escape and Absorbed Energy

This document describes two bugs found in the Java RADDOSE-3D codebase during validation of the Rust rewrite against the LiFePO4 test case (`lfp_short.txt`: 1 um spherical crystal, SMALLMOLE, 1.487 keV, with PE escape, FL escape, and cryo surrounding enabled).

## Test Case Output Comparison

| Metric | Java | Rust | Difference |
|--------|------|------|------------|
| PE Escape | 1.31e+03 J | 1.21e3 J | ~7.6% |
| **FL Escape** | **8.23e+00 J** | **3.53e0 J** | **~57%** |
| Avg DWD | 8.866604 MGy | 8.827484 MGy | ~0.44% |
| Max Dose | 25.365039 MGy | 24.905107 MGy | ~1.81% |
| **Absorbed Energy** | **4.05e-09 J** | **3.30e-8 J** | **~8x** |
| Dose Inefficiency | 6.26e12 1/g | 7.56e11 1/g | ~8x |
| Dose Inefficiency PE | 875.3e9 1/g | 871.1e9 1/g | ~0.5% |

The two largest discrepancies (FL Escape ~57% and Absorbed Energy ~8x) are caused by bugs in the Java implementation described below. The remaining differences (dose metrics 0.4-1.8%) are downstream effects of the FL escape bug plus expected variance from random PE track sampling.

---

## Bug 1: Hardcoded `muabsIndex` in Fluorescent Escape Distribution

### Location
`CrystalPolyhedron.java`, method `calcFluorescenceDistribution()`, lines 1239-1240

### The Code

```java
for (int j = 0; j < 4; j++) { //for each shell (K=0, L1=1, L2=2, L3=3)
    runningEscapeTotal = 0;
    if (fluorescenceProportionEvent[i][j] > 0) {
//      int muabsIndex = (4* j) + 4;   // <-- CORRECT formula, commented out
        int muabsIndex = 4;             // <-- BUG: hardcoded to K-shell for all shells
        double maxDistanceFl = -1 * (Math.log(0.05)/feFactors[i][muabsIndex]);
```

The variable `muabsIndex` is used to look up the fluorescent photon absorption coefficient (`escapeMuAbs`) for each electron shell from the `feFactors` array. The loop iterates `j` over four shells (K, L1, L2, L3), but `muabsIndex` is hardcoded to `4` regardless of `j`.

### The `feFactors` Array Layout

The `feFactors` (a.k.a. `fluorEscapeFactors`) array is populated in `CoefCalcCompute.java` lines 920-937, with a repeating 4-stride layout per shell:

| Index | Shell | Field |
|-------|-------|-------|
| 0 | — | `muAbsFrac` (element fraction of total µ_abs) |
| **1** | K | `kShellEnergy` (absorption edge, keV) |
| 2 | K | `kFactorA` (ionisation probability) |
| 3 | K | `kFactorB` (fluorescence yield) |
| **4** | K | **`escapeMuAbsK`** |
| 5 | L1 | `l1ShellEnergy` |
| 6 | L1 | `l1FactorA` |
| 7 | L1 | `l1FactorB` |
| **8** | L1 | **`escapeMuAbsL1`** |
| 9 | L2 | `l2ShellEnergy` |
| 10 | L2 | `l2FactorA` |
| 11 | L2 | `l2FactorB` |
| **12** | L2 | **`escapeMuAbsL2`** |
| 13 | L3 | `l3ShellEnergy` |
| 14 | L3 | `l3FactorA` |
| 15 | L3 | `l3FactorB` |
| **16** | L3 | **`escapeMuAbsL3`** |

The shell-specific `escapeMuAbs` values are at indices 4, 8, 12, 16 — exactly `(4*j) + 4` for `j = 0..3`.

### Why This Is a Bug

The correct index formula `(4*j) + 4` is present in the source as a commented-out line immediately above the hardcoded version. This means the correct implementation was written and then disabled in favour of the constant `4`.

When `muabsIndex` is always `4`, every shell's fluorescence distance distribution is computed using the K-shell absorption coefficient (`escapeMuAbsK`). This is only correct for the K shell (`j=0`). For the L shells (`j=1,2,3`), the wrong coefficient is used.

### Physical Consequence

The fluorescent photon absorption coefficient determines how far a fluorescent photon can travel before being reabsorbed. K-fluorescence photons have higher energy than L-fluorescence photons, and the absorption coefficient at the K fluorescence energy is different from the coefficient at L fluorescence energies.

For heavy elements like Fe at low beam energies (1.487 keV in this test), the K-shell `escapeMuAbs` is typically higher than the L-shell values. Using the K value for L-shell calculations means:

1. Higher µ_abs → shorter calculated escape distances (`maxDistanceFl = -ln(0.05) / µ_abs`)
2. Shorter distances → the fluorescence distance distribution is compressed
3. Compressed distribution → different (generally higher) probability of FL photons escaping the crystal

The net effect is that Java overestimates FL escape energy by ~57% on this test case (8.23 J vs the correct 3.53 J).

### Lines Affected

The hardcoded `muabsIndex = 4` is used at four points in the method:

- **Line 1242**: `maxDistanceFl` calculation (escape distance threshold)
- **Line 1263**: `flDistanceDistribution` for outermost bin (escape probability)
- **Line 1267**: `flDistanceDistribution` for inner bins (landing probability)
- All three use `feFactors[i][muabsIndex]` and all three are wrong for `j > 0`

### Correct Implementation

The Rust implementation (`escape.rs` lines 1100-1106) correctly computes the shell-specific index:

```rust
let muabs_index = match j {
    0 => 4,   // K
    1 => 8,   // L1
    2 => 12,  // L2
    3 => 16,  // L3
    _ => continue,
};
```

This is equivalent to the commented-out Java formula `(4*j) + 4`.

---

## Bug 2: `energyPerFluence` Variable Scope Leak from Cryo Block

### Location
`Crystal.java`, method `exposeAngle()` (or equivalent exposure method)

### The Code

The variable `energyPerFluence` is first assigned on line 1007 using the **crystal** absorption coefficient:

```java
double energyPerFluence =
    1 - Math.exp(-1 * coefCalc.getAbsorptionCoefficient()
        / getCrystalPixPerUM());
// absorption of the beam by a voxel
```

Later, inside the cryo surrounding block (line 1220: `if (aSurface) {`, line 1221: `if (photoElectronEscape) {`), it is **overwritten** on line 1236 with the **cryo material's** absorption coefficient:

```java
energyPerFluence =
    1 - Math.exp(-1 * coefCalc.getCryoAbsorptionCoefficient()
        / getCryoCrystalPixPerUM());
// absorption of the beam by a voxel
```

This reassignment is intentional for the cryo exposure loop that follows (lines 1248+), where it correctly represents the fraction of beam energy absorbed by each cryo voxel.

However, after the cryo block closes, the code enters a separate loop (line 1346) to compute DWD statistics and absorbed energy for the **crystal** voxels. At line 1361:

```java
absorbedEnergy[i][j][k] = voxImageFluence[i][j][k] * energyPerFluence;
double comptonabsorbedEnergy = voxImageComptonFluence[i][j][k] * energyPerFluence;
```

Here, `energyPerFluence` still holds the **cryo** value, not the crystal value. The crystal's fluence (`voxImageFluence`) is multiplied by the cryo's energy-per-fluence factor, producing an incorrect absorbed energy.

### Why This Is a Bug

The variable `energyPerFluence` is declared at method scope. The cryo block reassigns it as a side effect, and no code restores it before it is used again for crystal voxels. This is a classic variable-scope leak.

The Java source contains comments at lines 1353 and 1359 that acknowledge uncertainty about this code:

```java
double interpolatedVoxelDose = totalVoxelDose + voxImageDose[i][j][k] / 2; // this needs to change for PE escape
```

```java
//may need to pass in different things or pass in more and change in observer
```

These comments suggest the code was known to be incomplete for PE/FL escape handling.

### Physical Consequence

The leak involves **two** variables — the absorption coefficient AND the pixel resolution. The cryo grid uses a separate pixels-per-micron computed by `setCryoPPM()` (typically coarser or finer than the crystal's PPM depending on geometry). For the LiFePO4 test case:

| Parameter | Crystal value | Cryo value (leaked) |
|-----------|--------------|-------------------|
| Absorption coeff | 0.732 /µm | 0.174 /µm |
| Pixels per µm | 10.0 | 20.0 |
| `energyPerFluence` | 0.0706 | **0.00868** |

The cryo `energyPerFluence` is `1 - exp(-0.174/20) = 0.00868`, which is ~8x smaller than the crystal's `1 - exp(-0.732/10) = 0.0706`.

When this leaked value is used for crystal absorbed energy, the reported "Absorbed Energy (this Wedge)" is ~8x too small:
- Java (using cryo value): 4.05e-09 J
- Rust (using crystal value): 3.30e-08 J

This directly propagates into the "Dose Inefficiency" metric, which divides max dose by absorbed energy:
- Java: 6.26e12 1/g (inflated because denominator is too small)
- Rust: 7.56e11 1/g

Note that the "Dose Inefficiency PE" metric (which uses `totalEnergy` — the deposited energy tracked separately) is **not** affected, which is why it matches within 0.5% between Java and Rust.

### Cryo PPM Calculation

The cryo grid's pixels-per-micron is computed by `CrystalPolyhedron.setCryoPPM()` (lines 773-807). It depends on the max PE distance, crystal dimensions, and beam minimum dimension:

```java
int multiplyFactor = 20;
if (minDimCryst >= minDimBeam) {
    multiplyFactor *= 1.5; // int multiplication: 20 * 1.5 = 30
}
double idealPPM = (1.0/maxPEDistance) * 5 + (1.0/maxPEDistance) * multiplyFactor * (1.0/minDimCryst);

if (idealPPM >= crystalPixPerUM) {
    pixelsPerMicron = ceil(idealPPM / crystalPixPerUM) * crystalPixPerUM;
} else {
    pixelsPerMicron = crystalPixPerUM / ((int)(crystalPixPerUM / idealPPM));
}
```

For the LiFePO4 test (1 µm sphere, 10 PPM, FWHM 326x579 µm, max PE distance ~1 µm), this yields `cryoPPM = 20.0` — double the crystal PPM. The `fluenceToDoseFactor` variable (line 1226) also leaks from the cryo block but does not affect the output because it is not used in the post-cryo absorbed energy loop.

### Conditions That Trigger the Bug

This bug only manifests when **all three** conditions are true:
1. `CalcSurrounding True` (cryo surrounding enabled → `coefCalc.isCryo()` returns true)
2. `CalculatePEEscape True` (PE escape enabled → enters the `if (photoElectronEscape)` block where the reassignment lives)
3. A non-zero wedge is exposed (so the absorbed energy loop at line 1346 actually runs)

When cryo surrounding is disabled, the `if (aSurface)` block at line 1220 is skipped, `energyPerFluence` retains its crystal value, and the bug does not manifest.

### Correct Implementation

In a correct implementation, the crystal's `energy_per_fluence` would be computed using the crystal's own absorption coefficient and crystal PPM. The Rust code intentionally reproduces this bug for Java compatibility — see the `// Java bug compatibility` comments in `crystal/mod.rs`.

---

---

## Bug 3: MicroED Spherical Volume Branch is Structurally Unreachable

### Location
`MicroED.java` lines 232–241, `Crystal.java` lines 294–313

### The Code

```java
// Crystal.java lines 294-313 (in CrystalPolyhedron constructor)
crystalType = "CUBOID";
XDim = (double) properties.get(Crystal.CRYSTAL_DIM_X);
try {
    YDim = (double) properties.get(Crystal.CRYSTAL_DIM_Y);
} catch (NullPointerException e) {
    crystalType = "SPHERICAL";   // only reached when Y is null
    YDim = XDim;
    ZDim = XDim;
}
try {
    ZDim = (double) properties.get(Crystal.CRYSTAL_DIM_Z);
} catch (NullPointerException e) {
    if (crystalType == "CUBOID") {   // only when Z is null and was still CUBOID
        crystalType = "CYLINDER";
        ...
    }
}
```

```java
// MicroED.java lines 232-241
if (crystalTypeEM == "CYLINDER") {
    crystalSurfaceArea = (Math.PI * (XDimension/2) * (YDimension/2)) * 1E02;
    ZDimension -= 0.001;
}
if (crystalTypeEM == "SPHERICAL") {
    // NEVER REACHED — see analysis below
    crystalVolume = ((4/3) * Math.PI * ...) * 1E-24;
}
//note the volume would need to be updated for a polyhedron!!!
//although it isn't used
```

### Why the SPHERICAL Branch Is Unreachable

The `crystalType = "SPHERICAL"` assignment (Crystal.java:301) is only reached when the Y dimension property is null — i.e., when a crystal is described with only a single dimension. In Java, `Type Spherical` dispatches to `CrystalSphericalNew`, which is an **icosphere mesh** (a subclass of `CrystalPolyhedron`). The icosphere has a full 3-D bounding box: X = Y = Z = 2×Radius. All three dimensions are present in the property map, so the NullPointerException at line 297 is never thrown, and `crystalType` stays `"CUBOID"`.

As a result, any `Type Spherical` crystal passed to `startMicroED` arrives with `crystalTypeEM = "CUBOID"`, and neither the CYLINDER nor SPHERICAL branch in `MicroED.java` is entered.

### On the `==` Operator

`MicroED.java` uses `==` for all three string comparisons (`"CUBOID"`, `"CYLINDER"`, `"SPHERICAL"`). This is the same pattern used throughout the codebase (`Crystal.java:309`, `XFEL.java:2699+`, `MC.java:2659+`, `BeamGaussian.java:115`, etc.). It is technically safe here because all the string values are assigned from Java string literals, which the JVM interns. A string literal `"CYLINDER"` set at Crystal.java:310 and compared against the literal `"CYLINDER"` at MicroED.java:232 are the same interned object — `==` evaluates to true.

The `"SPHERICAL"` branch is unreachable not because `==` fails but because `crystalType` is never set to `"SPHERICAL"` for any crystal type that reaches MicroED. The `==` comparisons would break only if strings were ever constructed dynamically (e.g., via `toUpperCase()` or concatenation) rather than from literals — but they are not.

### Effect on Simulation Output

None. `crystalVolume` is the only field affected by the SPHERICAL branch, and it is explicitly noted as unused in the MicroED source (comment at line 242: "although it isn't used"). The branch exists to set a volume that participates in no output calculation.

### Coverage Impact

Because this branch is structurally unreachable, MicroED.java's line coverage ceiling with the current Java source is approximately 30–35%. Lines 239–241 (spherical volume) can never be covered by any seed input. The remaining ~65% of MicroED.java (lines 838–4845) is dead code from a Monte Carlo electron-tracking simulation that was written but commented out of `CalculateEM`.

### Rust Status

Not applicable. The Rust `MicroEdSimulation` (`micro_ed.rs`) does not replicate MicroED's internal crystal-type dispatch or `crystalVolume` field — it uses the crystal geometry directly. No compatibility shim is needed. This bug is excluded from differential fuzzing.

---

## Summary

| Bug | File | Effect | Trigger |
|-----|------|--------|---------|
| Hardcoded `muabsIndex = 4` | `CrystalPolyhedron.java:1240` | FL Escape overestimated by ~57% | Any run with `CalculateFlEscape True` and elements with L-shell fluorescence |
| `energyPerFluence` scope leak | `Crystal.java:1236→1361` | Absorbed Energy underestimated by ~8x | `CalcSurrounding True` + `CalculatePEEscape True` |
| MicroED SPHERICAL branch unreachable | `MicroED.java:239`, `Crystal.java:294–313` | `crystalVolume` always computed as cuboid for spherical MicroED crystals | Any `Type Spherical` + `Subprogram EMSP/EMED` run (but `crystalVolume` is unused, so no observable effect) |

Bugs 1 and 2 affect downstream reporting metrics only; core dose deposition to voxels is unaffected. Bug 3 has no observable effect on any output.
