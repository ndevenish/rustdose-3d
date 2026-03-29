---
name: Gumbel PDF negative beta issue in PE escape
description: Analysis of how negative Gumbel scale parameter causes PE dose to vanish in Rust port, and why Java works despite mathematically invalid distribution
type: project
---

## Gumbel PDF Beta Sign Issue

In `escape.rs`, the Gumbel-Left PDF originally guarded with `if beta <= 0.0 { return 0.0; }`. This caused 88% of photoelectric dose to silently vanish for low PE energy + high density cases (e.g. LiFePO₄ at 1.487 keV, density 3.60 g/mL).

**Why:** Beta is a quadratic in PE energy: `beta = scale[0]*E² + scale[1]*E`. The scale coefficients depend on density. Beta goes negative when `E_pe < -scale[1]/scale[0]`. For the lfp.txt case, critical energy is 0.817 keV but PE energy is only 0.788 keV.

**How to apply:** The fix was changing `beta <= 0.0` to `beta == 0.0`. Java has no guard at all — with negative beta the PDF produces negative values, but the negatives cancel in the normalization (ratio of negative/negative = positive) within `calc_propn_dose_deposited`. The polynomial parameterization was not designed for this low-energy regime, but the sign-cancellation accidentally preserves correct behavior. This is a known quirk to preserve when porting from Java.

**Conditions for negative beta:** high crystal density AND low PE kinetic energy (low beam energy or heavy atoms with high binding energy subtraction, e.g. Fe L-edges).


● Gumbel PDF Beta Sign Issue: Analysis Report

  The Bug

  In escape.rs:694, the Rust port guarded the Gumbel-Left PDF with:
  if beta <= 0.0 { return 0.0; }

  This is mathematically reasonable — a Gumbel distribution with beta <= 0 is undefined as a probability distribution. But Java's
  calculateGumbelDistribution has no such guard, and the code downstream relies on this.

  When Does Beta Go Negative?

  Beta is computed as a quadratic in PE energy:

  beta = scale[0] * E_pe² + scale[1] * E_pe

  Where the scale coefficients depend on crystal density. For peEnergy <= 20 keV:

  scale[0] = 0.0029 * ρ² - 0.0081 * ρ + 0.0076
  scale[1] = -0.0085 * ρ + 0.01751

  Beta goes negative when scale[0] * E_pe + scale[1] < 0, i.e.:

  E_pe < -scale[1] / scale[0]

  For this test case (LiFePO₄ crystal, ρ = 3.60 g/mL):
  - scale[0] = 0.0029(12.96) - 0.0081(3.60) + 0.0076 = 0.01602
  - scale[1] = -0.0085(3.60) + 0.01751 = -0.01309
  - Critical energy: 0.01309 / 0.01602 = 0.817 keV

  With beam energy 1.487 keV and binding energy subtraction 0.699 keV, the PE energy is 0.788 keV — below the critical threshold. So beta =
  -0.000372.

  This will occur whenever:
  - High density (scale[1] becomes more negative due to -0.0085ρ term)
  - Low PE energy (low beam energy and/or high binding energy subtraction, common with heavy atoms like Fe)

  The lfp.txt case (LiFePO₄ at 1.487 keV) is a natural example: the Fe L-edges consume most of the photon energy as binding energy, leaving a very
  low PE kinetic energy.

  What Java Actually Does With Negative Beta

  The Gumbel-Left PDF is f(x) = exp(z - exp(z)) / beta where z = (x - mu) / beta.

  With beta = -0.000372:
  - Near the mode (x ≈ mu): z = 0, f(mu) = exp(-1)/beta ≈ -989 (large negative)
  - Away from mode: z becomes very large positive (since beta is negative, x < mu gives positive z), exp(z) overflows, z - exp(z) = -∞, f(x) = 0

  So the "distribution" is a narrow spike of negative values near the mode, zero elsewhere. This is not a valid probability density.

  Why Java Gets Correct Results Anyway

  In calcProportionVoxDoseDepositedByDist, the path length "distribution" is used in ratios:

  propn[l] += (distanceWidths[i] * distanceHeights[i] / distnIntegral)
            * (distanceWidths[i] * energyHeight / totEnergyIntegral[i]);

  - distanceHeights[i] — derived from the Gumbel PDF: negative
  - distnIntegral — integral of the Gumbel PDF: negative
  - Ratio distanceHeights / distnIntegral = negative / negative = positive

  The negative signs cancel perfectly. The code is effectively computing |f(x)| / ∫|f(x)|dx — a valid normalized distribution. This works because:

  1. The PDF shape is correct (same bell curve), just flipped in sign
  2. The normalization divides out the sign
  3. The energy distribution terms are always positive

  Is Java "Wrong"?

  Strictly speaking, yes — the Gumbel distribution is undefined for beta ≤ 0, and the physical model (photoelectron path length distribution fitted
   from Monte Carlo simulations) was only calibrated for positive beta. The polynomial coefficients were derived from fits at specific
  density/energy ranges, and the extrapolation to very low PE energies at high densities produces an unphysical negative scale parameter.

  However, the code still produces physically reasonable results because:

  1. The distribution shape (which bins get what fraction of dose) is dominated by the energy deposition profile (calculateEnergyDistn), not the
  path length distribution alone.
  2. The convolution structure means the Gumbel only controls which path lengths are weighted more — and the sign-cancellation preserves this
  weighting correctly.
  3. At such low PE energies (< 1 keV), the photoelectron range is sub-micron anyway. The exact distribution shape matters less because most PE
  dose stays within 1-2 voxels of the source.

  The real issue is that the polynomial parameterization of the Gumbel scale was not designed for this regime. A more robust approach would be to
  clamp beta to a small positive value (e.g. max(beta, 0.001)) rather than allowing it to go negative. But since the Java code has been validated
  against experimental data in this configuration, matching its behavior (including the sign-cancellation) is the correct approach for the Rust
  port.

  Impact on This Test Case

  ┌───────────┬───────────┬───────────────────┬──────────────────────┐
  │  Metric   │   Java    │ Rust (with guard) │ Rust (guard removed) │
  ├───────────┼───────────┼───────────────────┼──────────────────────┤
  │ PE Escape │ 1.39e5    │ 0.00              │ 1.57e5               │
  ├───────────┼───────────┼───────────────────┼──────────────────────┤
  │ Max Dose  │ 26.53 MGy │ 3.02 MGy          │ 26.54 MGy            │
  ├───────────┼───────────┼───────────────────┼──────────────────────┤
  │ Avg Dose  │ 18.29 MGy │ 2.33 MGy          │ 17.76 MGy            │
  └───────────┴───────────┴───────────────────┴──────────────────────┘

  With the guard, 88% of the photoelectric dose vanished silently — it was removed from the source voxel (as PE dose) but never deposited anywhere
  (because propn was all zeros). Removing the guard restored correct behavior.

  Remaining Differences

  The ~3% differences between Java and Rust after the fix stem from the crystal geometry, not escape physics: Java maps Type Spherical to
  CrystalSphericalNew (42-vertex icosphere mesh with ray-casting depth), while Rust uses CrystalSpherical (analytic sphere with closed-form depth).
   This causes different voxel occupancy patterns, different depth calculations, and a different absorbed energy baseline (Java 2.41e-8 J vs Rust
  3.80e-8 J).
