/// Normal energy distribution for pink beam modelling.
///
/// Ported from `NormalEnergyDistribution` and `SampleNormalEnergyDistribution`
/// in the Java source.
use statrs::distribution::{Continuous, ContinuousCDF, Normal};

/// Truncation bound: ±40σ from mean (same as Java).
const TRUNCATION_SIGMAS: f64 = 40.0;

/// Normal energy distribution for a pink beam.
pub struct NormalEnergyDistribution {
    dist: Normal,
    pub mean_energy: f64,
    pub energy_fwhm: f64,
    pub sigma: f64,
}

impl NormalEnergyDistribution {
    /// Create a new distribution with `mean` (keV) and `fwhm` (keV).
    pub fn new(mean: f64, fwhm: f64) -> Self {
        let sigma = fwhm_to_sigma(fwhm);
        let dist = Normal::new(mean, sigma).expect("Invalid Normal distribution parameters");
        NormalEnergyDistribution {
            dist,
            mean_energy: mean,
            energy_fwhm: fwhm,
            sigma,
        }
    }

    /// PDF at `value`.
    pub fn pdf(&self, value: f64) -> f64 {
        self.dist.pdf(value)
    }

    /// CDF at `value`.
    pub fn cdf(&self, value: f64) -> f64 {
        self.dist.cdf(value)
    }

    /// Inverse CDF (quantile function) at `p` ∈ (0, 1).
    pub fn inverse_cdf(&self, p: f64) -> f64 {
        self.dist.inverse_cdf(p)
    }

    /// Returns `true` if `value` is outside the truncation bounds.
    pub fn is_neglected(&self, value: f64) -> bool {
        let cut = TRUNCATION_SIGMAS * self.sigma;
        value < self.mean_energy - cut || value > self.mean_energy + cut
    }

    /// Lower truncation bound.
    pub fn lower_bound(&self) -> f64 {
        self.mean_energy - TRUNCATION_SIGMAS * self.sigma
    }

    /// Upper truncation bound.
    pub fn upper_bound(&self) -> f64 {
        self.mean_energy + TRUNCATION_SIGMAS * self.sigma
    }
}

/// Sample `n` energies from a truncated normal distribution.
///
/// Returns a sorted array of `n` energy values uniformly spaced across the
/// CDF of the distribution, ignoring tails beyond ±40σ.
///
/// Matches the Java `SampleNormalEnergyDistribution` algorithm.
pub fn sample_normal_energies(mean: f64, fwhm: f64, n: usize) -> Vec<f64> {
    let dist = NormalEnergyDistribution::new(mean, fwhm);

    // Determine how many extra CDF points we need to have exactly `n` after
    // discarding the neglected tails.
    let mut before = n;
    loop {
        let count_kept = (0..before)
            .filter(|&i| {
                let p = i as f64 / (before - 1) as f64;
                let e = dist.inverse_cdf(p);
                !dist.is_neglected(e)
            })
            .count();
        if count_kept == n {
            break;
        }
        before += 1;
    }
    before -= 1; // Java does this adjustment

    let mut energies = Vec::with_capacity(n);
    for i in 0..before {
        let p = i as f64 / (before - 1) as f64;
        let e = dist.inverse_cdf(p);
        if !dist.is_neglected(e) {
            energies.push(e);
        }
    }
    energies
}

fn fwhm_to_sigma(fwhm: f64) -> f64 {
    // σ = FWHM / (2 √(2 ln 2))  ≈  FWHM / 2.3548
    fwhm / (2.0 * (2.0_f64.ln() * 2.0).sqrt())
}
