use crate::constants::KEV_TO_JOULES;
use crate::container::Container;
use crate::parser::config::{BeamConfig, Collimation};

use statrs::distribution::{ContinuousCDF, Normal};

/// Gaussian beam profile.
#[derive(Debug, Clone)]
pub struct BeamGaussian {
    /// FWHM X in µm.
    pub fwhm_x: f64,
    /// FWHM Y in µm.
    pub fwhm_y: f64,
    /// Flux in photons/sec.
    total_flux: f64,
    /// Photon energy in keV.
    photon_energy: f64,
    /// Pulse energy in mJ (optional).
    pulse_energy: f64,
    /// Energy FWHM for pink beam (optional).
    energy_fwhm: Option<f64>,
    /// Horizontal collimation in µm (None if uncollimated).
    coll_x_um: Option<f64>,
    /// Vertical collimation in µm (None if uncollimated).
    coll_y_um: Option<f64>,
    /// Whether circular collimation.
    is_circular: bool,
    /// Normalization factor for collimated beam.
    norm_factor: f64,
    /// Scale factor: converts intensity to J/µm²/s.
    scale_factor: f64,
    /// Attenuated photons per second (after container).
    attenuated_photons_per_sec: f64,
    /// Gaussian distributions.
    g_x: Normal,
    g_y: Normal,
}

impl BeamGaussian {
    /// FWHM to sigma conversion factor.
    const FWHM_TO_SIGMA: f64 = 2.354_820_045_030_949_4; // 2*sqrt(2*ln(2))

    pub fn from_config(config: &BeamConfig) -> Result<Self, String> {
        let fwhm_x = config.fwhm_x.ok_or("Gaussian beam requires FWHM X")?;
        let fwhm_y = config.fwhm_y.ok_or("Gaussian beam requires FWHM Y")?;
        let photons_per_sec = config.flux.ok_or("Gaussian beam requires flux")?;
        let photon_energy = config.energy.ok_or("Gaussian beam requires energy")?;

        let sigma_x = fwhm_x / Self::FWHM_TO_SIGMA;
        let sigma_y = fwhm_y / Self::FWHM_TO_SIGMA;

        let g_x = Normal::new(0.0, sigma_x).map_err(|e| e.to_string())?;
        let g_y = Normal::new(0.0, sigma_y).map_err(|e| e.to_string())?;

        let (coll_x_um, coll_y_um, is_circular) = match &config.collimation {
            Some(Collimation::Rectangular { h, v }) => (Some(*h), Some(*v), false),
            Some(Collimation::Circular { h, v }) => (Some(*h), Some(*v), true),
            Some(Collimation::Horizontal { h }) => (Some(*h), None, false),
            Some(Collimation::Vertical { v }) => (None, Some(*v), false),
            _ => (None, None, false),
        };

        // Calculate normalization factor
        let norm_factor = Self::calc_norm_factor(
            &g_x,
            &g_y,
            coll_x_um,
            coll_y_um,
            is_circular,
            sigma_x,
            sigma_y,
        );

        let scale_factor = photon_energy * KEV_TO_JOULES * photons_per_sec / norm_factor;

        Ok(BeamGaussian {
            fwhm_x,
            fwhm_y,
            total_flux: photons_per_sec,
            photon_energy,
            pulse_energy: config.pulse_energy.unwrap_or(0.0),
            energy_fwhm: config.energy_fwhm,
            coll_x_um,
            coll_y_um,
            is_circular,
            norm_factor,
            scale_factor,
            attenuated_photons_per_sec: photons_per_sec,
            g_x,
            g_y,
        })
    }

    fn calc_norm_factor(
        g_x: &Normal,
        g_y: &Normal,
        coll_x_um: Option<f64>,
        coll_y_um: Option<f64>,
        is_circular: bool,
        sigma_x: f64,
        sigma_y: f64,
    ) -> f64 {
        if is_circular {
            Self::bivariate_gaussian_volume_circular(sigma_x, sigma_y, coll_x_um, coll_y_um)
        } else {
            Self::bivariate_gaussian_volume_rect(g_x, g_y, coll_x_um, coll_y_um)
        }
    }

    fn bivariate_gaussian_volume_rect(
        g_x: &Normal,
        g_y: &Normal,
        coll_x: Option<f64>,
        coll_y: Option<f64>,
    ) -> f64 {
        // CDF-based integration for rectangular collimation
        let x_frac = match coll_x {
            Some(h) => g_x.cdf(h / 2.0) - g_x.cdf(-h / 2.0),
            None => 1.0,
        };
        let y_frac = match coll_y {
            Some(v) => g_y.cdf(v / 2.0) - g_y.cdf(-v / 2.0),
            None => 1.0,
        };
        x_frac * y_frac
    }

    fn bivariate_gaussian_volume_circular(
        sx: f64,
        sy: f64,
        coll_x: Option<f64>,
        coll_y: Option<f64>,
    ) -> f64 {
        // Match Java's polar integration for circular/elliptical aperture.
        // Java uses 100 angular steps with trapezoidal rule, integrating
        // the analytical radial integral of the 2D Gaussian over the ellipse.
        let half_x = coll_x.unwrap_or(100.0) / 2.0;
        let half_y = coll_y.unwrap_or(100.0) / 2.0;

        let a_const = 1.0 / (2.0 * std::f64::consts::PI * sx * sy);
        let step_theta = 2.0 * std::f64::consts::PI / 100.0;

        let mut last_held_height = 0.0;
        let mut last_held_theta = 0.0;
        let mut overall_sum = 0.0;

        let mut theta = 0.0;
        while theta <= 2.0 * std::f64::consts::PI {
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            let r =
                (half_x * half_y) / ((half_y * cos_t).powi(2) + (half_x * sin_t).powi(2)).sqrt();
            let cos_squared = cos_t * cos_t;
            let sin_squared = sin_t * sin_t;
            let cos_s_plus_sin_s_term =
                cos_squared / (2.0 * sx * sx) + sin_squared / (2.0 * sy * sy);
            let one_over = a_const / (2.0 * cos_s_plus_sin_s_term);
            let other_term = 1.0 - (-r * r * cos_s_plus_sin_s_term).exp();
            let overall_term = one_over * other_term;

            if theta != 0.0 {
                let area = (theta - last_held_theta) * ((last_held_height + overall_term) / 2.0);
                overall_sum += area;
            }
            last_held_theta = theta;
            last_held_height = overall_term;

            theta += step_theta;
        }

        overall_sum
    }

    fn gaussian_2d_value(g_x: &Normal, g_y: &Normal, x: f64, y: f64) -> f64 {
        use statrs::distribution::Continuous;
        g_x.pdf(x) * g_y.pdf(y)
    }

    fn recalc_scale_factor(&mut self) {
        self.scale_factor =
            self.photon_energy * KEV_TO_JOULES * self.attenuated_photons_per_sec / self.norm_factor;
    }

    fn is_in_collimation(&self, x: f64, y: f64) -> bool {
        if self.is_circular {
            let hx = self.coll_x_um.unwrap_or(f64::INFINITY) / 2.0;
            let hy = self.coll_y_um.unwrap_or(f64::INFINITY) / 2.0;
            (x / hx).powi(2) + (y / hy).powi(2) <= 1.0
        } else {
            let in_x = self.coll_x_um.is_none_or(|h| x.abs() <= h / 2.0);
            let in_y = self.coll_y_um.is_none_or(|v| y.abs() <= v / 2.0);
            in_x && in_y
        }
    }
}

impl super::Beam for BeamGaussian {
    fn beam_intensity(&self, coord_x: f64, coord_y: f64, off_axis_um: f64) -> f64 {
        let x = coord_x - off_axis_um;
        if !self.is_in_collimation(x, coord_y) {
            return 0.0;
        }
        Self::gaussian_2d_value(&self.g_x, &self.g_y, x, coord_y) * self.scale_factor
    }

    fn description(&self) -> String {
        let collimation = if let (Some(cx), Some(cy)) = (self.coll_x_um, self.coll_y_um) {
            format!("{:.1}x{:.1} um ", cx, cy)
        } else {
            String::new()
        };
        format!(
            "Gaussian beam, {}with {:.2} by {:.2} FWHM (x by y) and {:.1e} photons per second at {:.2} keV.",
            collimation, self.fwhm_x, self.fwhm_y, self.total_flux, self.photon_energy
        )
    }

    fn photons_per_sec(&self) -> f64 {
        self.attenuated_photons_per_sec
    }

    fn photon_energy(&self) -> f64 {
        self.photon_energy
    }

    fn pulse_energy(&self) -> f64 {
        self.pulse_energy
    }

    fn apply_container_attenuation(&mut self, container: &dyn Container) {
        let fraction = container.attenuation_fraction();
        self.attenuated_photons_per_sec = self.total_flux * (1.0 - fraction);
        if fraction > 0.0 {
            println!(
                "Beam photons per second after container attenuation is {:.2e} photons per second",
                self.attenuated_photons_per_sec
            );
        }
        self.recalc_scale_factor();
    }

    fn beam_minimum_dimension(&self) -> f64 {
        self.fwhm_x.min(self.fwhm_y)
    }

    fn beam_area(&self) -> f64 {
        match (&self.coll_x_um, &self.coll_y_um) {
            (Some(h), Some(v)) if self.is_circular => std::f64::consts::PI * (h / 2.0) * (v / 2.0),
            (Some(h), Some(v)) => h * v,
            _ => {
                // Use FWHM as effective area
                self.fwhm_x * self.fwhm_y
            }
        }
    }

    fn beam_x(&self) -> Option<f64> {
        self.coll_x_um.or(Some(self.fwhm_x))
    }

    fn beam_y(&self) -> Option<f64> {
        self.coll_y_um.or(Some(self.fwhm_y))
    }

    fn beam_type(&self) -> &str {
        "Gaussian"
    }

    fn is_circular(&self) -> bool {
        self.is_circular
    }

    fn energy_fwhm(&self) -> Option<f64> {
        self.energy_fwhm
    }

    fn sx(&self) -> f64 {
        self.fwhm_x / Self::FWHM_TO_SIGMA
    }

    fn sy(&self) -> f64 {
        self.fwhm_y / Self::FWHM_TO_SIGMA
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::beam::Beam;
    use crate::parser::config::{BeamConfig, BeamType, Collimation};

    fn make_gaussian() -> BeamGaussian {
        let config = BeamConfig {
            beam_type: Some(BeamType::Gaussian),
            energy: Some(12.0),
            flux: Some(1e10),
            fwhm_x: Some(20.0),
            fwhm_y: Some(20.0),
            collimation: Some(Collimation::Rectangular { h: 80.0, v: 80.0 }),
            ..Default::default()
        };
        BeamGaussian::from_config(&config).unwrap()
    }

    #[test]
    fn gaussian_peak_at_centre() {
        let b = make_gaussian();
        let center = b.beam_intensity(0.0, 0.0, 0.0);
        assert!(center > 0.0, "Gaussian center should be positive");

        // Off-center should be lower
        let off = b.beam_intensity(10.0, 0.0, 0.0);
        assert!(
            off < center,
            "Gaussian should fall off from center: center={}, off={}",
            center,
            off
        );
    }

    #[test]
    fn gaussian_falls_off_with_distance() {
        let b = make_gaussian();
        let mut prev = b.beam_intensity(0.0, 0.0, 0.0);
        for x in [5.0, 10.0, 15.0, 20.0] {
            let v = b.beam_intensity(x, 0.0, 0.0);
            assert!(v < prev, "Gaussian should decrease at x={}", x);
            prev = v;
        }
    }

    #[test]
    fn gaussian_circular_clips_correctly() {
        let config = BeamConfig {
            beam_type: Some(BeamType::Gaussian),
            energy: Some(12.0),
            flux: Some(1e10),
            fwhm_x: Some(20.0),
            fwhm_y: Some(20.0),
            collimation: Some(Collimation::Circular { h: 40.0, v: 40.0 }),
            ..Default::default()
        };
        let b = BeamGaussian::from_config(&config).unwrap();

        // Center should have intensity
        assert!(b.beam_intensity(0.0, 0.0, 0.0) > 0.0);
        // Corner of bounding box is outside circle
        assert_eq!(b.beam_intensity(19.9, 19.9, 0.0), 0.0);
        // On-axis within radius should work
        assert!(b.beam_intensity(19.0, 0.0, 0.0) > 0.0);
    }
}
