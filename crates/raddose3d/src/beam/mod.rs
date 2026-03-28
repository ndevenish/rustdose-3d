pub mod experimental;
pub mod gaussian;
pub mod tophat;

pub use experimental::BeamExperimental;
pub use gaussian::BeamGaussian;
pub use tophat::BeamTophat;

use crate::container::Container;
use crate::parser::config::BeamConfig;

/// Beam trait for X-ray/electron beams.
pub trait Beam: std::fmt::Debug + Send + Sync {
    /// Returns the mean intensity at position (x, y) in J/µm²/s.
    fn beam_intensity(&self, coord_x: f64, coord_y: f64, off_axis_um: f64) -> f64;

    /// Returns a short description of the beam.
    fn description(&self) -> String;

    /// Returns flux in photons per second.
    fn photons_per_sec(&self) -> f64;

    /// Returns photon energy in keV.
    fn photon_energy(&self) -> f64;

    /// Returns pulse energy in mJ (0 if not applicable).
    fn pulse_energy(&self) -> f64 {
        0.0
    }

    /// Generate beam array (only used by experimental beam).
    fn generate_beam_array(&mut self) {}

    /// Apply container attenuation to the beam flux.
    fn apply_container_attenuation(&mut self, container: &dyn Container);

    /// Returns the minimum beam dimension in µm.
    fn beam_minimum_dimension(&self) -> f64;

    /// Returns beam area in µm².
    fn beam_area(&self) -> f64;

    /// Returns exposure in electrons/Å² (for electron beams).
    fn exposure(&self) -> f64 {
        0.0
    }

    /// Returns beam horizontal extent in µm.
    fn beam_x(&self) -> Option<f64>;

    /// Returns beam vertical extent in µm.
    fn beam_y(&self) -> Option<f64>;

    /// Returns beam type name.
    fn beam_type(&self) -> &str;

    /// Is the beam circularly collimated?
    fn is_circular(&self) -> bool {
        false
    }

    /// Returns semi-angle (for electron beams).
    fn semi_angle(&self) -> f64 {
        0.0
    }

    /// Returns aperture radius (for electron beams).
    fn aperture_radius(&self) -> f64 {
        0.0
    }

    /// Image X dimension.
    fn image_x(&self) -> f64 {
        0.0
    }

    /// Image Y dimension.
    fn image_y(&self) -> f64 {
        0.0
    }

    /// Energy FWHM for pink beam (None if monochromatic).
    fn energy_fwhm(&self) -> Option<f64> {
        None
    }

    /// Sigma X of the beam profile.
    fn sx(&self) -> f64;

    /// Sigma Y of the beam profile.
    fn sy(&self) -> f64;

    /// Returns beam size X (collimation) in µm (alias for beam_x).
    fn beam_size_x(&self) -> f64 {
        self.beam_x()
            .unwrap_or(self.sx() * 2.0 * (2.0 * std::f64::consts::LN_2).sqrt() * 3.0)
    }

    /// Returns beam size Y (collimation) in µm (alias for beam_y).
    fn beam_size_y(&self) -> f64 {
        self.beam_y()
            .unwrap_or(self.sy() * 2.0 * (2.0 * std::f64::consts::LN_2).sqrt() * 3.0)
    }

    /// Returns collimation type string ("Rectangular" or "Circular").
    fn collimation_type(&self) -> &str {
        if self.is_circular() {
            "Circular"
        } else {
            "Rectangular"
        }
    }

    /// Set photons per femtosecond (for XFEL timing).
    fn set_photons_per_fs(&mut self, _photons_per_fs: f64) {}
}

/// Create a Beam from parsed configuration.
pub fn create_beam(config: &BeamConfig) -> Result<Box<dyn Beam>, String> {
    let beam_type = config
        .beam_type
        .as_deref()
        .unwrap_or("gaussian")
        .to_lowercase();

    match beam_type.as_str() {
        "gaussian" => Ok(Box::new(BeamGaussian::from_config(config)?)),
        "tophat" => Ok(Box::new(BeamTophat::from_config(config)?)),
        "experimental" => Ok(Box::new(BeamExperimental::from_config(config)?)),
        other => Err(format!("Unknown beam type: {}", other)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::config::Collimation;

    const BEAM_ENERGY_MARKER: f64 = 12.533743110;

    fn default_beam_config(beam_type: &str) -> BeamConfig {
        BeamConfig {
            beam_type: Some(beam_type.to_string()),
            energy: Some(BEAM_ENERGY_MARKER),
            flux: Some(10.0),
            collimation: Some(Collimation::Rectangular { h: 1.0, v: 1.0 }),
            fwhm_x: Some(1.0),
            fwhm_y: Some(1.0),
            pixel_size_x: Some(1.0),
            pixel_size_y: Some(1.0),
            ..Default::default()
        }
    }

    #[test]
    fn create_beam_tophat() {
        let b = create_beam(&default_beam_config("tophat")).unwrap();
        assert_eq!(b.beam_type(), "TopHat");
        assert!((b.photon_energy() - BEAM_ENERGY_MARKER).abs() < 1e-10);
    }

    #[test]
    fn create_beam_gaussian() {
        let b = create_beam(&default_beam_config("gaussian")).unwrap();
        assert_eq!(b.beam_type(), "Gaussian");
        assert!((b.photon_energy() - BEAM_ENERGY_MARKER).abs() < 1e-10);
    }

    #[test]
    fn create_beam_fails_on_invalid_type() {
        let result = create_beam(&default_beam_config("invalid"));
        assert!(result.is_err());
    }

    #[test]
    fn create_beam_fails_on_empty_type() {
        let result = create_beam(&default_beam_config(""));
        assert!(result.is_err());
    }
}
