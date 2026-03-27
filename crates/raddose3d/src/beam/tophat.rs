use crate::constants::KEV_TO_JOULES;
use crate::container::Container;
use crate::parser::config::{BeamConfig, Collimation};

/// Top-hat (uniform) beam profile.
#[derive(Debug, Clone)]
pub struct BeamTophat {
    /// Horizontal extent in µm.
    beam_x_um: f64,
    /// Vertical extent in µm.
    beam_y_um: f64,
    /// Flux in photons/sec.
    total_flux: f64,
    /// Photon energy in keV.
    photon_energy: f64,
    /// Pulse energy in mJ (optional).
    pulse_energy: f64,
    /// Energy FWHM for pink beam (optional).
    energy_fwhm: Option<f64>,
    /// Beam exposure in electrons/Å².
    exposure_val: f64,
    /// Semi-angle for electron beams.
    semi_angle: f64,
    /// Aperture radius for electron beams.
    aperture_radius: f64,
    /// Image dimensions.
    image_x: f64,
    image_y: f64,
    /// Whether circular collimation.
    is_circular: bool,
    /// Attenuated flux.
    attenuated_photons_per_sec: f64,
}

impl BeamTophat {
    pub fn from_config(config: &BeamConfig) -> Result<Self, String> {
        let photon_energy = config.energy.ok_or("TopHat beam requires energy")?;

        let (beam_x_um, beam_y_um, is_circular) = match &config.collimation {
            Some(Collimation::Rectangular { h, v }) => (*h, *v, false),
            Some(Collimation::Circular { h, v }) => (*h, *v, true),
            _ => return Err("TopHat beam requires collimation dimensions".to_string()),
        };

        let photons_per_sec = config.flux.unwrap_or(0.0);

        Ok(BeamTophat {
            beam_x_um,
            beam_y_um,
            total_flux: photons_per_sec,
            photon_energy,
            pulse_energy: config.pulse_energy.unwrap_or(0.0),
            energy_fwhm: config.energy_fwhm,
            exposure_val: config.exposure.unwrap_or(0.0),
            semi_angle: config.semi_angle.unwrap_or(0.0),
            aperture_radius: config.aperture_radius.unwrap_or(0.0),
            image_x: config.image_x.unwrap_or(0.0),
            image_y: config.image_y.unwrap_or(0.0),
            is_circular,
            attenuated_photons_per_sec: photons_per_sec,
        })
    }

    fn is_in_beam(&self, x: f64, y: f64) -> bool {
        if self.is_circular {
            let hx = self.beam_x_um / 2.0;
            let hy = self.beam_y_um / 2.0;
            (x / hx).powi(2) + (y / hy).powi(2) <= 1.0
        } else {
            x.abs() <= self.beam_x_um / 2.0 && y.abs() <= self.beam_y_um / 2.0
        }
    }

    fn uniform_intensity(&self) -> f64 {
        let area = self.beam_area_internal();
        if area <= 0.0 {
            return 0.0;
        }
        KEV_TO_JOULES * self.attenuated_photons_per_sec * self.photon_energy / area
    }

    fn beam_area_internal(&self) -> f64 {
        if self.is_circular {
            std::f64::consts::PI * (self.beam_x_um / 2.0) * (self.beam_y_um / 2.0)
        } else {
            self.beam_x_um * self.beam_y_um
        }
    }
}

impl super::Beam for BeamTophat {
    fn beam_intensity(&self, coord_x: f64, coord_y: f64, _off_axis_um: f64) -> f64 {
        if self.is_in_beam(coord_x, coord_y) {
            self.uniform_intensity()
        } else {
            0.0
        }
    }

    fn description(&self) -> String {
        format!(
            "Top hat beam, {:.1} by {:.1} (x by y) um with {:.2e} photons per second at {:.2}keV.",
            self.beam_x_um, self.beam_y_um, self.total_flux, self.photon_energy
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
    }

    fn beam_minimum_dimension(&self) -> f64 {
        self.beam_x_um.min(self.beam_y_um)
    }

    fn beam_area(&self) -> f64 {
        self.beam_area_internal()
    }

    fn beam_x(&self) -> Option<f64> {
        Some(self.beam_x_um)
    }

    fn beam_y(&self) -> Option<f64> {
        Some(self.beam_y_um)
    }

    fn beam_type(&self) -> &str {
        "TopHat"
    }

    fn is_circular(&self) -> bool {
        self.is_circular
    }

    fn exposure(&self) -> f64 {
        self.exposure_val
    }

    fn semi_angle(&self) -> f64 {
        self.semi_angle
    }

    fn aperture_radius(&self) -> f64 {
        self.aperture_radius
    }

    fn image_x(&self) -> f64 {
        self.image_x
    }

    fn image_y(&self) -> f64 {
        self.image_y
    }

    fn energy_fwhm(&self) -> Option<f64> {
        self.energy_fwhm
    }

    fn sx(&self) -> f64 {
        // Uniform distribution: sigma = width / sqrt(12)
        self.beam_x_um / 12.0_f64.sqrt()
    }

    fn sy(&self) -> f64 {
        self.beam_y_um / 12.0_f64.sqrt()
    }
}
