pub mod compute;
pub mod from_params;

pub use compute::CoefCalcCompute;
pub use from_params::CoefCalcFromParams;

use std::collections::HashSet;

/// Absorption/attenuation coefficient calculator trait.
///
/// Implements the calculation of X-ray absorption, attenuation, and scattering
/// coefficients for a crystal composition at a given photon energy.
pub trait CoefCalc: std::fmt::Debug + Send + Sync {
    /// Update coefficients for a given photon energy (keV).
    fn update_coefficients(&mut self, photon_energy: f64);

    /// Photoelectric absorption coefficient in µm⁻¹.
    fn absorption_coefficient(&self) -> f64;

    /// Total attenuation coefficient (photoelectric + inelastic + elastic) in µm⁻¹.
    fn attenuation_coefficient(&self) -> f64;

    /// Elastic (Rayleigh/coherent) scattering coefficient in µm⁻¹.
    fn elastic_coefficient(&self) -> f64;

    /// Inelastic (Compton) scattering coefficient in µm⁻¹.
    fn inelastic_coefficient(&self) -> f64;

    /// Crystal density in g/mL.
    fn density(&self) -> f64;

    /// Fluorescent escape factors for the current beam energy.
    /// Returns a 2D array: [element_index][factor_index].
    /// Factor layout per element (28 values):
    ///   [0]  mu_ratio (uj/upe)
    ///   [1]  K-shell ionization fraction
    ///   [2]  K fluorescence yield
    ///   [3]  K edge energy (keV)
    ///   [4]  K escape factor
    ///   [5..8]  L1 ionization, yield, energy, escape
    ///   [9..12]  L2 ...
    ///   [13..16] L3 ...
    ///   [17..27] M1-M5 ionization, binding energy pairs
    fn fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>>;

    /// Whether a cryo-solution surrounding is defined.
    fn is_cryo(&self) -> bool { false }

    /// Update cryo-solution coefficients.
    fn update_cryo_coefficients(&mut self, _photon_energy: f64) {}

    /// Cryo-solution absorption coefficient in µm⁻¹.
    fn cryo_absorption_coefficient(&self) -> f64 { 0.0 }

    /// Cryo-solution density in g/mL.
    fn cryo_density(&self) -> f64 { 0.0 }

    /// Cryo-solution inelastic coefficient.
    fn cryo_inelastic_coefficient(&self) -> f64 { 0.0 }

    /// Cryo fluorescent escape factors.
    fn cryo_fluorescent_escape_factors(&self, _beam_energy: f64) -> Vec<Vec<f64>> { vec![] }

    /// Present elements in crystal (or cryo solution).
    fn present_elements(&self, _cryo: bool) -> HashSet<String> { HashSet::new() }

    /// Solvent fraction.
    fn solvent_fraction(&self) -> f64 { 0.0 }
}

/// Create a CoefCalc from parsed crystal configuration.
pub fn create_coefcalc(
    config: &crate::parser::config::CrystalConfig,
) -> Result<Box<dyn CoefCalc>, String> {
    use crate::parser::config::CoefCalcType;

    let coefcalc_type = config.coefcalc.unwrap_or(CoefCalcType::Default);

    match coefcalc_type {
        CoefCalcType::Default | CoefCalcType::Average => {
            Ok(Box::new(CoefCalcFromParams::from_config(config)?))
        }
        other => Err(format!("CoefCalc type {:?} not yet implemented", other)),
    }
}
