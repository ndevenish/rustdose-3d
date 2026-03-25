pub mod average;
pub mod compute;
pub mod from_cif;
pub mod from_params;
pub mod from_pdb;
pub mod from_sequence;
pub mod from_sequence_saxs;
pub mod micro_ed;
pub mod saxs;
pub mod small_molecules;

pub use average::CoefCalcAverage;
pub use compute::CoefCalcCompute;
pub use from_cif::CoefCalcFromCIF;
pub use from_params::CoefCalcFromParams;
pub use from_pdb::CoefCalcFromPDB;
pub use from_sequence::CoefCalcFromSequence;
pub use from_sequence_saxs::CoefCalcFromSequenceSAXS;
pub use micro_ed::CoefCalcMicroED;
pub use saxs::CoefCalcSAXS;
pub use small_molecules::CoefCalcSmallMolecules;

use std::collections::HashSet;

/// Absorption/attenuation coefficient calculator trait.
///
/// Implements the calculation of X-ray absorption, attenuation, and scattering
/// coefficients for a crystal composition at a given photon energy.
pub trait CoefCalc: std::fmt::Debug + Send + Sync {
    /// Human-readable description of the coefficient source and values.
    /// Matches Java's `CoefCalcCompute.toString()`.
    fn description(&self) -> String {
        format!(
            "\nCrystal coefficients calculated with RADDOSE-3D. \n\
             Photelectric Coefficient: {:.2e} /um.\n\
             Inelastic Coefficient: {:.2e} /um.\n\
             Elastic Coefficient: {:.2e} /um.\n\
             Attenuation Coefficient: {:.2e} /um.\n\
             Density: {:.2} g/ml.\n",
            self.absorption_coefficient(),
            self.inelastic_coefficient(),
            self.elastic_coefficient(),
            self.attenuation_coefficient(),
            self.density(),
        )
    }

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
    fn fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>>;

    /// Whether a cryo-solution surrounding is defined.
    fn is_cryo(&self) -> bool {
        false
    }

    /// Update cryo-solution coefficients.
    fn update_cryo_coefficients(&mut self, _photon_energy: f64) {}

    /// Cryo-solution absorption coefficient in µm⁻¹.
    fn cryo_absorption_coefficient(&self) -> f64 {
        0.0
    }

    /// Cryo-solution density in g/mL.
    fn cryo_density(&self) -> f64 {
        0.0
    }

    /// Cryo-solution inelastic coefficient.
    fn cryo_inelastic_coefficient(&self) -> f64 {
        0.0
    }

    /// Cryo fluorescent escape factors.
    fn cryo_fluorescent_escape_factors(&self, _beam_energy: f64) -> Vec<Vec<f64>> {
        vec![]
    }

    /// Present elements in crystal (or cryo solution).
    fn present_elements(&self, _cryo: bool) -> HashSet<String> {
        HashSet::new()
    }

    /// Solvent fraction.
    fn solvent_fraction(&self) -> f64 {
        0.0
    }
}

/// Create a CoefCalc from parsed crystal configuration.
pub fn create_coefcalc(
    config: &crate::parser::config::CrystalConfig,
) -> Result<Box<dyn CoefCalc>, String> {
    use crate::parser::config::CoefCalcType;

    let coefcalc_type = config.coefcalc.unwrap_or(CoefCalcType::Default);

    match coefcalc_type {
        CoefCalcType::Default => Ok(Box::new(CoefCalcFromParams::from_config(config)?)),
        CoefCalcType::Average => Ok(Box::new(CoefCalcAverage)),
        CoefCalcType::Cif => {
            let path = config.cif.as_deref().ok_or("CIF mode requires cif path")?;
            Ok(Box::new(CoefCalcFromCIF::from_file(path)?))
        }
        CoefCalcType::SmallMole => Ok(Box::new(CoefCalcSmallMolecules::from_config(config)?)),
        CoefCalcType::Sequence => Ok(Box::new(CoefCalcFromSequence::from_config(config)?)),
        CoefCalcType::Saxs => Ok(Box::new(CoefCalcSAXS::from_config(config)?)),
        CoefCalcType::SaxsSeq => Ok(Box::new(CoefCalcFromSequenceSAXS::from_config(config)?)),
        CoefCalcType::Pdb => Ok(Box::new(CoefCalcFromPDB::from_config(config)?)),
        // MicroED CoefCalcType not yet in the parser; map via RdFortran placeholder or leave:
        CoefCalcType::RdFortran => {
            // Legacy RD v2 subprocess not implemented; fall back to FromParams
            eprintln!("Warning: RDFortran/RDv2 CoefCalc not implemented, using Default mode.");
            Ok(Box::new(CoefCalcFromParams::from_config(config)?))
        }
    }
}
