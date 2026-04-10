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

    // -----------------------------------------------------------------------
    // MC / XFEL / MicroED extensions — default stubs
    // -----------------------------------------------------------------------

    /// Per-element photoelectric probabilities at the given energy.
    fn photo_electric_probs_element(&self, _energy: f64) -> std::collections::HashMap<String, f64> {
        std::collections::HashMap::new()
    }

    /// Per-element Compton probabilities at the given energy.
    fn compton_probs_element(&self, _energy: f64) -> std::collections::HashMap<String, f64> {
        std::collections::HashMap::new()
    }

    /// Per-element photoelectric probabilities for the cryo surrounding.
    fn photo_electric_probs_element_surrounding(
        &self,
        _energy: f64,
    ) -> std::collections::HashMap<String, f64> {
        std::collections::HashMap::new()
    }

    /// Stopping power (keV/nm) at the given energy.
    fn stopping_power(&self, _energy: f64, _cryo: bool) -> f64 {
        0.0
    }

    /// Populate cross-section coefficients (called once before MC loop).
    fn populate_cross_section_coefficients(&mut self) {}

    /// Electron elastic mean free path (nm).
    /// Takes `&mut self` because it populates per-element cross-section state
    /// as a side effect (read by `elastic_probs()`).
    fn electron_elastic_mfpl(&mut self, _energy: f64, _cryo: bool) -> f64 {
        0.0
    }

    /// Elastic scattering probabilities keyed by element symbol.
    fn elastic_probs(&self, _cryo: bool) -> std::collections::HashMap<String, f64> {
        std::collections::HashMap::new()
    }

    /// Bethe ionisation cross-section (nm) for inner shells.
    fn bethe_ionisation_x_section(&self, _energy: f64, _cryo: bool) -> f64 {
        0.0
    }

    /// GOS inelastic mean free path (nm).
    fn gos_inel(&self, _cryo: bool, _energy: f64) -> f64 {
        0.0
    }

    /// GOS inner-shell lambda (nm).
    fn gos_inner_lambda(&self, _cryo: bool) -> f64 {
        0.0
    }

    /// GOS outer-shell lambda (nm).
    fn gos_outer_lambda(&self, _cryo: bool) -> f64 {
        0.0
    }

    /// GOS shell ionisation probabilities.
    fn gos_shell_probs(
        &self,
        _cryo: bool,
        _lambda: f64,
    ) -> std::collections::HashMap<String, Vec<f64>> {
        std::collections::HashMap::new()
    }

    /// Plasma (plasmon) mean free path (nm).
    fn plasma_mfpl(&self, _energy: f64, _cryo: bool) -> f64 {
        0.0
    }

    /// Average inelastic energy loss (keV).
    fn avg_inelastic_energy(&self, _energy: f64) -> f64 {
        0.0
    }

    /// Plasma frequency energy (eV).
    fn plasma_energy(&self, _cryo: bool) -> f64 {
        0.0
    }

    /// FSE (fast secondary electron) lambda given cross section.
    fn fse_lambda(&self, _cross_section: f64, _cryo: bool) -> f64 {
        0.0
    }

    /// All inner-shell ionisation probabilities.
    fn all_shell_probs(&self, _cryo: bool) -> std::collections::HashMap<String, Vec<f64>> {
        std::collections::HashMap::new()
    }

    /// GOS outer-shell probabilities.
    fn gos_outer_shell_probs(
        &self,
        _cryo: bool,
        _lambda: f64,
    ) -> std::collections::HashMap<String, f64> {
        std::collections::HashMap::new()
    }

    /// Number of valence electrons per subshell for an element.
    fn num_valence_electrons_subshells(&self, _element: &str) -> Vec<i32> {
        Vec::new()
    }

    /// GOS variable (flat per-element probability array, indexed by shell offset).
    fn gos_variable(&self, _cryo: bool) -> std::collections::HashMap<String, Vec<f64>> {
        std::collections::HashMap::new()
    }

    /// Plasmon variable (single probability value).
    fn plasmon_variable(&self, _cryo: bool) -> f64 {
        0.0
    }

    /// Return the energy adjustment parameter.
    fn return_adjustment(&self) -> f64 {
        0.0
    }

    /// Wk for a molecule shell.
    fn wk_molecule(&self, _a: f64, _element: &str, _shell: usize, _cryo: bool) -> f64 {
        0.0
    }

    /// Wcb for all outer shells.
    fn wcb_all(&self, _cryo: bool) -> f64 {
        0.0
    }

    /// Number of simulated electrons (for MC photon count).
    fn number_simulated_electrons(&self) -> u64 {
        0
    }

    /// Total number of atoms in the crystal given volume (cm³).
    fn total_atoms_in_crystal(&self, _volume: f64) -> f64 {
        0.0
    }

    /// Total atoms of a specific element given volume (cm³).
    fn total_atoms_in_crystal_element(&self, _volume: f64, _element: &str) -> f64 {
        0.0
    }

    /// Relative shell probabilities given per-element absolute probs.
    fn get_relative_shell_probs_for_element(
        &self,
        _element_abs_probs: &std::collections::HashMap<String, f64>,
        _energy: f64,
    ) -> std::collections::HashMap<String, Vec<f64>> {
        std::collections::HashMap::new()
    }

    /// Recoil energy for a distant GOS collision (keV).
    fn recoil_energy_distant(&self, _energy: f64, _wak: f64, _qak: f64) -> f64 {
        0.0
    }

    /// Per-element relative ionisation probabilities for each shell.
    fn relative_shell_probs(
        &self,
        _energy: f64,
        _cryo: bool,
    ) -> std::collections::HashMap<String, Vec<f64>> {
        std::collections::HashMap::new()
    }

    /// Binding energy (keV) of a shell for a given element.
    fn shell_binding_energy(&self, _element: &str, _shell: usize) -> f64 {
        0.0
    }

    /// Atomic number of a given element symbol.
    fn atomic_number_of(&self, _element: &str) -> usize {
        0
    }

    /// Fluorescence yield for a given element shell.
    fn shell_fluorescence_yield(&self, _element: &str, _shell: usize) -> f64 {
        0.0
    }

    /// GOS outer-shell ionisation probabilities (simple form) keyed by element.
    fn gos_outer_shell_probs_simple(
        &self,
        _cryo: bool,
        _lambda: f64,
    ) -> std::collections::HashMap<String, f64> {
        std::collections::HashMap::new()
    }

    /// Elastic coefficient for macro-molecular contribution.
    fn elastic_coefficient_macro(&self) -> f64 {
        self.elastic_coefficient()
    }

    /// Set photons per femtosecond (used by RD3D check).
    fn set_photons_per_fs(&mut self, _photons_per_fs: f64) {}

    /// Get density in g/mL (alias for density()).
    fn get_density(&self) -> f64 {
        self.density()
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
        CoefCalcType::MicroED => Ok(Box::new(CoefCalcMicroED::from_config(config)?)),
        // Legacy RD v2 subprocess not implemented; fall back to FromParams
        CoefCalcType::RdFortran => {
            println!("Warning: RDFortran/RDv2 CoefCalc not implemented, using Default mode.");
            Ok(Box::new(CoefCalcFromParams::from_config(config)?))
        }
    }
}
