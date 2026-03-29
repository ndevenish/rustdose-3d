use crate::parser::config::CrystalConfig;

use super::compute::CoefCalcCompute;
use super::CoefCalc;

/// CoefCalcSmallMolecules: coefficient calculator for small-molecule crystals.
///
/// Unlike protein crystals, small-molecule crystals do NOT fill the remaining
/// unit cell volume with water. Atoms are specified per monomer.
#[derive(Debug)]
pub struct CoefCalcSmallMolecules {
    pub compute: CoefCalcCompute,
}

impl CoefCalcSmallMolecules {
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let mut compute = CoefCalcCompute::new();

        let cell_a = config.cell_a.ok_or("SmallMolecules requires cell a")?;
        let cell_b = config.cell_b.ok_or("SmallMolecules requires cell b")?;
        let cell_c = config.cell_c.ok_or("SmallMolecules requires cell c")?;
        let cell_alpha = config.cell_alpha.unwrap_or(90.0);
        let cell_beta = config.cell_beta.unwrap_or(90.0);
        let cell_gamma = config.cell_gamma.unwrap_or(90.0);

        compute.calculate_cell_volume(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);

        let num_monomers = config.num_monomers.unwrap_or(1).max(1);
        compute.num_monomers = num_monomers;

        // Add small molecule atoms per monomer × num_monomers
        for ec in &config.small_mole_atoms {
            compute.increment_macro(&ec.symbol.to_uppercase(), ec.count * num_monomers as f64);
        }

        // Heavy solvent concentrations
        let solv_names: Vec<String> = config
            .heavy_solution_conc
            .iter()
            .map(|e| e.symbol.clone())
            .collect();
        let solv_concs: Vec<f64> = config.heavy_solution_conc.iter().map(|e| e.count).collect();
        if !solv_names.is_empty() {
            compute.add_solvent_concentrations(&solv_names, &solv_concs);
        }

        // Solvent fraction (if given, add solvent water; otherwise no water fill)
        let sf = config.solvent_fraction.unwrap_or(-1.0);
        if sf > 0.0 {
            compute.sol_fraction = sf;
            compute.calculate_solvent_water(sf);
        }
        // Note: fillRestWithWater = false for small molecules

        // Cryo/surrounding setup
        let calc_surrounding = config
            .calc_surrounding
            .as_deref()
            .is_some_and(|s| s.eq_ignore_ascii_case("true"));

        if calc_surrounding {
            let cryo_solution_names: Vec<String> = config
                .surrounding_heavy_conc
                .iter()
                .map(|e| e.symbol.clone())
                .collect();
            let cryo_solution_concs: Vec<f64> = config
                .surrounding_heavy_conc
                .iter()
                .map(|e| e.count)
                .collect();

            let oil_element_names: Vec<String> = config
                .oil_elements
                .iter()
                .map(|e| e.symbol.clone())
                .collect();
            let oil_element_nums: Vec<f64> = config.oil_elements.iter().map(|e| e.count).collect();

            let oil_density = config.oil_density.unwrap_or(0.0);

            compute.add_cryo_concentrations(
                &cryo_solution_names,
                &cryo_solution_concs,
                config.oil_based.as_deref(),
                &oil_element_names,
                &oil_element_nums,
                oil_density,
            );
            compute.calculate_cryo_density();
        }

        compute.calculate_density();
        Ok(CoefCalcSmallMolecules { compute })
    }
}

impl CoefCalc for CoefCalcSmallMolecules {
    fn update_coefficients(&mut self, photon_energy: f64) {
        let (photo, coherent, compton, total) =
            self.compute.calculate_coefficients_all(photon_energy);
        self.compute.abs_coeff_photo = photo;
        self.compute.elas_coeff = coherent;
        self.compute.abs_coeff_comp = compton;
        self.compute.att_coeff = total;
        let (_p, coh_m, _c, _t) = self.compute.calculate_coefficients_macro(photon_energy);
        self.compute.elas_coeff_macro = coh_m;
    }

    fn absorption_coefficient(&self) -> f64 {
        self.compute.abs_coeff_photo
    }
    fn attenuation_coefficient(&self) -> f64 {
        self.compute.att_coeff
    }
    fn elastic_coefficient(&self) -> f64 {
        self.compute.elas_coeff
    }
    fn inelastic_coefficient(&self) -> f64 {
        self.compute.abs_coeff_comp
    }
    fn density(&self) -> f64 {
        self.compute.crystal_density
    }
    fn fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>> {
        self.compute.calc_fluorescent_escape_factors(beam_energy)
    }
    fn solvent_fraction(&self) -> f64 {
        self.compute.sol_fraction
    }

    fn is_cryo(&self) -> bool {
        self.compute.is_cryo()
    }

    fn update_cryo_coefficients(&mut self, photon_energy: f64) {
        self.compute.update_cryo_coefficients(photon_energy);
    }

    fn cryo_absorption_coefficient(&self) -> f64 {
        self.compute.cryo_abs_coeff_photo
    }

    fn cryo_density(&self) -> f64 {
        self.compute.cryo_density
    }

    fn cryo_inelastic_coefficient(&self) -> f64 {
        self.compute.cryo_abs_coeff_comp
    }

    fn cryo_fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>> {
        self.compute
            .calc_cryo_fluorescent_escape_factors(beam_energy)
    }

    fn present_elements(&self, cryo: bool) -> std::collections::HashSet<String> {
        if cryo {
            self.compute.cryo_elements.clone()
        } else {
            self.compute.present_elements.clone()
        }
    }
}
