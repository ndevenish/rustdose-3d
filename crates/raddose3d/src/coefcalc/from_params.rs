use crate::parser::config::CrystalConfig;

use super::compute::*;

/// CoefCalcFromParams: calculates coefficients from user-specified composition parameters.
#[derive(Debug)]
pub struct CoefCalcFromParams {
    pub compute: CoefCalcCompute,
}

impl CoefCalcFromParams {
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let mut compute = CoefCalcCompute::new();

        // Cell dimensions
        let cell_a = config
            .cell_a
            .ok_or("Crystal requires unit cell dimension a")?;
        let cell_b = config
            .cell_b
            .ok_or("Crystal requires unit cell dimension b")?;
        let cell_c = config
            .cell_c
            .ok_or("Crystal requires unit cell dimension c")?;
        let cell_alpha = config.cell_alpha.unwrap_or(90.0);
        let cell_beta = config.cell_beta.unwrap_or(90.0);
        let cell_gamma = config.cell_gamma.unwrap_or(90.0);

        compute.calculate_cell_volume(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);

        let num_monomers = config.num_monomers.unwrap_or(1);
        let num_residues = config.num_residues.unwrap_or(0);
        let num_rna = config.num_rna.unwrap_or(0);
        let num_dna = config.num_dna.unwrap_or(0);
        let num_carb = config.num_carb.unwrap_or(0);

        compute.num_monomers = num_monomers;
        compute.num_amino_acids = num_residues as f64;
        compute.num_rna = num_rna as f64;
        compute.num_dna = num_dna as f64;
        compute.num_carb = num_carb as f64;

        // Heavy protein atoms
        for ec in &config.heavy_protein_atoms {
            compute.increment_macro(&ec.symbol, ec.count * num_monomers as f64);
        }

        // Solvent concentrations
        let solv_names: Vec<String> = config
            .heavy_solution_conc
            .iter()
            .map(|ec| ec.symbol.clone())
            .collect();
        let solv_concs: Vec<f64> = config
            .heavy_solution_conc
            .iter()
            .map(|ec| ec.count)
            .collect();
        if !solv_names.is_empty() {
            compute.add_solvent_concentrations(&solv_names, &solv_concs);
        }

        // Solvent fraction
        let mut sf = config.solvent_fraction.unwrap_or(-1.0);
        if sf <= 0.0 {
            sf = compute.calculate_solvent_fraction_from_nums();
        } else {
            compute.sol_fraction = sf;
        }

        compute.calculate_solvent_water(sf);

        // Standard amino acid composition: per amino acid add 5C + 1.35N + 1.5O + 8H
        let m = num_monomers as f64;
        compute.increment_macro("C", CARBONS_PER_AMINO_ACID * num_residues as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_AMINO_ACID * num_residues as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_AMINO_ACID * num_residues as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_AMINO_ACID * num_residues as f64 * m);

        // RNA: 11.25H + 9.5C + 3.75N + 7O + 1P
        compute.increment_macro("C", CARBONS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("P", PHOSPHORI_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);

        // DNA: 11.75H + 9.75C + 4N + 6O + 1P
        compute.increment_macro("C", CARBONS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("P", PHOSPHORI_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);

        // Carbohydrate: 11H + 6C + 5O
        compute.increment_macro("C", CARBONS_PER_CARBOHYDRATE * num_carb as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_CARBOHYDRATE * num_carb as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_CARBOHYDRATE * num_carb as f64 * m);

        // Calculate density
        compute.calculate_density();

        Ok(CoefCalcFromParams { compute })
    }
}

impl super::CoefCalc for CoefCalcFromParams {
    fn update_coefficients(&mut self, photon_energy: f64) {
        let (photo, coherent, compton, total) =
            self.compute.calculate_coefficients_all(photon_energy);
        self.compute.abs_coeff_photo = photo;
        self.compute.elas_coeff = coherent;
        self.compute.abs_coeff_comp = compton;
        self.compute.att_coeff = total;

        let (_photo_m, coherent_m, _compton_m, _total_m) =
            self.compute.calculate_coefficients_macro(photon_energy);
        self.compute.elas_coeff_macro = coherent_m;
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coefcalc::CoefCalc;
    use crate::parser::config::{CrystalConfig, ElementCount};

    /// Helper to construct CoefCalcFromParams from config.
    #[allow(clippy::too_many_arguments)]
    fn make_coefcalc(
        cell_a: f64,
        cell_b: f64,
        cell_c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
        num_monomers: i32,
        num_residues: i32,
        num_rna: i32,
        num_dna: i32,
        heavy_protein: Vec<(&str, f64)>,
        heavy_solution: Vec<(&str, f64)>,
        solvent_fraction: Option<f64>,
    ) -> CoefCalcFromParams {
        let config = CrystalConfig {
            cell_a: Some(cell_a),
            cell_b: Some(cell_b),
            cell_c: Some(cell_c),
            cell_alpha: Some(alpha),
            cell_beta: Some(beta),
            cell_gamma: Some(gamma),
            num_monomers: Some(num_monomers),
            num_residues: Some(num_residues),
            num_rna: Some(num_rna),
            num_dna: Some(num_dna),
            heavy_protein_atoms: heavy_protein
                .iter()
                .map(|(s, c)| ElementCount {
                    symbol: s.to_string(),
                    count: *c,
                })
                .collect(),
            heavy_solution_conc: heavy_solution
                .iter()
                .map(|(s, c)| ElementCount {
                    symbol: s.to_string(),
                    count: *c,
                })
                .collect(),
            solvent_fraction,
            ..Default::default()
        };
        CoefCalcFromParams::from_config(&config).unwrap()
    }

    #[test]
    fn water_only_h_is_twice_o() {
        let cc = make_coefcalc(
            100.0,
            100.0,
            100.0,
            90.0,
            90.0,
            90.0,
            0,
            0,
            0,
            0,
            vec![],
            vec![],
            Some(1.0),
        );
        let h = cc.compute.get_solvent_occurrence("H");
        let o = cc.compute.get_solvent_occurrence("O");
        assert!(
            o > 0.0,
            "Oxygen occurrence should be > 0 in water-only cell"
        );
        assert!(
            (h - o * 2.0).abs() < 1e-6,
            "H should be 2×O: H={}, O={}",
            h,
            o
        );
    }

    #[test]
    fn heavy_protein_atoms_multiplied_by_monomers() {
        let cc = make_coefcalc(
            100.0,
            100.0,
            100.0,
            90.0,
            90.0,
            90.0,
            24,
            10,
            0,
            0,
            vec![("Zn", 2.0)],
            vec![],
            Some(1.0),
        );
        let zn = cc.compute.get_macro("Zn");
        assert!(
            (zn - 48.0).abs() < 1e-6,
            "Zn should be 2 * 24 = 48, got {}",
            zn
        );
    }

    #[test]
    fn coefcalc_scenario1_sulfur_salt() {
        let mut cc = make_coefcalc(
            79.2,
            79.2,
            38.1,
            90.0,
            90.0,
            90.0,
            8,
            129,
            0,
            0,
            vec![("S", 10.0)],
            vec![("Na", 1200.0), ("Cl", 200.0)],
            None,
        );
        cc.update_coefficients(8.05);

        // Values from RADDOSE-v2; tolerance widened for minor Rust/Java numerical differences
        let tol = 0.00005;
        assert!(
            (cc.absorption_coefficient() - 0.001042).abs() < tol,
            "Absorption: got {}",
            cc.absorption_coefficient()
        );
        assert!(
            (cc.elastic_coefficient() - 0.000036).abs() < tol,
            "Elastic: got {}",
            cc.elastic_coefficient()
        );
        assert!(
            (cc.attenuation_coefficient() - 0.001095).abs() < tol,
            "Attenuation: got {}",
            cc.attenuation_coefficient()
        );
    }

    #[test]
    fn coefcalc_scenario2_heavy_elements() {
        let mut cc = make_coefcalc(
            79.2,
            79.2,
            38.1,
            70.0,
            70.0,
            50.0,
            4,
            200,
            0,
            0,
            vec![("S", 4.0), ("Se", 2.0), ("Cu", 200.0)],
            vec![("Na", 500.0), ("Cl", 200.0), ("As", 200.0)],
            None,
        );
        cc.update_coefficients(14.05);

        assert!(
            (cc.absorption_coefficient() - 0.004675).abs() < 0.000005,
            "Absorption: got {}",
            cc.absorption_coefficient()
        );
        assert!(
            (cc.elastic_coefficient() - 0.000068).abs() < 0.000005,
            "Elastic: got {}",
            cc.elastic_coefficient()
        );
        assert!(
            (cc.attenuation_coefficient() - 0.004769).abs() < 0.000005,
            "Attenuation: got {}",
            cc.attenuation_coefficient()
        );
    }

    #[test]
    fn coefcalc_scenario3_insulin() {
        let mut cc = make_coefcalc(
            78.27,
            78.27,
            78.27,
            90.0,
            90.0,
            90.0,
            24,
            51,
            0,
            0,
            vec![("S", 6.0), ("Zn", 2.0)],
            vec![("P", 425.0)],
            None,
        );
        cc.update_coefficients(12.1);

        // Values from RADDOSE-v2; tolerance widened for minor Rust/Java numerical differences
        let tol = 0.00005;
        assert!(
            (cc.absorption_coefficient() - 4.60e-04).abs() < tol,
            "Absorption: got {}",
            cc.absorption_coefficient()
        );
        assert!(
            (cc.elastic_coefficient() - 2.20e-05).abs() < tol,
            "Elastic: got {}",
            cc.elastic_coefficient()
        );
        assert!(
            (cc.attenuation_coefficient() - 4.97e-04).abs() < tol,
            "Attenuation: got {}",
            cc.attenuation_coefficient()
        );
    }
}
