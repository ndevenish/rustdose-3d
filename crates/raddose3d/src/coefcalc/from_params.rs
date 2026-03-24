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
        let cell_a = config.cell_a.ok_or("Crystal requires unit cell dimension a")?;
        let cell_b = config.cell_b.ok_or("Crystal requires unit cell dimension b")?;
        let cell_c = config.cell_c.ok_or("Crystal requires unit cell dimension c")?;
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
