use crate::parser::config::CrystalConfig;

use super::compute::{
    CoefCalcCompute, CARBONS_PER_AMINO_ACID, CARBONS_PER_CARBOHYDRATE, CARBONS_PER_DNA_NUCLEOTIDE,
    CARBONS_PER_RNA_NUCLEOTIDE, HYDROGENS_PER_AMINO_ACID, HYDROGENS_PER_CARBOHYDRATE,
    HYDROGENS_PER_DNA_NUCLEOTIDE, HYDROGENS_PER_RNA_NUCLEOTIDE, NITROGENS_PER_AMINO_ACID,
    NITROGENS_PER_DNA_NUCLEOTIDE, NITROGENS_PER_RNA_NUCLEOTIDE, OXYGENS_PER_AMINO_ACID,
    OXYGENS_PER_CARBOHYDRATE, OXYGENS_PER_DNA_NUCLEOTIDE, OXYGENS_PER_RNA_NUCLEOTIDE,
    PHOSPHORI_PER_DNA_NUCLEOTIDE, PHOSPHORI_PER_RNA_NUCLEOTIDE,
};
use super::CoefCalc;

/// CoefCalcMicroED: electron diffraction (MicroED) coefficient calculator.
///
/// Same composition logic as CoefCalcFromParams but sets `is_em = true`
/// to signal that electron rather than X-ray cross-sections should be used.
#[derive(Debug)]
pub struct CoefCalcMicroED {
    pub compute: CoefCalcCompute,
    pub is_em: bool,
}

impl CoefCalcMicroED {
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let mut compute = CoefCalcCompute::new();

        let cell_a = config.cell_a.ok_or("MicroED requires cell a")?;
        let cell_b = config.cell_b.ok_or("MicroED requires cell b")?;
        let cell_c = config.cell_c.ok_or("MicroED requires cell c")?;
        let cell_alpha = config.cell_alpha.unwrap_or(90.0);
        let cell_beta = config.cell_beta.unwrap_or(90.0);
        let cell_gamma = config.cell_gamma.unwrap_or(90.0);

        compute.calculate_cell_volume(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);

        let num_monomers = config.num_monomers.unwrap_or(1).max(1);
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
            compute.increment_macro(&ec.symbol.to_uppercase(), ec.count * num_monomers as f64);
        }

        // Solvent concentrations
        let solv_names: Vec<String> = config
            .heavy_solution_conc
            .iter()
            .map(|e| e.symbol.clone())
            .collect();
        let solv_concs: Vec<f64> = config.heavy_solution_conc.iter().map(|e| e.count).collect();
        if !solv_names.is_empty() {
            compute.add_solvent_concentrations(&solv_names, &solv_concs);
        }

        let mut sf = config.solvent_fraction.unwrap_or(-1.0);
        if sf <= 0.0 {
            sf = compute.calculate_solvent_fraction_from_nums();
        } else {
            compute.sol_fraction = sf;
        }
        compute.calculate_solvent_water(sf);

        // Standard composition (same as CoefCalcFromParams)
        let m = num_monomers as f64;
        compute.increment_macro("C", CARBONS_PER_AMINO_ACID * num_residues as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_AMINO_ACID * num_residues as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_AMINO_ACID * num_residues as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_AMINO_ACID * num_residues as f64 * m);

        compute.increment_macro("C", CARBONS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("P", PHOSPHORI_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);

        compute.increment_macro("C", CARBONS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("P", PHOSPHORI_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);

        compute.increment_macro("C", CARBONS_PER_CARBOHYDRATE * num_carb as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_CARBOHYDRATE * num_carb as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_CARBOHYDRATE * num_carb as f64 * m);

        compute.calculate_density();
        Ok(CoefCalcMicroED {
            compute,
            is_em: true,
        })
    }
}

impl CoefCalc for CoefCalcMicroED {
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
}
