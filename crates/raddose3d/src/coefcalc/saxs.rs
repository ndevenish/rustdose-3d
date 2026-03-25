use crate::parser::config::CrystalConfig;

use super::compute::{
    CoefCalcCompute,
    CARBONS_PER_AMINO_ACID, NITROGENS_PER_AMINO_ACID, OXYGENS_PER_AMINO_ACID,
    HYDROGENS_PER_AMINO_ACID,
    CARBONS_PER_RNA_NUCLEOTIDE, NITROGENS_PER_RNA_NUCLEOTIDE, OXYGENS_PER_RNA_NUCLEOTIDE,
    HYDROGENS_PER_RNA_NUCLEOTIDE, PHOSPHORI_PER_RNA_NUCLEOTIDE,
    CARBONS_PER_DNA_NUCLEOTIDE, NITROGENS_PER_DNA_NUCLEOTIDE, OXYGENS_PER_DNA_NUCLEOTIDE,
    HYDROGENS_PER_DNA_NUCLEOTIDE, PHOSPHORI_PER_DNA_NUCLEOTIDE,
    CARBONS_PER_CARBOHYDRATE, OXYGENS_PER_CARBOHYDRATE, HYDROGENS_PER_CARBOHYDRATE,
    AVOGADRO_NUM,
};
use super::CoefCalc;

/// Average molecular masses used to estimate num_monomers from protein concentration.
const AVG_RESIDUE_MASS: f64 = 110.0;
const AVG_RNA_MASS: f64 = 339.5;
const AVG_DNA_MASS: f64 = 327.0;
const AVG_CARB_MASS: f64 = 180.0;
/// Conversion factor: Å³ → litres.
const ANGSTROM_TO_LITRE: f64 = 1e-27;
/// Default unit cell dimension when not specified (SAXS uses large volume).
const DEFAULT_CELL_LENGTH: f64 = 1000.0;

/// CoefCalcSAXS: calculates coefficients for a SAXS experiment.
///
/// The number of monomers is derived from the protein concentration (g/L)
/// rather than being given explicitly.
#[derive(Debug)]
pub struct CoefCalcSAXS {
    pub compute: CoefCalcCompute,
}

impl CoefCalcSAXS {
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let mut compute = CoefCalcCompute::new();

        let cell_a = config.cell_a.unwrap_or(DEFAULT_CELL_LENGTH);
        let cell_b = config.cell_b.unwrap_or(DEFAULT_CELL_LENGTH);
        let cell_c = config.cell_c.unwrap_or(DEFAULT_CELL_LENGTH);
        let cell_alpha = config.cell_alpha.unwrap_or(90.0);
        let cell_beta  = config.cell_beta.unwrap_or(90.0);
        let cell_gamma = config.cell_gamma.unwrap_or(90.0);

        let vol = compute.calculate_cell_volume_ret(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);

        let num_residues = config.num_residues.unwrap_or(0);
        let num_rna      = config.num_rna.unwrap_or(0);
        let num_dna      = config.num_dna.unwrap_or(0);
        let num_carb     = config.num_carb.unwrap_or(0);
        let protein_conc = config.protein_conc.ok_or("SAXS requires protein_conc")?;

        compute.num_amino_acids = num_residues as f64;
        compute.num_rna = num_rna as f64;
        compute.num_dna = num_dna as f64;
        compute.num_carb = num_carb as f64;

        let num_monomers = calculate_num_monomers_saxs(
            num_residues, num_rna, num_dna, num_carb, protein_conc, vol,
        );
        compute.num_monomers = num_monomers;

        // Heavy protein atoms
        for ec in &config.heavy_protein_atoms {
            compute.increment_macro(&ec.symbol.to_uppercase(), ec.count * num_monomers as f64);
        }

        // Solvent concentrations
        let solv_names: Vec<String> = config.heavy_solution_conc.iter().map(|e| e.symbol.clone()).collect();
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

        let m = num_monomers as f64;
        compute.increment_macro("C", CARBONS_PER_AMINO_ACID    * num_residues as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_AMINO_ACID  * num_residues as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_AMINO_ACID    * num_residues as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_AMINO_ACID  * num_residues as f64 * m);

        compute.increment_macro("C", CARBONS_PER_RNA_NUCLEOTIDE   * num_rna as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_RNA_NUCLEOTIDE   * num_rna as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);
        compute.increment_macro("P", PHOSPHORI_PER_RNA_NUCLEOTIDE * num_rna as f64 * m);

        compute.increment_macro("C", CARBONS_PER_DNA_NUCLEOTIDE   * num_dna as f64 * m);
        compute.increment_macro("N", NITROGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_DNA_NUCLEOTIDE   * num_dna as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);
        compute.increment_macro("P", PHOSPHORI_PER_DNA_NUCLEOTIDE * num_dna as f64 * m);

        compute.increment_macro("C", CARBONS_PER_CARBOHYDRATE   * num_carb as f64 * m);
        compute.increment_macro("O", OXYGENS_PER_CARBOHYDRATE   * num_carb as f64 * m);
        compute.increment_macro("H", HYDROGENS_PER_CARBOHYDRATE * num_carb as f64 * m);

        compute.calculate_density();
        Ok(CoefCalcSAXS { compute })
    }
}

/// Calculate num_monomers from protein concentration (g/L) and cell volume (Å³).
pub fn calculate_num_monomers_saxs(
    num_residues: i32, num_rna: i32, num_dna: i32, num_carb: i32,
    protein_conc: f64, cell_volume_angstrom3: f64,
) -> i32 {
    let total_mass = AVG_RESIDUE_MASS * num_residues as f64
        + AVG_RNA_MASS * num_rna as f64
        + AVG_DNA_MASS * num_dna as f64
        + AVG_CARB_MASS * num_carb as f64;
    let molarity = protein_conc / total_mass;
    let volume_litres = ANGSTROM_TO_LITRE * cell_volume_angstrom3;
    let num = (molarity * volume_litres * AVOGADRO_NUM).round();
    if num < 1.0 {
        eprintln!("WARNING: calculated monomers < 1; increase unit cell size. Using 1.");
        return 1;
    }
    eprintln!("Calculated number of monomers in cell volume: {}", num as i32);
    num as i32
}

impl CoefCalc for CoefCalcSAXS {
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

    fn absorption_coefficient(&self) -> f64 { self.compute.abs_coeff_photo }
    fn attenuation_coefficient(&self) -> f64 { self.compute.att_coeff }
    fn elastic_coefficient(&self) -> f64 { self.compute.elas_coeff }
    fn inelastic_coefficient(&self) -> f64 { self.compute.abs_coeff_comp }
    fn density(&self) -> f64 { self.compute.crystal_density }
    fn fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>> {
        self.compute.calc_fluorescent_escape_factors(beam_energy)
    }
    fn solvent_fraction(&self) -> f64 { self.compute.sol_fraction }
}
