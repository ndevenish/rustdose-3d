use super::compute::{
    CoefCalcCompute, AVOGADRO_NUM, CARBONS_PER_CARBOHYDRATE, HYDROGENS_PER_CARBOHYDRATE,
    OXYGENS_PER_CARBOHYDRATE,
};
use super::from_sequence::parse_sequence_file;
use super::CoefCalc;
use crate::parser::config::CrystalConfig;

const ANGSTROM_TO_LITRE: f64 = 1e-27;
const DEFAULT_CELL_LENGTH: f64 = 1000.0;

/// CoefCalcFromSequenceSAXS: SAXS mode using a sequence file for composition.
///
/// Like CoefCalcFromSequence but the number of monomers is calculated from
/// the protein concentration (g/L) and the actual molecular weight from the
/// sequence, rather than using average residue masses.
#[derive(Debug)]
pub struct CoefCalcFromSequenceSAXS {
    pub compute: CoefCalcCompute,
    pub total_molecular_weight: f64,
}

impl CoefCalcFromSequenceSAXS {
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let cell_a = config.cell_a.unwrap_or(DEFAULT_CELL_LENGTH);
        let cell_b = config.cell_b.unwrap_or(DEFAULT_CELL_LENGTH);
        let cell_c = config.cell_c.unwrap_or(DEFAULT_CELL_LENGTH);
        let cell_alpha = config.cell_alpha.unwrap_or(90.0);
        let cell_beta = config.cell_beta.unwrap_or(90.0);
        let cell_gamma = config.cell_gamma.unwrap_or(90.0);
        let seq_file = config
            .seq_file
            .as_deref()
            .ok_or("SequenceSAXS requires seq_file")?;
        let protein_conc = config
            .protein_conc
            .ok_or("SequenceSAXS requires protein_conc")?;
        let num_carb = config.num_carb.unwrap_or(0);

        let mut compute = CoefCalcCompute::new();
        let vol = compute
            .calculate_cell_volume_ret(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);

        // Parse sequence to get composition and molecular weight
        println!("Parsing sequence file: {seq_file}");
        let mut total_mw = 0.0;
        parse_sequence_file(seq_file, &mut compute, &mut total_mw)?;
        println!("Number of Amino Acids: {:.0}", compute.num_amino_acids);
        if compute.num_dna > 0.0 {
            println!("Number of DNA Residues: {:.0}", compute.num_dna);
        }
        if compute.num_rna > 0.0 {
            println!("Number of RNA Residues: {:.0}", compute.num_rna);
        }
        println!("Total molecular weight: {:.2} g/mol", total_mw);

        // Calculate num_monomers from protein concentration and actual MW
        let num_monomers = calculate_num_monomers_from_mw(protein_conc, total_mw, vol);
        compute.num_monomers = num_monomers;

        // Heavy protein atoms
        for ec in &config.heavy_protein_atoms {
            compute.increment_macro(&ec.symbol.to_uppercase(), ec.count);
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

        // Add carbohydrates (pre-multiply)
        compute.increment_macro("C", CARBONS_PER_CARBOHYDRATE * num_carb as f64);
        compute.increment_macro("O", OXYGENS_PER_CARBOHYDRATE * num_carb as f64);
        compute.increment_macro("H", HYDROGENS_PER_CARBOHYDRATE * num_carb as f64);
        compute.num_carb = num_carb as f64;

        // Multiply by num_monomers
        compute.multiply_atoms(num_monomers as f64);

        compute.calculate_density();
        Ok(CoefCalcFromSequenceSAXS {
            compute,
            total_molecular_weight: total_mw,
        })
    }
}

fn calculate_num_monomers_from_mw(
    protein_conc: f64,
    total_mw: f64,
    cell_volume_angstrom3: f64,
) -> i32 {
    let molarity = protein_conc / total_mw;
    let volume_litres = ANGSTROM_TO_LITRE * cell_volume_angstrom3;
    let num = (molarity * volume_litres * AVOGADRO_NUM).round();
    if num < 1.0 {
        println!("WARNING: calculated monomers < 1; increase unit cell size. Using 1.");
        return 1;
    }
    println!(
        "Calculated number of monomers in cell volume: {}",
        num as i32
    );
    num as i32
}

impl CoefCalc for CoefCalcFromSequenceSAXS {
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
