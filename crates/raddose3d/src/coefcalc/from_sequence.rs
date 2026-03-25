use crate::parser::config::CrystalConfig;
use crate::residue::{get_residue_by_one_letter, TYPE_PROTEIN, TYPE_RNA, TYPE_DNA};

use super::compute::{
    CoefCalcCompute,
    CARBONS_PER_CARBOHYDRATE, OXYGENS_PER_CARBOHYDRATE, HYDROGENS_PER_CARBOHYDRATE,
};
use super::CoefCalc;

/// CoefCalcFromSequence: compute coefficients from a FASTA-style sequence file.
///
/// Sequence file format:
/// - Comment lines start with `>`. If the comment contains "ISDNA" or "ISRNA",
///   subsequent residues are treated as DNA/RNA; otherwise as protein.
/// - Sequence lines contain one-letter residue codes.
#[derive(Debug)]
pub struct CoefCalcFromSequence {
    pub compute: CoefCalcCompute,
    pub total_molecular_weight: f64,
}

impl CoefCalcFromSequence {
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let cell_a = config.cell_a.ok_or("Sequence requires cell a")?;
        let cell_b = config.cell_b.ok_or("Sequence requires cell b")?;
        let cell_c = config.cell_c.ok_or("Sequence requires cell c")?;
        let cell_alpha = config.cell_alpha.unwrap_or(90.0);
        let cell_beta = config.cell_beta.unwrap_or(90.0);
        let cell_gamma = config.cell_gamma.unwrap_or(90.0);
        let seq_file = config.seq_file.as_deref().ok_or("Sequence requires seq_file")?;
        let num_monomers = config.num_monomers.unwrap_or(1).max(1);
        let num_carb = config.num_carb.unwrap_or(0);

        let mut compute = CoefCalcCompute::new();
        compute.calculate_cell_volume(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);
        compute.num_monomers = num_monomers;

        // Heavy protein atoms (per monomer)
        for ec in &config.heavy_protein_atoms {
            compute.increment_macro(&ec.symbol.to_uppercase(), ec.count);
        }

        // Solvent concentrations
        let solv_names: Vec<String> = config.heavy_solution_conc.iter().map(|e| e.symbol.clone()).collect();
        let solv_concs: Vec<f64> = config.heavy_solution_conc.iter().map(|e| e.count).collect();
        if !solv_names.is_empty() {
            compute.add_solvent_concentrations(&solv_names, &solv_concs);
        }

        let mut total_mw = 0.0;
        // Parse sequence file
        eprintln!("Parsing sequence file: {seq_file}");
        parse_sequence_file(seq_file, &mut compute, &mut total_mw)?;
        eprintln!("Number of Amino Acids: {:.0}", compute.num_amino_acids);
        if compute.num_dna > 0.0 {
            eprintln!("Number of DNA Residues: {:.0}", compute.num_dna);
        }
        if compute.num_rna > 0.0 {
            eprintln!("Number of RNA Residues: {:.0}", compute.num_rna);
        }

        let mut sf = config.solvent_fraction.unwrap_or(-1.0);
        if sf <= 0.0 {
            sf = compute.calculate_solvent_fraction_from_nums();
        } else {
            compute.sol_fraction = sf;
        }
        compute.calculate_solvent_water(sf);

        // Add carbohydrates
        let m = num_monomers as f64;
        compute.increment_macro("C", CARBONS_PER_CARBOHYDRATE    * num_carb as f64);
        compute.increment_macro("O", OXYGENS_PER_CARBOHYDRATE    * num_carb as f64);
        compute.increment_macro("H", HYDROGENS_PER_CARBOHYDRATE  * num_carb as f64);
        compute.num_carb = num_carb as f64;

        // Multiply everything by num_monomers
        compute.multiply_atoms(m);

        compute.calculate_density();
        Ok(CoefCalcFromSequence { compute, total_molecular_weight: total_mw })
    }
}

/// Parse a FASTA-style sequence file, accumulating atom counts into `compute`.
pub fn parse_sequence_file(
    path: &str,
    compute: &mut CoefCalcCompute,
    total_mw: &mut f64,
) -> Result<(), String> {
    let content =
        std::fs::read_to_string(path).map_err(|e| format!("Cannot read sequence file: {e}"))?;
    let mut residue_type = TYPE_PROTEIN;
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            let upper = line.to_uppercase();
            if upper.contains("ISDNA") {
                residue_type = TYPE_DNA;
            } else if upper.contains("ISRNA") {
                residue_type = TYPE_RNA;
            } else {
                residue_type = TYPE_PROTEIN;
            }
        } else {
            for ch in line.chars() {
                let s = ch.to_string();
                if let Some(res) = get_residue_by_one_letter(&s, residue_type) {
                    *total_mw += res.molecular_weight;
                    match residue_type {
                        TYPE_PROTEIN => compute.num_amino_acids += 1.0,
                        TYPE_RNA     => compute.num_rna += 1.0,
                        TYPE_DNA     => compute.num_dna += 1.0,
                        _ => {}
                    }
                    compute.increment_macro("H",  res.hydrogens);
                    compute.increment_macro("O",  res.oxygens);
                    compute.increment_macro("C",  res.carbons);
                    compute.increment_macro("N",  res.nitrogens);
                    compute.increment_macro("P",  res.phosphoruses);
                    compute.increment_macro("S",  res.sulphurs);
                    if res.seleniums > 0.0 {
                        compute.increment_macro("SE", res.seleniums);
                    }
                } else {
                    eprintln!("Warning: unrecognised residue '{ch}'");
                }
            }
        }
    }
    Ok(())
}

impl CoefCalc for CoefCalcFromSequence {
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
