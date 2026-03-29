use crate::parser::config::CrystalConfig;
use crate::residue::get_residue_by_three_letter;

use super::compute::CoefCalcCompute;
use super::CoefCalc;

// PDB column positions (0-based, Rust slicing = [start..end])
const CRYST1_A: (usize, usize) = (6, 15);
const CRYST1_B: (usize, usize) = (15, 24);
const CRYST1_C: (usize, usize) = (24, 33);
const CRYST1_ALPHA: (usize, usize) = (33, 40);
const CRYST1_BETA: (usize, usize) = (40, 47);
const CRYST1_GAMMA: (usize, usize) = (47, 54);

const ATOM_RESIDUE_NAME: (usize, usize) = (17, 20);
const OCCUPANCY: (usize, usize) = (54, 60);
const ELEMENT_SYMBOL: (usize, usize) = (76, 78);

const SEQRES_START: usize = 19;
const SEQRES_RESI_LEN: usize = 3;

const SMTRY1_POS: (usize, usize) = (13, 19);
const MATRX_FLAG: (usize, usize) = (59, 60);

/// CoefCalcFromPDB: compute coefficients from a PDB file or 4-letter RCSB code.
///
/// Downloads from `https://files.rcsb.org/view/<CODE>.pdb` if the input is not
/// a `.pdb` file path.
#[derive(Debug)]
pub struct CoefCalcFromPDB {
    pub compute: CoefCalcCompute,
}

impl CoefCalcFromPDB {
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let pdb_code = config
            .pdb
            .as_deref()
            .ok_or("PDB mode requires pdb code or path")?;

        // Solvent concentrations (added before parsing)
        let mut compute = CoefCalcCompute::new();
        let solv_names: Vec<String> = config
            .heavy_solution_conc
            .iter()
            .map(|e| e.symbol.clone())
            .collect();
        let solv_concs: Vec<f64> = config.heavy_solution_conc.iter().map(|e| e.count).collect();
        if !solv_names.is_empty() {
            compute.add_solvent_concentrations(&solv_names, &solv_concs);
        }

        // Load PDB content
        let content = if pdb_code.ends_with(".pdb") {
            std::fs::read_to_string(pdb_code)
                .map_err(|e| format!("Cannot read PDB file '{pdb_code}': {e}"))?
        } else {
            download_pdb(&pdb_code.to_uppercase())?
        };

        let mut cs_sym = 0usize;
        let mut ncs_sym = 1usize;
        let mut occupancy_warned = false;
        let mut found_cryst1 = false;

        for line in content.lines() {
            if line.len() < 6 {
                continue;
            }
            let directive = &line[..6];

            match directive {
                "CRYST1" => {
                    if parse_cryst1(line, &mut compute) {
                        found_cryst1 = true;
                    }
                }
                "HETATM" => {
                    parse_hetatm(line, &mut compute, &mut occupancy_warned);
                }
                "SEQRES" => {
                    parse_seqres(line, &mut compute);
                }
                "REMARK" => {
                    if line.len() >= SMTRY1_POS.1 && &line[SMTRY1_POS.0..SMTRY1_POS.1] == "SMTRY1" {
                        cs_sym += 1;
                    }
                }
                "MTRIX1" => {
                    if line.len() > MATRX_FLAG.1 {
                        if &line[MATRX_FLAG.0..MATRX_FLAG.1] == "1" {
                            // Coordinates already present — don't count this NCS op
                        } else {
                            ncs_sym += 1;
                        }
                    }
                }
                _ => {}
            }
        }

        if !found_cryst1 {
            println!("Warning: no CRYST1 line found in PDB.");
        }

        println!("Crystallographic symmetry operators: {cs_sym}");
        println!("Non-crystallographic symmetry operators: {ncs_sym}");

        let num_monomers = cs_sym * ncs_sym;
        compute.num_monomers = num_monomers as i32;
        compute.multiply_atoms(num_monomers as f64);
        println!("Number of monomers: {num_monomers}");

        let sf = compute.calculate_solvent_fraction_from_nums();
        compute.calculate_solvent_water(sf);
        compute.calculate_density();

        Ok(CoefCalcFromPDB { compute })
    }
}

fn parse_cryst1(line: &str, compute: &mut CoefCalcCompute) -> bool {
    let get = |r: (usize, usize)| line.get(r.0..r.1).map(str::trim);
    let a: f64 = get(CRYST1_A).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let b: f64 = get(CRYST1_B).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let c: f64 = get(CRYST1_C).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let alpha: f64 = get(CRYST1_ALPHA)
        .and_then(|s| s.parse().ok())
        .unwrap_or(90.0);
    let beta: f64 = get(CRYST1_BETA)
        .and_then(|s| s.parse().ok())
        .unwrap_or(90.0);
    let gamma: f64 = get(CRYST1_GAMMA)
        .and_then(|s| s.parse().ok())
        .unwrap_or(90.0);
    if a == 0.0 {
        return false;
    }
    println!("PDB unit cell: {a} {b} {c} {alpha} {beta} {gamma}");
    compute.calculate_cell_volume(a, b, c, alpha, beta, gamma);
    true
}

fn parse_hetatm(line: &str, compute: &mut CoefCalcCompute, warned: &mut bool) {
    if line.len() < ELEMENT_SYMBOL.1 {
        return;
    }
    let residue_name = line
        .get(ATOM_RESIDUE_NAME.0..ATOM_RESIDUE_NAME.1)
        .map(str::trim)
        .unwrap_or("")
        .to_uppercase();
    if residue_name == "HOH" {
        return;
    }

    let occ_str = line
        .get(OCCUPANCY.0..OCCUPANCY.1)
        .map(str::trim)
        .unwrap_or("");
    let element_sym = line
        .get(ELEMENT_SYMBOL.0..ELEMENT_SYMBOL.1)
        .map(str::trim)
        .unwrap_or("")
        .to_uppercase();

    let occ: f64 = if occ_str.is_empty() {
        if !*warned {
            println!("Warning: occupancy missing, assuming 1.0");
            *warned = true;
        }
        1.0
    } else {
        occ_str.parse().unwrap_or(1.0)
    };

    if element_sym.is_empty() {
        return;
    }
    compute.increment_macro(&element_sym, occ);
    // Also track in hetero_atom_occurrence
    *compute
        .hetero_atom_occurrence
        .entry(element_sym)
        .or_insert(0.0) += occ;
}

fn parse_seqres(line: &str, compute: &mut CoefCalcCompute) {
    if line.len() <= SEQRES_START {
        return;
    }
    let mut seq = &line[SEQRES_START..];
    while seq.len() >= SEQRES_RESI_LEN {
        let three = seq[..SEQRES_RESI_LEN].trim();
        seq = &seq[SEQRES_RESI_LEN..];
        // skip leading space after residue
        if seq.starts_with(' ') {
            seq = &seq[1..];
        }
        if three.is_empty() {
            continue;
        }
        if let Some(res) = get_residue_by_three_letter(three) {
            match res.residue_type {
                crate::residue::TYPE_PROTEIN => compute.num_amino_acids += 1.0,
                crate::residue::TYPE_RNA => compute.num_rna += 1.0,
                crate::residue::TYPE_DNA => compute.num_dna += 1.0,
                _ => {}
            }
            compute.increment_macro("H", res.hydrogens);
            compute.increment_macro("O", res.oxygens);
            compute.increment_macro("C", res.carbons);
            compute.increment_macro("N", res.nitrogens);
            compute.increment_macro("P", res.phosphoruses);
            compute.increment_macro("S", res.sulphurs);
            if res.seleniums > 0.0 {
                compute.increment_macro("SE", res.seleniums);
            }
        } else {
            println!("Warning: could not decipher PDB 3-letter code: '{three}'");
        }
    }
}

fn download_pdb(code: &str) -> Result<String, String> {
    let url = format!("https://files.rcsb.org/view/{code}.pdb");
    println!("Downloading PDB from {url}");
    ureq::get(&url)
        .call()
        .map_err(|e| format!("Failed to download PDB '{code}': {e}"))?
        .body_mut()
        .read_to_string()
        .map_err(|e| format!("Failed to read PDB response: {e}"))
}

impl CoefCalc for CoefCalcFromPDB {
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
