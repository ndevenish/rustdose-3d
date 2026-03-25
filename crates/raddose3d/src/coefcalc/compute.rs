use std::collections::{HashMap, HashSet};

use crate::element::{ElementDatabase, CrossSection};

/// Constants for amino acid/nucleotide/carbohydrate composition.
pub const PROTEIN_DENSITY: f64 = 1.35;
pub const RNA_DENSITY: f64 = 1.3;
pub const DNA_DENSITY: f64 = 1.35;
pub const CARBOHYDRATE_DENSITY: f64 = 1.54;
pub const HETATM_DENSITY: f64 = 1.35;

pub const AMINO_ACID_AVE_MASS: f64 = 110.0;
pub const DNA_NUCLEOTIDE_MASS: f64 = 312.0;
pub const RNA_NUCLEOTIDE_MASS: f64 = 321.0;
pub const CARBOHYDRATE_AVE_MASS: f64 = 162.0;

pub const HYDROGENS_PER_AMINO_ACID: f64 = 8.0;
pub const CARBONS_PER_AMINO_ACID: f64 = 5.0;
pub const NITROGENS_PER_AMINO_ACID: f64 = 1.35;
pub const OXYGENS_PER_AMINO_ACID: f64 = 1.5;

pub const HYDROGENS_PER_RNA_NUCLEOTIDE: f64 = 11.25;
pub const CARBONS_PER_RNA_NUCLEOTIDE: f64 = 9.5;
pub const NITROGENS_PER_RNA_NUCLEOTIDE: f64 = 3.75;
pub const OXYGENS_PER_RNA_NUCLEOTIDE: f64 = 7.0;
pub const PHOSPHORI_PER_RNA_NUCLEOTIDE: f64 = 1.0;

pub const HYDROGENS_PER_DNA_NUCLEOTIDE: f64 = 11.75;
pub const CARBONS_PER_DNA_NUCLEOTIDE: f64 = 9.75;
pub const NITROGENS_PER_DNA_NUCLEOTIDE: f64 = 4.0;
pub const OXYGENS_PER_DNA_NUCLEOTIDE: f64 = 6.0;
pub const PHOSPHORI_PER_DNA_NUCLEOTIDE: f64 = 1.0;

pub const HYDROGENS_PER_CARBOHYDRATE: f64 = 11.0;
pub const CARBONS_PER_CARBOHYDRATE: f64 = 6.0;
pub const OXYGENS_PER_CARBOHYDRATE: f64 = 5.0;

pub const ATOMIC_MASS_UNIT: f64 = 1.66e-24;
pub const AVOGADRO_NUM: f64 = 6.022e23;
pub const ANGSTROMS_TO_ML: f64 = 1e-24;
pub const MASS_TO_CELL_VOLUME: f64 = 1e27;
pub const WATER_CONCENTRATION: f64 = 51666.0;
pub const UNITS_PER_MILLI_UNIT: f64 = 1000.0;
pub const UNITS_PER_DECI_UNIT: f64 = 10.0;
pub const NUM_FLUOR_ESCAPE_FACTORS: usize = 27;

/// CoefCalcCompute: core absorption coefficient calculator.
#[derive(Debug)]
pub struct CoefCalcCompute {
    // Coefficients
    pub abs_coeff_photo: f64,
    pub abs_coeff_comp: f64,
    pub att_coeff: f64,
    pub elas_coeff: f64,
    pub elas_coeff_macro: f64,
    pub crystal_density: f64,
    pub molecular_weight: f64,

    // Cryo coefficients
    pub cryo_abs_coeff_photo: f64,
    pub cryo_abs_coeff_comp: f64,
    pub cryo_att_coeff: f64,
    pub cryo_elas_coeff: f64,
    pub cryo_density: f64,

    // Cell volume in Å³
    pub cell_volume: f64,

    // Solvent fraction
    pub sol_fraction: f64,

    // Composition
    pub num_amino_acids: f64,
    pub num_rna: f64,
    pub num_dna: f64,
    pub num_carb: f64,
    pub num_monomers: i32,

    // Element occurrences: element name → count
    pub macromolecular_occurrence: HashMap<String, f64>,
    pub solvent_occurrence: HashMap<String, f64>,
    pub solvent_concentration: HashMap<String, f64>,
    pub hetero_atom_occurrence: HashMap<String, f64>,
    pub cryo_occurrence: HashMap<String, f64>,
    pub cryo_concentration: HashMap<String, f64>,

    // Present elements
    pub present_elements: HashSet<String>,
    pub cryo_elements: HashSet<String>,
}

impl CoefCalcCompute {
    pub fn new() -> Self {
        CoefCalcCompute {
            abs_coeff_photo: 0.0,
            abs_coeff_comp: 0.0,
            att_coeff: 0.0,
            elas_coeff: 0.0,
            elas_coeff_macro: 0.0,
            crystal_density: 0.0,
            molecular_weight: 0.0,
            cryo_abs_coeff_photo: 0.0,
            cryo_abs_coeff_comp: 0.0,
            cryo_att_coeff: 0.0,
            cryo_elas_coeff: 0.0,
            cryo_density: 0.0,
            cell_volume: 0.0,
            sol_fraction: 0.0,
            num_amino_acids: 0.0,
            num_rna: 0.0,
            num_dna: 0.0,
            num_carb: 0.0,
            num_monomers: 1,
            macromolecular_occurrence: HashMap::new(),
            solvent_occurrence: HashMap::new(),
            solvent_concentration: HashMap::new(),
            hetero_atom_occurrence: HashMap::new(),
            cryo_occurrence: HashMap::new(),
            cryo_concentration: HashMap::new(),
            present_elements: HashSet::new(),
            cryo_elements: HashSet::new(),
        }
    }

    /// Total atoms of an element (macromolecular + solvent).
    pub fn total_atoms(&self, element_name: &str) -> f64 {
        self.macromolecular_occurrence
            .get(element_name)
            .copied()
            .unwrap_or(0.0)
            + self
                .solvent_occurrence
                .get(element_name)
                .copied()
                .unwrap_or(0.0)
    }

    /// Increment macromolecular occurrence.
    pub fn increment_macro(&mut self, element_name: &str, count: f64) {
        *self
            .macromolecular_occurrence
            .entry(element_name.to_string())
            .or_insert(0.0) += count;
    }

    /// Increment solvent occurrence.
    pub fn increment_solvent(&mut self, element_name: &str, count: f64) {
        *self
            .solvent_occurrence
            .entry(element_name.to_string())
            .or_insert(0.0) += count;
    }

    /// Set solvent occurrence.
    pub fn set_solvent_occurrence(&mut self, element_name: &str, count: f64) {
        self.solvent_occurrence
            .insert(element_name.to_string(), count);
    }

    /// Get solvent occurrence.
    pub fn get_solvent_occurrence(&self, element_name: &str) -> f64 {
        self.solvent_occurrence
            .get(element_name)
            .copied()
            .unwrap_or(0.0)
    }

    /// Calculate unit cell volume and return the value.
    pub fn calculate_cell_volume_ret(
        &mut self,
        a: f64, b: f64, c: f64,
        alpha_deg: f64, beta_deg: f64, gamma_deg: f64,
    ) -> f64 {
        self.calculate_cell_volume(a, b, c, alpha_deg, beta_deg, gamma_deg);
        self.cell_volume
    }

    /// Calculate unit cell volume from dimensions and angles.
    pub fn calculate_cell_volume(
        &mut self,
        a: f64, b: f64, c: f64,
        alpha_deg: f64, beta_deg: f64, gamma_deg: f64,
    ) {
        let alpha = alpha_deg.to_radians();
        let beta = beta_deg.to_radians();
        let gamma = gamma_deg.to_radians();

        let ult = 1.0
            + 2.0 * alpha.cos() * beta.cos() * gamma.cos()
            - alpha.cos().powi(2)
            - beta.cos().powi(2)
            - gamma.cos().powi(2);

        if ult < 0.0 {
            eprintln!("Warning: error calculating unit cell volume - please check inputs.");
        }

        self.cell_volume = a * b * c * ult.abs().sqrt();
        eprintln!("Cell volume: {:.2} Angstroms cubed.", self.cell_volume);
    }

    /// Calculate crystal density from composition.
    pub fn calculate_density(&mut self) {
        let db = ElementDatabase::instance();

        // Collect all present elements
        self.present_elements.clear();
        for name in self.macromolecular_occurrence.keys() {
            self.present_elements.insert(name.clone());
        }
        for name in self.solvent_occurrence.keys() {
            self.present_elements.insert(name.clone());
        }

        let mut mass = 0.0;
        self.molecular_weight = 0.0;

        for name in &self.present_elements {
            if let Some(e) = db.get(name) {
                let total = self.total_atoms(name);
                mass += total * e.atomic_weight_in_grams();
                self.molecular_weight += total * e.atomic_weight();
            }
        }

        self.crystal_density = mass * MASS_TO_CELL_VOLUME / (self.cell_volume * UNITS_PER_MILLI_UNIT);
    }

    /// Calculate solvent fraction from composition numbers.
    pub fn calculate_solvent_fraction_from_nums(&mut self) -> f64 {
        let protein_mass = ATOMIC_MASS_UNIT * AMINO_ACID_AVE_MASS * self.num_amino_acids
            * self.num_monomers as f64
            / (self.cell_volume * PROTEIN_DENSITY * ANGSTROMS_TO_ML);

        let rna_mass = ATOMIC_MASS_UNIT * RNA_NUCLEOTIDE_MASS * self.num_rna
            * self.num_monomers as f64
            / (self.cell_volume * RNA_DENSITY * ANGSTROMS_TO_ML);

        let dna_mass = ATOMIC_MASS_UNIT * DNA_NUCLEOTIDE_MASS * self.num_dna
            * self.num_monomers as f64
            / (self.cell_volume * DNA_DENSITY * ANGSTROMS_TO_ML);

        let carb_mass = ATOMIC_MASS_UNIT * CARBOHYDRATE_AVE_MASS * self.num_carb
            * self.num_monomers as f64
            / (self.cell_volume * CARBOHYDRATE_DENSITY * ANGSTROMS_TO_ML);

        let sf = 1.0 - protein_mass - rna_mass - dna_mass - carb_mass;

        if sf < 0.0 {
            eprintln!("Warning: Solvent mass calculated as a negative number...");
        }

        eprintln!("Solvent fraction determined as {:.2}%.", sf * 100.0);
        self.sol_fraction = sf;
        sf
    }

    /// Convert solvent concentrations to atom counts and add water.
    pub fn calculate_solvent_water(&mut self, solvent_fraction: f64) {
        let _db = ElementDatabase::instance();
        let mut non_water_atoms = 0.0;

        // Convert concentrations to atom counts
        let conc_snapshot: Vec<(String, f64)> = self
            .solvent_concentration
            .iter()
            .map(|(k, v)| (k.clone(), *v))
            .collect();

        for (name, conc) in &conc_snapshot {
            let atom_count =
                conc * AVOGADRO_NUM * self.cell_volume * solvent_fraction * 1e-3 * 1e-27;
            self.increment_solvent(name, atom_count);
            non_water_atoms += atom_count;
        }

        // Calculate water molecules
        let water_molecules = (WATER_CONCENTRATION * AVOGADRO_NUM / UNITS_PER_MILLI_UNIT
            * self.cell_volume
            * (1.0 / MASS_TO_CELL_VOLUME)
            * solvent_fraction
            - non_water_atoms)
            .max(0.0);

        // Add water: 2H + 1O per molecule
        let h_current = self.get_solvent_occurrence("H");
        self.set_solvent_occurrence("H", h_current + water_molecules * 2.0);

        let o_current = self.get_solvent_occurrence("O");
        self.set_solvent_occurrence("O", o_current + water_molecules);
    }

    /// Add solvent heavy atom concentrations.
    pub fn add_solvent_concentrations(&mut self, names: &[String], concs: &[f64]) {
        for (name, &conc) in names.iter().zip(concs.iter()) {
            *self
                .solvent_concentration
                .entry(name.clone())
                .or_insert(0.0) += conc;
        }
    }

    /// Calculate coefficients for all elements at a given energy.
    /// Returns (photoelectric, coherent, compton, total) in µm⁻¹.
    pub fn calculate_coefficients_all(&self, energy: f64) -> (f64, f64, f64, f64) {
        let db = ElementDatabase::instance();
        let mut photo = 0.0;
        let mut coherent = 0.0;
        let mut total = 0.0;
        let mut compton = 0.0;

        for name in &self.present_elements {
            if let Some(e) = db.get(name) {
                let cs = e.get_abs_coefficients(energy);
                let atoms = self.total_atoms(name);

                photo += atoms * cs[&CrossSection::Photoelectric] / self.cell_volume
                    / UNITS_PER_DECI_UNIT;
                coherent += atoms * cs[&CrossSection::Coherent] / self.cell_volume
                    / UNITS_PER_DECI_UNIT;
                total += atoms * cs[&CrossSection::Total] / self.cell_volume
                    / UNITS_PER_DECI_UNIT;
                compton += atoms * cs[&CrossSection::Compton] / self.cell_volume
                    / UNITS_PER_DECI_UNIT;
            }
        }

        (
            photo / UNITS_PER_MILLI_UNIT,
            coherent / UNITS_PER_MILLI_UNIT,
            compton / UNITS_PER_MILLI_UNIT,
            total / UNITS_PER_MILLI_UNIT,
        )
    }

    /// Calculate coefficients for macromolecular atoms only.
    pub fn calculate_coefficients_macro(&self, energy: f64) -> (f64, f64, f64, f64) {
        let db = ElementDatabase::instance();
        let mut photo = 0.0;
        let mut coherent = 0.0;
        let mut total = 0.0;
        let mut compton = 0.0;

        for name in &self.present_elements {
            if let Some(e) = db.get(name) {
                let cs = e.get_abs_coefficients(energy);
                let atoms = self
                    .macromolecular_occurrence
                    .get(name.as_str())
                    .copied()
                    .unwrap_or(0.0);

                photo +=
                    atoms * cs[&CrossSection::Photoelectric] / self.cell_volume / UNITS_PER_DECI_UNIT;
                coherent +=
                    atoms * cs[&CrossSection::Coherent] / self.cell_volume / UNITS_PER_DECI_UNIT;
                total +=
                    atoms * cs[&CrossSection::Total] / self.cell_volume / UNITS_PER_DECI_UNIT;
                compton +=
                    atoms * cs[&CrossSection::Compton] / self.cell_volume / UNITS_PER_DECI_UNIT;
            }
        }

        (
            photo / UNITS_PER_MILLI_UNIT,
            coherent / UNITS_PER_MILLI_UNIT,
            compton / UNITS_PER_MILLI_UNIT,
            total / UNITS_PER_MILLI_UNIT,
        )
    }

    /// Multiply all macromolecular atom occurrences by a factor (used by PDB/Sequence).
    pub fn multiply_atoms(&mut self, factor: f64) {
        for val in self.macromolecular_occurrence.values_mut() {
            *val *= factor;
        }
    }

    /// Set macromolecular occurrence (absolute, not increment).
    pub fn set_macro(&mut self, element_name: &str, count: f64) {
        self.macromolecular_occurrence
            .insert(element_name.to_string(), count);
    }

    /// Get macromolecular occurrence.
    pub fn get_macro(&self, element_name: &str) -> f64 {
        self.macromolecular_occurrence
            .get(element_name)
            .copied()
            .unwrap_or(0.0)
    }

    /// Add cryo-solution concentrations (stub — populates cryo_occurrence for density calc).
    pub fn add_cryo_concentrations(
        &mut self,
        cryo_names: &[String],
        cryo_concs: &[f64],
    ) {
        for (name, &conc) in cryo_names.iter().zip(cryo_concs.iter()) {
            *self
                .cryo_concentration
                .entry(name.clone())
                .or_insert(0.0) += conc;
        }
        // Populate cryo_occurrence for density calculation
        let vol = self.cell_volume;
        let sf = self.sol_fraction.max(0.0);
        for (name, &conc) in self.cryo_concentration.clone().iter() {
            let count = conc * AVOGADRO_NUM * vol * sf * 1e-3 * 1e-27;
            *self.cryo_occurrence.entry(name.clone()).or_insert(0.0) += count;
            self.cryo_elements.insert(name.clone());
        }
    }

    /// Calculate fluorescent escape factors for present elements.
    pub fn calc_fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>> {
        let db = ElementDatabase::instance();
        let mut factors = Vec::new();

        for name in &self.present_elements {
            if let Some(e) = db.get(name) {
                let mut row = vec![0.0; NUM_FLUOR_ESCAPE_FACTORS];

                // mu_ratio: element absorption / total absorption
                let el_cs = e.get_abs_coefficients(beam_energy);
                let el_photo = self.total_atoms(name)
                    * el_cs[&CrossSection::Photoelectric]
                    / self.cell_volume
                    / UNITS_PER_DECI_UNIT
                    / UNITS_PER_MILLI_UNIT;

                if self.abs_coeff_photo > 0.0 {
                    row[0] = el_photo / self.abs_coeff_photo;
                }

                let k_edge = e.k_edge().unwrap_or(0.0);
                let l1_edge = e.l1_edge().unwrap_or(0.0);
                let l2_edge = e.l2_edge().unwrap_or(0.0);
                let l3_edge = e.l3_edge().unwrap_or(0.0);
                let m1_edge_val = e.m1_edge();
                let m2_edge = e.m2_edge().unwrap_or(0.0);
                let m3_edge = e.m3_edge().unwrap_or(0.0);
                let m4_edge = e.m4_edge().unwrap_or(0.0);
                let m5_edge = e.m5_edge().unwrap_or(0.0);

                // K shell
                if beam_energy > k_edge && k_edge > 0.0 {
                    row[1] = e.k_ionisation_prob();
                    row[2] = e.k_fluorescence_yield().unwrap_or(0.0);
                    row[3] = k_edge;
                    if let Some(k_fl_avg) = e.k_fl_average() {
                        let (k_photo, _, _, _) = self.calculate_coefficients_all(k_fl_avg);
                        row[4] = k_photo;
                    }
                }

                // L1 shell
                if beam_energy > l1_edge && l1_edge > 0.0 && e.atomic_number() >= 12 {
                    row[5] = e.l1_ionisation_prob() * (1.0 - row[1]);
                    row[6] = e.l1_fluorescence_yield().unwrap_or(0.0);
                    row[7] = l1_edge;
                    if let Some(l_fl_avg) = e.l_fl_average() {
                        let (l_photo, _, _, _) = self.calculate_coefficients_all(l_fl_avg);
                        row[8] = l_photo;
                    }
                }

                // L2 shell
                if beam_energy > l2_edge && l2_edge > 0.0 && e.atomic_number() >= 12 {
                    row[9] = e.l2_ionisation_prob() * (1.0 - row[1] - row[5]);
                    row[10] = e.l2_fluorescence_yield().unwrap_or(0.0);
                    row[11] = l2_edge;
                    if let Some(l_fl_avg) = e.l_fl_average() {
                        let (l_photo, _, _, _) = self.calculate_coefficients_all(l_fl_avg);
                        row[12] = l_photo;
                    }
                }

                // L3 shell
                if beam_energy > l3_edge && l3_edge > 0.0 && e.atomic_number() >= 12 {
                    row[13] = e.l3_ionisation_prob() * (1.0 - row[1] - row[5] - row[9]);
                    row[14] = e.l3_fluorescence_yield().unwrap_or(0.0);
                    row[15] = l3_edge;
                    if let Some(l_fl_avg) = e.l_fl_average() {
                        let (l_photo, _, _, _) = self.calculate_coefficients_all(l_fl_avg);
                        row[16] = l_photo;
                    }
                }

                // M shells (Z >= 73 for heavy elements)
                if beam_energy > m1_edge_val && m1_edge_val > 0.0 && e.atomic_number() >= 73 {
                    row[17] = e.m1_ionisation_prob()
                        * (1.0 - row[1] - row[5] - row[9] - row[13]);
                    row[18] = m1_edge_val;
                }
                if beam_energy > m2_edge && m2_edge > 0.0 && e.atomic_number() >= 73 {
                    row[19] = e.m2_ionisation_prob()
                        * (1.0 - row[1] - row[5] - row[9] - row[13] - row[17]);
                    row[20] = m2_edge;
                }
                if beam_energy > m3_edge && m3_edge > 0.0 && e.atomic_number() >= 73 {
                    row[21] = e.m3_ionisation_prob()
                        * (1.0 - row[1] - row[5] - row[9] - row[13] - row[17] - row[19]);
                    row[22] = m3_edge;
                }
                if beam_energy > m4_edge && m4_edge > 0.0 && e.atomic_number() >= 73 {
                    row[23] = e.m4_ionisation_prob()
                        * (1.0 - row[1] - row[5] - row[9] - row[13] - row[17] - row[19] - row[21]);
                    row[24] = m4_edge;
                }
                if beam_energy > m5_edge && m5_edge > 0.0 && e.atomic_number() >= 73 {
                    row[25] = e.m5_ionisation_prob()
                        * (1.0 - row[1] - row[5] - row[9] - row[13] - row[17] - row[19] - row[21] - row[23]);
                    row[26] = m5_edge;
                }

                factors.push(row);
            }
        }

        factors
    }
}
