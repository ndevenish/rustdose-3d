use std::collections::{HashMap, HashSet};

use crate::constants::KEV_TO_JOULES;
use crate::element::{CrossSection, ElementDatabase, ElementDatabaseEM};

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

    /// Number of simulated photons/electrons for MC/XFEL (0 = use default 1,000,000).
    pub num_simulated_electrons: u64,

    /// Molecular weight of the cryo-surrounding (sum of cryo occurrence × atomic weight).
    pub molecular_weight_surrounding: f64,

    /// Per-element elastic cross-section × number density (nm⁻¹), last computed for crystal.
    pub elastic_x_sections: HashMap<String, f64>,
    /// Sum of all per-element elastic cross-section × number density (nm⁻¹), crystal.
    pub elastic_x_section_tot: f64,
    /// Same as above but for surrounding/cryo.
    pub elastic_x_sections_surrounding: HashMap<String, f64>,
    pub elastic_x_section_tot_surrounding: f64,

    // GOS inelastic state: element → [maxShells=10][4]  (long, trans, close, tot)
    pub gos_inelastic: HashMap<String, Vec<[f64; 4]>>,
    pub gos_inelastic_surrounding: HashMap<String, Vec<[f64; 4]>>,
    /// Plasmon (conduction-band) inelastic: [long, trans, close, tot]
    pub cb_inel: [f64; 4],
    pub cb_inel_surrounding: [f64; 4],
    /// Sturnheimer a parameter from last `gos_inel` call.
    pub sturnheimer_adjustment: f64,
}

impl Default for CoefCalcCompute {
    fn default() -> Self {
        Self::new()
    }
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
            num_simulated_electrons: 0,
            molecular_weight_surrounding: 0.0,
            elastic_x_sections: HashMap::new(),
            elastic_x_section_tot: 0.0,
            elastic_x_sections_surrounding: HashMap::new(),
            elastic_x_section_tot_surrounding: 0.0,
            gos_inelastic: HashMap::new(),
            gos_inelastic_surrounding: HashMap::new(),
            cb_inel: [0.0; 4],
            cb_inel_surrounding: [0.0; 4],
            sturnheimer_adjustment: 1.5,
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
        a: f64,
        b: f64,
        c: f64,
        alpha_deg: f64,
        beta_deg: f64,
        gamma_deg: f64,
    ) -> f64 {
        self.calculate_cell_volume(a, b, c, alpha_deg, beta_deg, gamma_deg);
        self.cell_volume
    }

    /// Calculate unit cell volume from dimensions and angles.
    pub fn calculate_cell_volume(
        &mut self,
        a: f64,
        b: f64,
        c: f64,
        alpha_deg: f64,
        beta_deg: f64,
        gamma_deg: f64,
    ) {
        let alpha = alpha_deg.to_radians();
        let beta = beta_deg.to_radians();
        let gamma = gamma_deg.to_radians();

        let ult = 1.0 + 2.0 * alpha.cos() * beta.cos() * gamma.cos()
            - alpha.cos().powi(2)
            - beta.cos().powi(2)
            - gamma.cos().powi(2);

        if ult < 0.0 {
            println!("Warning: error calculating unit cell volume - please check inputs.");
        }

        self.cell_volume = a * b * c * ult.abs().sqrt();
        println!("Cell volume: {:.2} Angstroms cubed.", self.cell_volume);
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

        self.crystal_density =
            mass * MASS_TO_CELL_VOLUME / (self.cell_volume * UNITS_PER_MILLI_UNIT);
    }

    /// Calculate solvent fraction from composition numbers.
    pub fn calculate_solvent_fraction_from_nums(&mut self) -> f64 {
        let protein_mass = ATOMIC_MASS_UNIT
            * AMINO_ACID_AVE_MASS
            * self.num_amino_acids
            * self.num_monomers as f64
            / (self.cell_volume * PROTEIN_DENSITY * ANGSTROMS_TO_ML);

        let rna_mass =
            ATOMIC_MASS_UNIT * RNA_NUCLEOTIDE_MASS * self.num_rna * self.num_monomers as f64
                / (self.cell_volume * RNA_DENSITY * ANGSTROMS_TO_ML);

        let dna_mass =
            ATOMIC_MASS_UNIT * DNA_NUCLEOTIDE_MASS * self.num_dna * self.num_monomers as f64
                / (self.cell_volume * DNA_DENSITY * ANGSTROMS_TO_ML);

        let carb_mass =
            ATOMIC_MASS_UNIT * CARBOHYDRATE_AVE_MASS * self.num_carb * self.num_monomers as f64
                / (self.cell_volume * CARBOHYDRATE_DENSITY * ANGSTROMS_TO_ML);

        // Hetatm mass: light atoms (atomic number < 29) only, matching Java behaviour.
        let db = ElementDatabase::instance();
        let mut hetatm_total_mass = 0.0;
        for (sym, count) in &self.hetero_atom_occurrence {
            if let Some(elem) = db.get(sym) {
                if elem.atomic_number() < 29 {
                    hetatm_total_mass += count * elem.atomic_weight_in_grams();
                }
            }
        }
        let hetatm_mass = hetatm_total_mass / (self.cell_volume * HETATM_DENSITY * ANGSTROMS_TO_ML);

        let sf = 1.0 - protein_mass - rna_mass - dna_mass - carb_mass - hetatm_mass;

        if sf < 0.0 {
            println!("Warning: Solvent mass calculated as a negative number...");
        }

        println!("Solvent fraction determined as {:.2}%.", sf * 100.0);
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

                photo += atoms * cs[&CrossSection::Photoelectric]
                    / self.cell_volume
                    / UNITS_PER_DECI_UNIT;
                coherent +=
                    atoms * cs[&CrossSection::Coherent] / self.cell_volume / UNITS_PER_DECI_UNIT;
                total += atoms * cs[&CrossSection::Total] / self.cell_volume / UNITS_PER_DECI_UNIT;
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

    /// Per-element cumulative photoelectric absorption probabilities.
    /// Returns HashMap<element_name, cumulative_prob> where probs are sorted by element name.
    /// Mirrors Java's CoefCalcCompute.getPhotoElectricProbsElement().
    pub fn calc_photo_electric_probs_element(
        &self,
        energy: f64,
    ) -> std::collections::HashMap<String, f64> {
        let db = ElementDatabase::instance();
        let total_photo = self.abs_coeff_photo;
        let mut result = std::collections::HashMap::new();
        if total_photo <= 0.0 {
            return result;
        }
        let mut elements: Vec<&String> = self.present_elements.iter().collect();
        elements.sort();
        let mut running = 0.0;
        for name in elements {
            if let Some(e) = db.get(name) {
                let cs = e.get_abs_coefficients(energy);
                let atoms = self.total_atoms(name);
                let photo = atoms * cs[&CrossSection::Photoelectric]
                    / self.cell_volume
                    / UNITS_PER_DECI_UNIT
                    / UNITS_PER_MILLI_UNIT;
                running += photo / total_photo;
                result.insert(name.clone(), running);
            }
        }
        result
    }

    /// Per-element cumulative Compton scattering probabilities.
    pub fn calc_compton_probs_element(
        &self,
        energy: f64,
    ) -> std::collections::HashMap<String, f64> {
        let db = ElementDatabase::instance();
        let total_compton = self.abs_coeff_comp;
        let mut result = std::collections::HashMap::new();
        if total_compton <= 0.0 {
            return result;
        }
        let mut elements: Vec<&String> = self.present_elements.iter().collect();
        elements.sort();
        let mut running = 0.0;
        for name in elements {
            if let Some(e) = db.get(name) {
                let cs = e.get_abs_coefficients(energy);
                let atoms = self.total_atoms(name);
                let compton = atoms * cs[&CrossSection::Compton]
                    / self.cell_volume
                    / UNITS_PER_DECI_UNIT
                    / UNITS_PER_MILLI_UNIT;
                running += compton / total_compton;
                result.insert(name.clone(), running);
            }
        }
        result
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

                photo += atoms * cs[&CrossSection::Photoelectric]
                    / self.cell_volume
                    / UNITS_PER_DECI_UNIT;
                coherent +=
                    atoms * cs[&CrossSection::Coherent] / self.cell_volume / UNITS_PER_DECI_UNIT;
                total += atoms * cs[&CrossSection::Total] / self.cell_volume / UNITS_PER_DECI_UNIT;
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
        for val in self.hetero_atom_occurrence.values_mut() {
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

    /// Add cryo-solution concentrations.
    /// Matches Java `CoefCalcCompute.addCryoConcentrations()`.
    ///
    /// Parameters:
    /// - `cryo_solution_names` / `cryo_solution_concs`: heavy atom concentrations in surrounding
    /// - `oil_based`: if "true", use density-based path (no water fill)
    /// - `oil_element_names` / `oil_element_nums`: molecular formula of surrounding (density-based)
    /// - `oil_density`: density of surrounding material in g/mL (density-based)
    pub fn add_cryo_concentrations(
        &mut self,
        cryo_solution_names: &[String],
        cryo_solution_concs: &[f64],
        oil_based: Option<&str>,
        oil_element_names: &[String],
        oil_element_nums: &[f64],
        oil_density: f64,
    ) {
        let mut non_water_atoms = 0.0;

        let calc_water = !oil_based
            .map(|s| s.eq_ignore_ascii_case("true"))
            .unwrap_or(false);

        // Add cryo solution concentrations (heavy atom concentrations in surrounding)
        for (name, &conc) in cryo_solution_names.iter().zip(cryo_solution_concs.iter()) {
            *self.cryo_concentration.entry(name.clone()).or_insert(0.0) += conc;
        }

        // If oil-based (density-based), compute element concentrations from molecular formula
        if !calc_water && !oil_element_names.is_empty() {
            let db = ElementDatabase::instance();
            let mut g_per_mol = 0.0;
            for (name, &num) in oil_element_names.iter().zip(oil_element_nums.iter()) {
                if let Some(e) = db.get(name) {
                    g_per_mol += e.atomic_weight() * num;
                }
            }
            // concentration in µmol/L → multiply by 1E6 to get per-litre
            let oil_conc = (oil_density / g_per_mol) * 1e6;
            for (name, &num) in oil_element_names.iter().zip(oil_element_nums.iter()) {
                let element_conc = oil_conc * num;
                *self.cryo_concentration.entry(name.clone()).or_insert(0.0) += element_conc;
            }
        }

        // Convert concentrations to atom counts per unit cell
        let vol = self.cell_volume;
        for (name, &conc) in self.cryo_concentration.clone().iter() {
            let atom_count = conc * AVOGADRO_NUM * vol * 1e-3 * 1e-27;
            *self.cryo_occurrence.entry(name.clone()).or_insert(0.0) += atom_count;
            non_water_atoms += atom_count;
        }

        // Add water if not oil-based
        if calc_water {
            let water_molecules = WATER_CONCENTRATION * AVOGADRO_NUM / UNITS_PER_MILLI_UNIT
                * vol
                * (1.0 / MASS_TO_CELL_VOLUME)
                - non_water_atoms;

            *self.cryo_occurrence.entry("H".to_string()).or_insert(0.0) += water_molecules * 2.0;
            *self.cryo_occurrence.entry("O".to_string()).or_insert(0.0) += water_molecules;
        }

        // Populate cryo_elements set
        self.cryo_elements = self.cryo_occurrence.keys().cloned().collect();
    }

    /// Calculate cryo-solution density from its composition.
    /// Matches Java `CoefCalcCompute.calculateCryoDensity()`.
    pub fn calculate_cryo_density(&mut self) {
        let db = ElementDatabase::instance();
        let mut mass = 0.0;
        self.molecular_weight_surrounding = 0.0;

        for name in &self.cryo_elements.clone() {
            if let Some(e) = db.get(name) {
                let occ = self.cryo_occurrence.get(name).copied().unwrap_or(0.0);
                mass += occ * e.atomic_weight_in_grams();
                self.molecular_weight_surrounding += occ * e.atomic_weight();
            }
        }

        self.cryo_density = mass * MASS_TO_CELL_VOLUME / (self.cell_volume * UNITS_PER_MILLI_UNIT);
    }

    /// Whether cryo elements have been populated.
    pub fn is_cryo(&self) -> bool {
        !self.cryo_elements.is_empty()
    }

    /// Calculate cryo-solution absorption/attenuation coefficients.
    /// Matches Java `CoefCalcCompute.calculateCryoCoefficientsAll()`.
    pub fn calculate_cryo_coefficients_all(&self, energy: f64) -> (f64, f64, f64, f64) {
        let db = ElementDatabase::instance();
        let mut photo = 0.0;
        let mut coherent = 0.0;
        let mut total = 0.0;
        let mut compton = 0.0;

        for name in &self.cryo_elements {
            if let Some(e) = db.get(name) {
                let cs = e.get_abs_coefficients(energy);
                let occ = self.cryo_occurrence.get(name).copied().unwrap_or(0.0);

                photo +=
                    occ * cs[&CrossSection::Photoelectric] / self.cell_volume / UNITS_PER_DECI_UNIT;
                coherent +=
                    occ * cs[&CrossSection::Coherent] / self.cell_volume / UNITS_PER_DECI_UNIT;
                total += occ * cs[&CrossSection::Total] / self.cell_volume / UNITS_PER_DECI_UNIT;
                compton +=
                    occ * cs[&CrossSection::Compton] / self.cell_volume / UNITS_PER_DECI_UNIT;
            }
        }

        (
            photo / UNITS_PER_MILLI_UNIT,
            coherent / UNITS_PER_MILLI_UNIT,
            compton / UNITS_PER_MILLI_UNIT,
            total / UNITS_PER_MILLI_UNIT,
        )
    }

    /// Calculate coefficients for a single cryo element.
    /// Matches Java `CoefCalcCompute.calculateCoefficientsCryoElement()`.
    fn calculate_coefficients_cryo_element(
        &self,
        energy: f64,
        element_name: &str,
    ) -> (f64, f64, f64, f64) {
        let db = ElementDatabase::instance();
        let mut photo = 0.0;
        let mut coherent = 0.0;
        let mut total = 0.0;
        let mut compton = 0.0;

        if let Some(e) = db.get(element_name) {
            let cs = e.get_abs_coefficients(energy);
            let occ = self
                .cryo_occurrence
                .get(element_name)
                .copied()
                .unwrap_or(0.0);

            photo +=
                occ * cs[&CrossSection::Photoelectric] / self.cell_volume / UNITS_PER_DECI_UNIT;
            coherent += occ * cs[&CrossSection::Coherent] / self.cell_volume / UNITS_PER_DECI_UNIT;
            total += occ * cs[&CrossSection::Total] / self.cell_volume / UNITS_PER_DECI_UNIT;
            compton += occ * cs[&CrossSection::Compton] / self.cell_volume / UNITS_PER_DECI_UNIT;
        }

        (
            photo / UNITS_PER_MILLI_UNIT,
            coherent / UNITS_PER_MILLI_UNIT,
            compton / UNITS_PER_MILLI_UNIT,
            total / UNITS_PER_MILLI_UNIT,
        )
    }

    /// Update cryo coefficients for a given photon energy.
    /// Matches Java `CoefCalcCompute.updateCryoCoefficients(double)`.
    pub fn update_cryo_coefficients(&mut self, photon_energy: f64) {
        let (photo, coherent, compton, total) = self.calculate_cryo_coefficients_all(photon_energy);
        self.cryo_abs_coeff_photo = photo;
        self.cryo_elas_coeff = coherent;
        self.cryo_abs_coeff_comp = compton;
        self.cryo_att_coeff = total;
    }

    /// Calculate fluorescent escape factors for cryo elements.
    /// Matches Java `CoefCalcCompute.getCryoFluorescentEscapeFactors()`.
    /// Uses the same Rust index layout as `calc_fluorescent_escape_factors()`.
    pub fn calc_cryo_fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>> {
        let db = ElementDatabase::instance();
        let mut factors = Vec::new();

        for name in &self.cryo_elements {
            if let Some(e) = db.get(name) {
                let mut row = vec![0.0; NUM_FLUOR_ESCAPE_FACTORS];

                // mu_ratio: element cryo absorption / total cryo absorption
                let (el_photo, _, _, _) =
                    self.calculate_coefficients_cryo_element(beam_energy, name);

                if self.cryo_abs_coeff_photo > 0.0 {
                    row[0] = el_photo / self.cryo_abs_coeff_photo;
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
                let mut k_factor_a = 0.0;
                if beam_energy > k_edge && k_edge > 0.0 {
                    k_factor_a = e.k_ionisation_prob();
                    let k_factor_b = e.k_fluorescence_yield().unwrap_or(0.0);
                    // Escape muabs: compute crystal coefficients at K fluorescence energy
                    let escape_muabs_k = if let Some(k_fl_avg) = e.k_fl_average() {
                        let (k_photo, _, _, _) = self.calculate_coefficients_all(k_fl_avg);
                        k_photo
                    } else {
                        0.0
                    };
                    row[1] = k_factor_a;
                    row[2] = k_factor_b;
                    row[3] = k_edge;
                    row[4] = escape_muabs_k;
                }

                // L1 shell
                let mut l1_factor_a = 0.0;
                if beam_energy > l1_edge && l1_edge > 0.0 && e.atomic_number() >= 12 {
                    l1_factor_a = e.l1_ionisation_prob() * (1.0 - k_factor_a);
                    row[5] = l1_factor_a;
                    row[6] = e.l1_fluorescence_yield().unwrap_or(0.0);
                    row[7] = l1_edge;
                    if let Some(l_fl_avg) = e.l_fl_average() {
                        let (l_photo, _, _, _) = self.calculate_coefficients_all(l_fl_avg);
                        row[8] = l_photo;
                    }
                }

                // L2 shell
                let mut l2_factor_a = 0.0;
                if beam_energy > l2_edge && l2_edge > 0.0 && e.atomic_number() >= 12 {
                    l2_factor_a = e.l2_ionisation_prob() * (1.0 - k_factor_a - l1_factor_a);
                    row[9] = l2_factor_a;
                    row[10] = e.l2_fluorescence_yield().unwrap_or(0.0);
                    row[11] = l2_edge;
                    if let Some(l_fl_avg) = e.l_fl_average() {
                        let (l_photo, _, _, _) = self.calculate_coefficients_all(l_fl_avg);
                        row[12] = l_photo;
                    }
                }

                // L3 shell
                let mut l3_factor_a = 0.0;
                if beam_energy > l3_edge && l3_edge > 0.0 && e.atomic_number() >= 12 {
                    l3_factor_a =
                        e.l3_ionisation_prob() * (1.0 - k_factor_a - l1_factor_a - l2_factor_a);
                    row[13] = l3_factor_a;
                    row[14] = e.l3_fluorescence_yield().unwrap_or(0.0);
                    row[15] = l3_edge;
                    if let Some(l_fl_avg) = e.l_fl_average() {
                        let (l_photo, _, _, _) = self.calculate_coefficients_all(l_fl_avg);
                        row[16] = l_photo;
                    }
                }

                // M shells (Z >= 73)
                let mut m1_factor_a = 0.0;
                if beam_energy > m1_edge_val && m1_edge_val > 0.0 && e.atomic_number() >= 73 {
                    m1_factor_a = e.m1_ionisation_prob()
                        * (1.0 - k_factor_a - l1_factor_a - l2_factor_a - l3_factor_a);
                    row[17] = m1_factor_a;
                    row[18] = m1_edge_val;
                }
                let mut m2_factor_a = 0.0;
                if beam_energy > m2_edge && m2_edge > 0.0 && e.atomic_number() >= 73 {
                    m2_factor_a = e.m2_ionisation_prob()
                        * (1.0
                            - k_factor_a
                            - l1_factor_a
                            - l2_factor_a
                            - l3_factor_a
                            - m1_factor_a);
                    row[19] = m2_factor_a;
                    row[20] = m2_edge;
                }
                let mut m3_factor_a = 0.0;
                if beam_energy > m3_edge && m3_edge > 0.0 && e.atomic_number() >= 73 {
                    m3_factor_a = e.m3_ionisation_prob()
                        * (1.0
                            - k_factor_a
                            - l1_factor_a
                            - l2_factor_a
                            - l3_factor_a
                            - m1_factor_a
                            - m2_factor_a);
                    row[21] = m3_factor_a;
                    row[22] = m3_edge;
                }
                let mut m4_factor_a = 0.0;
                if beam_energy > m4_edge && m4_edge > 0.0 && e.atomic_number() >= 73 {
                    m4_factor_a = e.m4_ionisation_prob()
                        * (1.0
                            - k_factor_a
                            - l1_factor_a
                            - l2_factor_a
                            - l3_factor_a
                            - m1_factor_a
                            - m2_factor_a
                            - m3_factor_a);
                    row[23] = m4_factor_a;
                    row[24] = m4_edge;
                }
                if beam_energy > m5_edge && m5_edge > 0.0 && e.atomic_number() >= 73 {
                    let m5_factor_a = e.m5_ionisation_prob()
                        * (1.0
                            - k_factor_a
                            - l1_factor_a
                            - l2_factor_a
                            - l3_factor_a
                            - m1_factor_a
                            - m2_factor_a
                            - m3_factor_a
                            - m4_factor_a);
                    row[25] = m5_factor_a;
                    row[26] = m5_edge;
                }

                factors.push(row);
            }
        }

        factors
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
                let el_photo = self.total_atoms(name) * el_cs[&CrossSection::Photoelectric]
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
                    row[17] = e.m1_ionisation_prob() * (1.0 - row[1] - row[5] - row[9] - row[13]);
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
                        * (1.0
                            - row[1]
                            - row[5]
                            - row[9]
                            - row[13]
                            - row[17]
                            - row[19]
                            - row[21]
                            - row[23]);
                    row[26] = m5_edge;
                }

                factors.push(row);
            }
        }

        factors
    }

    /// Per-element relative shell ionisation probabilities at the given beam energy.
    /// Returns HashMap<element_name, Vec<cumulative_prob>> with 9 entries (K, L1-L3, M1-M5).
    /// Mirrors Java MC.getRelativeShellProbs().
    pub fn calc_relative_shell_probs(
        &self,
        energy: f64,
        cryo: bool,
    ) -> std::collections::HashMap<String, Vec<f64>> {
        let db = ElementDatabase::instance();
        let elements: &std::collections::HashSet<String> = if cryo {
            &self.cryo_elements
        } else {
            &self.present_elements
        };
        let mut result = std::collections::HashMap::new();

        for name in elements {
            if let Some(e) = db.get(name) {
                let mut shell_probs = vec![0.0_f64; 9];
                let mut running = 0.0_f64;

                let k_prob = if let Some(k) = e.k_edge() {
                    if energy > k {
                        let p = e.k_ionisation_prob();
                        running += p;
                        p
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[0] = running;

                let l1_prob = if e.atomic_number() >= 12 {
                    if let Some(l1) = e.l1_edge() {
                        if energy > l1 {
                            let p = e.l1_ionisation_prob() * (1.0 - k_prob);
                            running += p;
                            p
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[1] = running;

                let l2_prob = if e.atomic_number() >= 12 {
                    if let Some(l2) = e.l2_edge() {
                        if energy > l2 {
                            let p = e.l2_ionisation_prob() * (1.0 - k_prob - l1_prob);
                            running += p;
                            p
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[2] = running;

                let l3_prob = if e.atomic_number() >= 12 {
                    if let Some(l3) = e.l3_edge() {
                        if energy > l3 {
                            let p = e.l3_ionisation_prob() * (1.0 - k_prob - l1_prob - l2_prob);
                            running += p;
                            p
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[3] = running;

                let m1_prob = if e.atomic_number() >= 73 {
                    if energy > e.m1_edge() {
                        let p =
                            e.m1_ionisation_prob() * (1.0 - k_prob - l1_prob - l2_prob - l3_prob);
                        running += p;
                        p
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[4] = running;

                let m2_prob = if e.atomic_number() >= 73 {
                    if let Some(m2) = e.m2_edge() {
                        if energy > m2 {
                            let p = e.m2_ionisation_prob()
                                * (1.0 - k_prob - l1_prob - l2_prob - l3_prob - m1_prob);
                            running += p;
                            p
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[5] = running;

                let m3_prob = if e.atomic_number() >= 73 {
                    if let Some(m3) = e.m3_edge() {
                        if energy > m3 {
                            let p = e.m3_ionisation_prob()
                                * (1.0 - k_prob - l1_prob - l2_prob - l3_prob - m1_prob - m2_prob);
                            running += p;
                            p
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[6] = running;

                let m4_prob = if e.atomic_number() >= 73 {
                    if let Some(m4) = e.m4_edge() {
                        if energy > m4 {
                            let p = e.m4_ionisation_prob()
                                * (1.0
                                    - k_prob
                                    - l1_prob
                                    - l2_prob
                                    - l3_prob
                                    - m1_prob
                                    - m2_prob
                                    - m3_prob);
                            running += p;
                            p
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                shell_probs[7] = running;

                if e.atomic_number() >= 73 {
                    if let Some(m5) = e.m5_edge() {
                        if energy > m5 {
                            let p = e.m5_ionisation_prob()
                                * (1.0
                                    - k_prob
                                    - l1_prob
                                    - l2_prob
                                    - l3_prob
                                    - m1_prob
                                    - m2_prob
                                    - m3_prob
                                    - m4_prob);
                            running += p;
                        }
                    }
                }
                shell_probs[8] = running;

                result.insert(name.clone(), shell_probs);
            }
        }
        result
    }

    /// Shell binding energy (keV) for a given element and shell index.
    /// Mirrors Java MC.getShellBindingEnergy().
    /// Shell indices: 0=K, 1=L1, 2=L2, 3=L3, 4=M1, 5=M2, 6=M3, 7=M4, 8=M5
    pub fn calc_shell_binding_energy(&self, element: &str, shell: usize) -> f64 {
        let db = ElementDatabase::instance();
        if let Some(e) = db.get(element) {
            match shell {
                0 => e.k_edge().unwrap_or(0.0),
                1 => e.l1_edge().unwrap_or(0.0),
                2 => e.l2_edge().unwrap_or(0.0),
                3 => e.l3_edge().unwrap_or(0.0),
                4 => e.m1_edge(),
                5 => e.m2_edge().unwrap_or(0.0),
                6 => e.m3_edge().unwrap_or(0.0),
                7 => e.m4_edge().unwrap_or(0.0),
                8 => e.m5_edge().unwrap_or(0.0),
                _ => 0.0,
            }
        } else {
            0.0
        }
    }

    // ── MC electron transport: stopping power, elastic MFPL, elastic probs ───

    /// Bethe stopping power (keV/nm) including Joy-Luo correction and radiative
    /// component.  Matches Java `CoefCalcCompute.calcStoppingPower()`.
    pub fn calc_stopping_power(&self, energy_kev: f64, surrounding: bool) -> f64 {
        const M_E: f64 = 9.109_383_56e-31; // electron mass (kg)
        const C: f64 = 299_792_458.0; // speed of light (m/s)
        const C2: f64 = C * C;

        let db = ElementDatabase::instance();
        let (elements, density, mol_weight) = if surrounding {
            (
                &self.cryo_elements,
                self.cryo_density,
                self.molecular_weight_surrounding,
            )
        } else {
            (
                &self.present_elements,
                self.crystal_density,
                self.molecular_weight,
            )
        };

        if mol_weight <= 0.0 || density <= 0.0 {
            return 0.0;
        }

        let v0 = energy_kev * KEV_TO_JOULES;
        let beta_sq = 1.0 - (M_E * C2 / (v0 + M_E * C2)).powi(2);
        if beta_sq <= 0.0 {
            return 0.0;
        }
        let gamma = 1.0 / (1.0 - beta_sq).sqrt();
        let ke = (gamma - 1.0) * M_E * C2;
        let ke_mev = (ke / KEV_TO_JOULES) / 1000.0;

        let mut sum_z: f64 = 0.0;
        let mut sum_a: f64 = 0.0;
        let mut mean_ln_i: f64 = 0.0;
        let mut mean_z: f64 = 0.0;

        for name in elements {
            let e = match db.get(name) {
                Some(e) => e,
                None => continue,
            };
            let a = e.atomic_weight();
            let z = e.atomic_number() as f64;
            let num_atoms = if surrounding {
                self.cryo_occurrence.get(name).copied().unwrap_or(0.0)
            } else {
                self.total_atoms(name)
            };

            let mol_weight_fraction = (num_atoms * a) / mol_weight;
            sum_z += z * num_atoms;
            sum_a += a * num_atoms;
            mean_z += mol_weight_fraction * z;

            // Mean excitation energy with condensed-phase correction
            let mut j = e.mean_ionisation_potential(); // eV
            let zi = e.atomic_number();
            if zi != 1 && zi != 6 && zi != 7 && zi != 8 && zi != 9 && zi != 17 {
                j *= 1.13; // gas → solid/liquid phase
            }
            // Joy-Luo low-energy correction
            let k = 0.7344 * z.powf(0.0367);
            let j_star = j / (1.0 + k * (j / (energy_kev * 1000.0)));

            mean_ln_i += mol_weight_fraction * (z / a) * j_star.ln();
        }

        if sum_a <= 0.0 || sum_z <= 0.0 {
            return 0.0;
        }

        let z_over_a = sum_z / sum_a;
        mean_ln_i /= z_over_a;
        let mean_i = mean_ln_i.exp(); // eV

        // Joy-Luo correction on mean I
        let mean_j = mean_i / (1.0 + 0.85 * mean_i / (energy_kev * 1000.0));
        let mean_j_joules = (mean_j / 1000.0) * KEV_TO_JOULES;

        let energy_joules = energy_kev * KEV_TO_JOULES;
        let delta = 0.1832;

        // Modified Bethe formula (Fbeta)
        let f_beta = ((M_E * C2 * energy_joules * beta_sq) / (2.0 * (1.0 - beta_sq))).ln()
            - (2.0 * (1.0 - beta_sq).sqrt() - 1.0 + beta_sq) * 2.0_f64.ln()
            + 1.0
            - beta_sq
            + (1.0 / 8.0) * (1.0 - (1.0 - beta_sq).sqrt());

        let mut stopping_power =
            (0.153536 / beta_sq) * z_over_a * (f_beta - 2.0 * mean_j_joules.ln() - delta);
        stopping_power = stopping_power * 1000.0 * density / 1e7;

        // Radiative stopping power
        let radiative = (ke_mev * mean_z / 800.0) * stopping_power;
        stopping_power += radiative;

        stopping_power // keV/nm
    }

    /// Compute per-element elastic cross-section contributions and MFPL.
    /// Matches Java `getElectronElasticCrossSection()` + `getElectronElasticMFPL()`.
    ///
    /// As a side effect, populates `elastic_x_sections`/`elastic_x_section_tot`
    /// (or their `_surrounding` variants) which are later read by `elastic_probs_calc()`.
    ///
    /// Returns the elastic mean free path length in nm.
    pub fn calc_electron_elastic_mfpl(&mut self, electron_energy: f64, surrounding: bool) -> f64 {
        let db = ElementDatabase::instance();
        let db_em = ElementDatabaseEM::instance();

        let (elements, density, mol_weight) = if surrounding {
            (
                self.cryo_elements.clone(),
                self.cryo_density,
                self.molecular_weight_surrounding,
            )
        } else {
            (
                self.present_elements.clone(),
                self.crystal_density,
                self.molecular_weight,
            )
        };

        if mol_weight <= 0.0 {
            return 0.0;
        }

        const M_E: f64 = 9.109_383_56e-31;
        const C2: f64 = 299_792_458.0 * 299_792_458.0;

        let v0 = electron_energy * KEV_TO_JOULES;
        let beta_sq = 1.0 - (M_E * C2 / (v0 + M_E * C2)).powi(2);

        let mut part_lambda = 0.0;
        let mut x_section_tot_per_element = 0.0;
        let mut x_sections: HashMap<String, f64> = HashMap::new();

        let cell_vol_nm3 = self.cell_volume / 1000.0; // Å³ → nm³

        for name in &elements {
            let e_xray = match db.get(name) {
                Some(e) => e,
                None => continue,
            };
            let z = e_xray.atomic_number() as f64;
            let a = e_xray.atomic_weight();

            // Langmore & Smith formula (nm²/atom)
            let mut elastic_xs = if beta_sq > 0.0 {
                (1.4e-6 * z.powf(1.5) / beta_sq) * (1.0 - (0.26 * z) / (137.0 * beta_sq.sqrt()))
            } else {
                0.0
            };

            // mol weight fraction and number density (atoms/nm³)
            let (mol_weight_fraction, n_per_vol) = if surrounding {
                let occ = self.cryo_occurrence.get(name).copied().unwrap_or(0.0);
                (
                    (occ * a) / mol_weight,
                    if cell_vol_nm3 > 0.0 {
                        occ / cell_vol_nm3
                    } else {
                        0.0
                    },
                )
            } else {
                let occ = self.total_atoms(name);
                (
                    (occ * a) / mol_weight,
                    if cell_vol_nm3 > 0.0 {
                        occ / cell_vol_nm3
                    } else {
                        0.0
                    },
                )
            };

            // Override with ELSEPA tabulated value when available (0.05–300 keV)
            if (0.05..=300.0).contains(&electron_energy) {
                if let Some(e_em) = db_em.get(name) {
                    let xs = e_em.elastic_coefficient(electron_energy);
                    if xs > 0.0 {
                        elastic_xs = xs;
                    }
                }
            }

            let xs_n = elastic_xs * n_per_vol; // nm⁻¹
            x_section_tot_per_element += xs_n;
            x_sections.insert(name.clone(), xs_n);

            part_lambda += (mol_weight_fraction * elastic_xs) / a;
        }

        // Store side-effect state
        if surrounding {
            self.elastic_x_sections_surrounding = x_sections;
            self.elastic_x_section_tot_surrounding = x_section_tot_per_element;
        } else {
            self.elastic_x_sections = x_sections;
            self.elastic_x_section_tot = x_section_tot_per_element;
        }

        // MFPL = 1 / (N_A × (density/1e21) × partLambda)
        let denom = AVOGADRO_NUM * (density / 1e21) * part_lambda;
        if denom > 0.0 {
            1.0 / denom
        } else {
            0.0
        }
    }

    /// Cumulative elastic scattering probabilities per element.
    /// Matches Java `getElasticProbs()`.  Must be called after
    /// `calc_electron_elastic_mfpl()` which populates the side-effect state.
    pub fn elastic_probs_calc(&self, surrounding: bool) -> HashMap<String, f64> {
        let (xs_map, xs_tot) = if surrounding {
            (
                &self.elastic_x_sections_surrounding,
                self.elastic_x_section_tot_surrounding,
            )
        } else {
            (&self.elastic_x_sections, self.elastic_x_section_tot)
        };

        let mut probs = HashMap::new();
        if xs_tot <= 0.0 {
            return probs;
        }

        let mut running_sum = 0.0;
        for (name, &xs) in xs_map {
            running_sum += xs / xs_tot;
            probs.insert(name.clone(), running_sum);
        }
        probs
    }

    // -----------------------------------------------------------------------
    // GOS (Generalised Oscillator Strength) inelastic scattering physics
    // Ported from Java CoefCalcCompute.java — Ritchie/Ashley GOS model.
    // -----------------------------------------------------------------------

    /// Returns total atoms of an element in the cryo-surrounding.
    fn total_atoms_surr(&self, element: &str) -> f64 {
        self.cryo_occurrence.get(element).copied().unwrap_or(0.0)
    }

    /// Number of valence electrons and inner shells for element Z.
    /// Returns `[valence, numInnerShells]`.  Matches Java `getNumValenceElectronsSubshells`.
    pub fn num_valence_electrons_subshells_z(z: i32) -> [i32; 2] {
        if z <= 2 {
            [z, 0]
        } else if z <= 10 {
            [z - 2, 1]
        } else if z == 11 {
            [z - 4, 2]
        } else if z == 12 {
            [z - 6, 3]
        } else if (13..=19).contains(&z) {
            [z - 10, 4]
        } else if (20..=22).contains(&z) {
            [z - 12, 5]
        } else if z == 23 {
            [z - 14, 6]
        } else if (24..=32).contains(&z) {
            [z - 18, 7]
        } else if (33..=34).contains(&z) {
            [z - 28, 9]
        } else {
            // Fallback for higher Z — use shell structure from getNumValenceElectrons
            if z <= 10 {
                [z - 2, 1]
            } else if z <= 28 {
                [z - 10, 2]
            } else if z <= 60 {
                [z - 28, 3]
            } else {
                [z - 60, 4]
            }
        }
    }

    /// Shell binding energy in eV using subshell indices.
    /// Index 0=K, 1=L1, 2=L2, 3=L3, 4=M1, 5=M2, 6=M3, 7=M4, 8=M5, 9=N1.
    /// Matches Java `getShellBindingSubshell`.
    pub fn shell_binding_subshell_kev(element: &str, shell_index: usize) -> f64 {
        use crate::element::database::ElementDatabase;
        let db = ElementDatabase::instance();
        let Some(e) = db.get(element) else { return 0.0 };
        let val = match shell_index {
            0 => e.k_edge(),
            1 => e.l1_binding(),
            2 => e.l2_edge(),
            3 => e.l3_edge(),
            4 => e.m1_binding(),
            5 => e.m2_edge(),
            6 => e.m3_edge(),
            7 => e.m4_edge(),
            8 => e.m5_edge(),
            9 => e.n1_binding(),
            _ => None,
        };
        val.unwrap_or(0.0)
    }

    /// Shell binding energy in eV using coarser (3-shell) index.
    /// Index 0=K, 1=L1, 2=M1.  Matches Java `getShellBinding`.
    pub fn shell_binding_kev(element: &str, shell_index: usize) -> f64 {
        use crate::element::database::ElementDatabase;
        let db = ElementDatabase::instance();
        let Some(e) = db.get(element) else { return 0.0 };
        match shell_index {
            0 => e.k_edge().unwrap_or(0.0),
            1 => e.l1_edge().unwrap_or(0.0),
            2 => e.m1_edge(),
            _ => 0.0,
        }
    }

    /// Plasma energy (eV) from electron density.  Matches Java `getPlasmaEnergyAll`.
    pub fn plasma_energy_all(&self, surrounding: bool) -> f64 {
        let hbar_sq: f64 = (6.62607004e-34 / (2.0 * std::f64::consts::PI)).powi(2);
        let m: f64 = 9.10938356e-31;
        let e_sq: f64 = (4.80320425e-10_f64).powi(2) / 1000.0; // esu²→kg cm³ s⁻²
        let mut nz: f64 = 0.0;
        let elements: Vec<String> = if surrounding {
            self.cryo_elements.iter().cloned().collect()
        } else {
            self.present_elements.iter().cloned().collect()
        };
        for name in &elements {
            let atom_num = if surrounding {
                self.total_atoms_surr(name)
            } else {
                self.total_atoms(name)
            };
            if let Some(z) = crate::element::database::ElementDatabase::instance()
                .get(name)
                .map(|e| e.atomic_number() as f64)
            {
                // electrons cm⁻³ — cell_volume in Å³, 1Å³ = 1e-24 cm³
                nz += (z * atom_num) / (self.cell_volume / 1e24);
            }
        }
        let mut plasma_energy = (4.0 * std::f64::consts::PI * hbar_sq * nz * e_sq / m).sqrt(); // J
        plasma_energy = (plasma_energy / KEV_TO_JOULES) * 1000.0; // eV
        plasma_energy
    }

    /// Wcb — bulk conduction-band plasmon energy (eV).  Matches Java `getWcbAll`.
    pub fn calc_wcb_all(&self, surrounding: bool) -> f64 {
        const FCB_CUTOFF: f64 = 0.0;
        let hbar_sq: f64 = (6.62607004e-34 / (2.0 * std::f64::consts::PI)).powi(2);
        let m: f64 = 9.10938356e-31;
        let e_sq: f64 = (4.80320425e-10_f64).powi(2) / 1000.0;
        let elements: Vec<String> = if surrounding {
            self.cryo_elements.iter().cloned().collect()
        } else {
            self.present_elements.iter().cloned().collect()
        };
        let mut nz: f64 = 0.0;
        let mut sum_z: f64 = 0.0;
        let mut sum_fcb: f64 = 0.0;
        for name in &elements {
            let atom_num = if surrounding {
                self.total_atoms_surr(name)
            } else {
                self.total_atoms(name)
            };
            let db = crate::element::database::ElementDatabase::instance();
            let Some(e) = db.get(name) else { continue };
            let z = e.atomic_number() as f64;
            nz += (z * atom_num) / (self.cell_volume / 1e24);
            sum_z += z * atom_num;
            let [_valence, num_inner] = Self::num_valence_electrons_subshells_z(e.atomic_number());
            for i in 0..=num_inner {
                let uk_kev = Self::shell_binding_kev(name, i as usize);
                let uk_ev = uk_kev * 1000.0;
                if uk_ev < FCB_CUTOFF {
                    // subshells[i] electrons from the simple shell structure
                    const SUBSHELLS: [i32; 10] = [2, 2, 2, 4, 2, 2, 4, 10, 2, 30];
                    let [valence, num_inner2] =
                        Self::num_valence_electrons_subshells_z(e.atomic_number());
                    let fk = if i == num_inner2 {
                        valence
                    } else {
                        SUBSHELLS[i as usize]
                    };
                    sum_fcb += fk as f64 * atom_num;
                }
            }
        }
        if sum_z == 0.0 {
            return 0.0;
        }
        let mut plasma_energy = (4.0 * std::f64::consts::PI * hbar_sq * nz * e_sq / m).sqrt(); // J
        plasma_energy = (plasma_energy / KEV_TO_JOULES) * 1000.0; // eV
        ((sum_fcb / sum_z).sqrt()) * plasma_energy
    }

    /// Wk for a molecule shell.  Matches Java `getWkMolecule`.
    pub fn wk_molecule_calc(
        &self,
        a: f64,
        element: &str,
        shell_index: usize,
        surrounding: bool,
    ) -> f64 {
        const SUBSHELLS: [i32; 10] = [2, 2, 2, 4, 2, 2, 4, 10, 2, 30];
        let db = crate::element::database::ElementDatabase::instance();
        let Some(e) = db.get(element) else { return 0.0 };
        let z = e.atomic_number();
        let elements: Vec<String> = if surrounding {
            self.cryo_elements.iter().cloned().collect()
        } else {
            self.present_elements.iter().cloned().collect()
        };
        let tot_num = if surrounding {
            self.total_atoms_surr(element)
        } else {
            self.total_atoms(element)
        };
        let mut sum_z: i64 = 0;
        for name in &elements {
            let atom_num = if surrounding {
                self.total_atoms_surr(name)
            } else {
                self.total_atoms(name)
            };
            if let Some(elem) = db.get(name) {
                sum_z += (elem.atomic_number() as i64) * (atom_num.round() as i64);
            }
        }
        if sum_z == 0 {
            return 0.0;
        }
        let [valence, num_inner] = Self::num_valence_electrons_subshells_z(z);
        let fk = if shell_index == num_inner as usize {
            valence
        } else {
            SUBSHELLS.get(shell_index).copied().unwrap_or(0)
        };
        let plasma_energy = self.plasma_energy_all(surrounding);
        let uk_kev = Self::shell_binding_subshell_kev(element, shell_index);
        let term1 = (a * uk_kev * 1000.0).powi(2);
        let term2 = (2.0 / 3.0) * ((fk as f64 * tot_num) / (sum_z as f64)) * plasma_energy.powi(2);
        (term1 + term2).sqrt()
    }

    /// Qminus (minimum momentum transfer, J).  Matches Java `getQminusModified`.
    fn get_qminus_modified(e_kev: f64, wak_ev: f64) -> f64 {
        let m: f64 = 9.10938356e-31;
        let c: f64 = 299_792_458.0;
        let csq = c * c;
        let e_j = e_kev * KEV_TO_JOULES;
        let wak_j = (wak_ev / 1000.0) * KEV_TO_JOULES;
        if wak_j > e_j {
            return 0.0;
        }
        let a = (e_j * (e_j + 2.0 * m * csq)).sqrt();
        let b = ((e_j - wak_j) * (e_j - wak_j + 2.0 * m * csq)).sqrt();
        ((a - b).powi(2) + (m * csq).powi(2)).sqrt() - m * csq
    }

    /// Wak: effective energy loss (eV).  Matches Java `getWak`.
    fn get_wak(e_kev: f64, wk_ev: f64, uk_ev: f64) -> f64 {
        if e_kev * 1000.0 > 3.0 * wk_ev - 2.0 * uk_ev {
            wk_ev
        } else {
            (e_kev * 1000.0 + 2.0 * uk_ev) / 3.0
        }
    }

    /// Qak: cut-off momentum transfer (eV).  Matches Java `getQak`.
    fn get_qak(e_kev: f64, wk_ev: f64, uk_ev: f64) -> f64 {
        if e_kev * 1000.0 > 3.0 * wk_ev - 2.0 * uk_ev {
            uk_ev
        } else {
            uk_ev * (e_kev * 1000.0 / (3.0 * wk_ev - 2.0 * uk_ev))
        }
    }

    /// Edash = E + Uk/1000.  Matches Java `getEdash`.
    fn get_edash(e_kev: f64, uk_ev: f64) -> f64 {
        e_kev + uk_ev / 1000.0
    }

    /// Relativistic close-collision parameter.  Matches Java `getClosea`.
    fn get_close_a(e_kev: f64) -> f64 {
        let m: f64 = 9.10938356e-31;
        let c: f64 = 299_792_458.0;
        let csq = c * c;
        let vo = e_kev * KEV_TO_JOULES;
        (vo / (vo + m * csq)).powi(2)
    }

    /// Analytical solution for the close-collision integral.
    /// Matches Java `solveCloseAnalytical`.
    fn solve_close_analytical(win_ev: f64, n: i32, energy_kev: f64) -> f64 {
        let e = energy_kev * KEV_TO_JOULES;
        let w = (win_ev / 1000.0) * KEV_TO_JOULES;
        if e <= w {
            return 0.0;
        }
        let a = Self::get_close_a(energy_kev);
        match n {
            0 => {
                -1.0 / w + 1.0 / (e - w) + ((1.0 - a) / e) * ((e - w) / w).ln() + a * w / e.powi(2)
            }
            1 => {
                w.ln()
                    + e / (e - w)
                    + (2.0 - a) * (e - w).ln()
                    + (a * w.powi(2)) / (2.0 * e.powi(2))
            }
            2 => {
                (2.0 - a) * w
                    + (2.0 * e.powi(2) - w.powi(2)) / (e - w)
                    + (3.0 - a) * e * (e - w).ln()
                    + (a * w.powi(3)) / (3.0 * e.powi(2))
            }
            _ => 0.0,
        }
    }

    /// Close-collision integral.  Matches Java `doCloseIntegral`.
    fn do_close_integral(e_kev: f64, n: i32, uk_ev: f64, qak_ev: f64) -> f64 {
        let e_dash = Self::get_edash(e_kev, uk_ev);
        let wmax = 1000.0 * e_dash / 2.0; // eV
        Self::solve_close_analytical(wmax, n, e_dash)
            - Self::solve_close_analytical(qak_ev, n, e_dash)
    }

    /// Distant-interaction GOS distribution p(W) for a shell.
    fn get_p_dis_w(
        e_kev: f64,
        w_j: f64,
        _a: f64,
        element: &str,
        shell: usize,
        _surrounding: bool,
        wk_cache: f64,
    ) -> f64 {
        let uk_kev = Self::shell_binding_kev(element, shell);
        let wak = Self::get_wak(e_kev, wk_cache, uk_kev * 1000.0);
        let wdis_ev = 3.0 * wak - 2.0 * uk_kev * 1000.0;
        let uk_j = (uk_kev / 1.0) * KEV_TO_JOULES; // uk in keV → J
        let wdis_j = (wdis_ev / 1000.0) * KEV_TO_JOULES;
        if w_j >= uk_j && w_j < wdis_j {
            (2.0 / (wdis_j - uk_j).powi(2)) * (wdis_j - w_j)
        } else {
            0.0
        }
    }

    /// Integral of W^(n-1) × p_dis(W) by trapezoid rule.  Matches Java `integrateDist`.
    #[allow(clippy::too_many_arguments)]
    fn integrate_dist(
        &self,
        e_kev: f64,
        uk_ev: f64,
        n: i32,
        shell: usize,
        element: &str,
        a: f64,
        surrounding: bool,
    ) -> f64 {
        const WBINS: usize = 1000;
        let wk = self.wk_molecule_calc(a, element, shell, surrounding);
        if wk <= 0.0 {
            return 0.0;
        }
        let wak = Self::get_wak(e_kev, wk, uk_ev);
        let wdis_ev = 3.0 * wak - 2.0 * uk_ev;
        let wdis_j = (wdis_ev / 1000.0) * KEV_TO_JOULES;
        // Java's integrateDist ignores the Uk parameter for the starting point and
        // uses getShellBinding(i, e) (coarse: K/L1/M1 by index) instead of
        // getShellBindingSubshell (fine subshell). Match that behaviour.
        let coarse_uk_kev = Self::shell_binding_kev(element, shell);
        let binding_j = coarse_uk_kev * KEV_TO_JOULES;
        let step = wdis_j / (WBINS as f64);
        if step <= 0.0 {
            return 0.0;
        }
        let mut w = binding_j;
        let mut prev_y = 0.0_f64;
        let mut sum = 0.0_f64;
        let mut count = 0;
        loop {
            let this_y = if w == 0.0 {
                0.0
            } else {
                let p = Self::get_p_dis_w(e_kev, w, a, element, shell, surrounding, wk);
                if n == 1 {
                    p
                } else {
                    w.powi(n - 1) * p
                }
            };
            if count > 0 {
                sum += step * (this_y + prev_y) / 2.0;
            }
            count += 1;
            w += step;
            prev_y = this_y;
            if w > wdis_j {
                break;
            }
        }
        sum
    }

    /// Plasmon distribution integral.  Matches Java `integrateDistPlasmon`.
    fn integrate_dist_plasmon(e_kev: f64, n: i32, wcb_ev: f64) -> f64 {
        let uk_ev = 0.0_f64;
        let e_dash = Self::get_edash(e_kev, uk_ev);
        let wmax = (e_dash / 2.0) * KEV_TO_JOULES; // J
        let w = (wcb_ev / 1000.0) * KEV_TO_JOULES; // J
        if w > 0.0 && w < wmax {
            if n == 1 {
                w.powi(0)
            } else {
                w.powi(n - 1)
            }
        } else {
            0.0
        }
    }

    /// Sum of Z×ln(I) for mean excitation energy, weighted by atoms.
    fn get_zlni(&self, surrounding: bool) -> f64 {
        let elements: Vec<String> = if surrounding {
            self.cryo_elements.iter().cloned().collect()
        } else {
            self.present_elements.iter().cloned().collect()
        };
        let db = crate::element::database::ElementDatabase::instance();
        let mut sum_z = 0.0_f64;
        let mut sum_a = 0.0_f64;
        let mut meln_i = 0.0_f64;
        for name in &elements {
            let atom_num = if surrounding {
                self.total_atoms_surr(name)
            } else {
                self.total_atoms(name)
            };
            let Some(e) = db.get(name) else { continue };
            let z = e.atomic_number() as f64;
            let a = e.atomic_weight();
            let mol_weight_fraction = (atom_num * a) / self.molecular_weight;
            sum_z += z * atom_num;
            sum_a += a * atom_num;
            let mut j = e.mean_ionisation_potential(); // eV
            if z as i32 != 1
                && z as i32 != 6
                && z as i32 != 7
                && z as i32 != 8
                && z as i32 != 9
                && z as i32 != 17
            {
                j *= 1.13;
            }
            meln_i += mol_weight_fraction * (z / a) * j.ln();
        }
        if sum_a == 0.0 {
            return 0.0;
        }
        meln_i /= sum_z / sum_a;
        sum_z * meln_i
    }

    /// Check if Sturnheimer `a` is consistent with mean excitation energy.
    /// Returns `[ZlnI, RHS]`.  Matches Java `checkMeanI`.
    fn check_mean_i(&self, _e_kev: f64, a: f64, surrounding: bool) -> [f64; 2] {
        const FCB_CUTOFF: f64 = 0.0;
        const SUBSHELLS: [i32; 10] = [2, 2, 2, 4, 2, 2, 4, 10, 2, 30];
        let zlni = self.get_zlni(surrounding);
        let elements: Vec<String> = if surrounding {
            self.cryo_elements.iter().cloned().collect()
        } else {
            self.present_elements.iter().cloned().collect()
        };
        let db = crate::element::database::ElementDatabase::instance();
        let mut sum_fcb = 0.0_f64;
        let mut sum_fk_ln_wk = 0.0_f64;
        for name in &elements {
            let atom_num = if surrounding {
                self.total_atoms_surr(name)
            } else {
                self.total_atoms(name)
            };
            let Some(e) = db.get(name) else { continue };
            let [valence, num_inner] = Self::num_valence_electrons_subshells_z(e.atomic_number());
            for i in 0..=(num_inner as usize) {
                let mut fk = SUBSHELLS.get(i).copied().unwrap_or(0);
                if i == num_inner as usize {
                    let uk_ev = Self::shell_binding_subshell_kev(name, i) * 1000.0;
                    if uk_ev > FCB_CUTOFF {
                        fk = valence;
                    } else {
                        sum_fcb += valence as f64;
                        fk = 0;
                    }
                }
                let wk = self.wk_molecule_calc(a, name, i, surrounding);
                if wk > 0.0 {
                    sum_fk_ln_wk += atom_num * fk as f64 * wk.ln();
                }
            }
        }
        let wcb = self.calc_wcb_all(surrounding);
        let fcb_ln_wcb = if wcb > 0.0 { sum_fcb * wcb.ln() } else { 0.0 };
        let rhs = fcb_ln_wcb + sum_fk_ln_wk;
        [zlni, rhs]
    }

    /// Find Sturnheimer `a` by bisection on `checkMeanI`.  Matches Java `getSturnheimera`.
    fn get_sturnheimer_a(&self, e_kev: f64, surrounding: bool) -> f64 {
        let mut a = 1.5_f64;
        let check = self.check_mean_i(e_kev, a, surrounding);
        if check[1] > check[0] {
            // a too high, decrease
            while a > 0.1 {
                a -= 0.1;
                let c = self.check_mean_i(e_kev, a, surrounding);
                if c[1] <= c[0] {
                    break;
                }
            }
        } else {
            // a too low, increase
            while a < 3.0 {
                let prev = self.check_mean_i(e_kev, a, surrounding);
                a += 0.1;
                let c = self.check_mean_i(e_kev, a, surrounding);
                if c[1] >= c[0] {
                    // interpolate back
                    let frac = (c[1] - c[0]) / (c[1] - prev[1]);
                    a -= 0.1 * frac;
                    break;
                }
            }
        }
        a
    }

    /// Populate GOS inelastic tables and return lambda (nm).
    /// `n=0` → MFP, `n=1` → stopping power, `n=2` → straggling.
    /// Matches Java `populateGOSInel`.
    #[allow(clippy::needless_range_loop)]
    pub fn populate_gos_inel(&mut self, e_kev: f64, n: i32, a: f64, surrounding: bool) -> f64 {
        const FCB_CUTOFF: f64 = 0.0;
        const SUBSHELLS: [i32; 10] = [2, 2, 2, 4, 2, 2, 4, 10, 2, 30];
        const MAX_SHELLS: usize = 10;

        let e_charge: f64 = 4.80320425e-10;
        let m: f64 = 9.10938356e-31;
        let c: f64 = 299_792_458.0;
        let csq = c * c;
        let vo = e_kev * KEV_TO_JOULES;
        let beta_sq = 1.0 - (m * csq / (vo + m * csq)).powi(2);
        let v_sq = beta_sq * csq;
        let constant = 2.0 * std::f64::consts::PI * (e_charge.powi(4) / 1.0e18) / (m * v_sq);

        let elements: Vec<String> = if surrounding {
            self.cryo_elements.iter().cloned().collect()
        } else {
            self.present_elements.iter().cloned().collect()
        };
        let db = crate::element::database::ElementDatabase::instance();

        let mut check_sum = 0.0_f64;
        let mut sum_fcb = 0.0_f64;
        let mut new_map: HashMap<String, Vec<[f64; 4]>> = HashMap::new();

        for name in &elements {
            let atom_num = if surrounding {
                self.total_atoms_surr(name)
            } else {
                self.total_atoms(name)
            };
            if atom_num <= 0.0 {
                continue;
            }
            let Some(e) = db.get(name) else { continue };
            let [valence, num_inner] = Self::num_valence_electrons_subshells_z(e.atomic_number());
            let mut inel_shell = vec![[0.0_f64; 4]; MAX_SHELLS];

            for i in 0..=(num_inner as usize) {
                let mut fk = SUBSHELLS.get(i).copied().unwrap_or(0);
                if i == num_inner as usize {
                    fk = valence;
                }
                let wk = self.wk_molecule_calc(a, name, i, surrounding);
                if wk <= 0.0 {
                    continue;
                }
                let uk_ev = Self::shell_binding_subshell_kev(name, i) * 1000.0;
                if uk_ev < FCB_CUTOFF {
                    sum_fcb += fk as f64 * atom_num;
                    fk = 0;
                }
                let wak = Self::get_wak(e_kev, wk, uk_ev);
                let qak = Self::get_qak(e_kev, wk, uk_ev);
                let qminus = Self::get_qminus_modified(e_kev, wak);
                let integral = self.integrate_dist(e_kev, uk_ev, n, i, name, a, surrounding);

                if e_kev * 1000.0 > uk_ev && uk_ev >= FCB_CUTOFF {
                    if qak > 0.0 && qminus > 0.0 {
                        let qak_j = (qak / 1000.0) * KEV_TO_JOULES;
                        let log_term =
                            (qak_j / qminus) * ((qminus + 2.0 * m * csq) / (qak_j + 2.0 * m * csq));
                        inel_shell[i][0] =
                            (1.0e18 * constant * (atom_num * fk as f64 * log_term.ln() * integral))
                                / (self.cell_volume / 1000.0);
                    }
                    inel_shell[i][1] = (1.0e18
                        * constant
                        * (atom_num
                            * fk as f64
                            * ((1.0 / (1.0 - beta_sq)).ln() - beta_sq)
                            * integral))
                        / (self.cell_volume / 1000.0);
                    let close_integral = Self::do_close_integral(e_kev, n, uk_ev, qak);
                    inel_shell[i][2] =
                        (1.0e18 * constant * (atom_num * fk as f64 * close_integral))
                            / (self.cell_volume / 1000.0);
                    inel_shell[i][3] = inel_shell[i][0] + inel_shell[i][1] + inel_shell[i][2];
                    check_sum += inel_shell[i][3];
                }
            }
            new_map.insert(name.clone(), inel_shell);
        }

        // Plasmon contribution
        let wcb = self.calc_wcb_all(surrounding);
        let qcb = wcb;
        let qminus_cb = Self::get_qminus_modified(e_kev, wcb);
        let integral_cb = Self::integrate_dist_plasmon(e_kev, n, wcb);
        let mut cb = [0.0_f64; 4];
        if wcb > 0.0 && e_kev * 1000.0 > wcb {
            let qcb_j = (qcb / 1000.0) * KEV_TO_JOULES;
            if qcb_j > 0.0 && qminus_cb > 0.0 {
                let log_term =
                    (qcb_j / qminus_cb) * ((qminus_cb + 2.0 * m * csq) / (qcb_j + 2.0 * m * csq));
                cb[0] = (1.0e18 * constant * (sum_fcb * log_term.ln() * integral_cb))
                    / (self.cell_volume / 1000.0);
            }
            cb[1] = (1.0e18
                * constant
                * (sum_fcb * ((1.0 / (1.0 - beta_sq)).ln() - beta_sq) * integral_cb))
                / (self.cell_volume / 1000.0);
            let close_integral_cb = Self::do_close_integral(e_kev, n, 0.0, qcb);
            cb[2] = (1.0e18 * constant * sum_fcb * close_integral_cb) / (self.cell_volume / 1000.0);
            cb[3] = cb[0] + cb[1] + cb[2];
            check_sum += cb[3];
        }

        if surrounding {
            self.gos_inelastic_surrounding = new_map;
            self.cb_inel_surrounding = cb;
        } else {
            self.gos_inelastic = new_map;
            self.cb_inel = cb;
        }

        if check_sum > 0.0 {
            1.0 / check_sum
        } else {
            0.0
        }
    }

    /// GOS inelastic MFP (nm).  Matches Java `getGOSInel`.
    pub fn calc_gos_inel(&mut self, surrounding: bool, e_kev: f64) -> f64 {
        let a = self.get_sturnheimer_a(e_kev, surrounding);
        self.sturnheimer_adjustment = a;
        self.populate_gos_inel(e_kev, 0, a, surrounding)
    }

    /// GOS inner-shell lambda (nm).  Matches Java `getGOSInnerLambda`.
    pub fn calc_gos_inner_lambda(&self, surrounding: bool) -> f64 {
        let map = if surrounding {
            &self.gos_inelastic_surrounding
        } else {
            &self.gos_inelastic
        };
        let db = crate::element::database::ElementDatabase::instance();
        let mut lambda_sum = 0.0_f64;
        for (name, shells) in map {
            let z = db.get(name).map(|e| e.atomic_number()).unwrap_or(0);
            let [_, num_inner] = Self::num_valence_electrons_subshells_z(z);
            for i in 0..(num_inner as usize) {
                if i < shells.len() {
                    lambda_sum += shells[i][3];
                }
            }
        }
        if lambda_sum > 0.0 {
            1.0 / lambda_sum
        } else {
            0.0
        }
    }

    /// GOS outer-shell lambda (nm).  Matches Java `getGOSOuterLambda`.
    pub fn calc_gos_outer_lambda(&self, surrounding: bool) -> f64 {
        let map = if surrounding {
            &self.gos_inelastic_surrounding
        } else {
            &self.gos_inelastic
        };
        let db = crate::element::database::ElementDatabase::instance();
        let mut lambda_sum = 0.0_f64;
        for (name, shells) in map {
            let z = db.get(name).map(|e| e.atomic_number()).unwrap_or(0);
            let [_, num_inner] = Self::num_valence_electrons_subshells_z(z);
            let outer = num_inner as usize;
            if outer < shells.len() {
                lambda_sum += shells[outer][3];
            }
        }
        if lambda_sum > 0.0 {
            1.0 / lambda_sum
        } else {
            0.0
        }
    }

    /// GOS shell ionisation probabilities (CDF per element per shell).
    /// Matches Java `getGOSShellProbs`.
    pub fn calc_gos_shell_probs(
        &self,
        surrounding: bool,
        tot_lambda: f64,
    ) -> HashMap<String, Vec<f64>> {
        let tot_inel = 1.0 / tot_lambda;
        let map = if surrounding {
            &self.gos_inelastic_surrounding
        } else {
            &self.gos_inelastic
        };
        let mut probs: HashMap<String, Vec<f64>> = HashMap::new();
        let mut running_sum = 0.0_f64;
        for (name, shells) in map {
            let mut shell_probs = vec![0.0_f64; 9];
            for (i, sh) in shells.iter().enumerate().take(9) {
                let frac = if tot_inel > 0.0 {
                    sh[3] / tot_inel
                } else {
                    0.0
                };
                running_sum += frac;
                shell_probs[i] = running_sum;
            }
            probs.insert(name.clone(), shell_probs);
        }
        probs
    }

    /// GOS outer-shell ionisation probabilities (CDF per element).
    /// Matches Java `getGOSOuterShellProbs`.
    pub fn calc_gos_outer_shell_probs(
        &self,
        surrounding: bool,
        tot_lambda: f64,
    ) -> HashMap<String, f64> {
        let tot_inel = 1.0 / tot_lambda;
        let map = if surrounding {
            &self.gos_inelastic_surrounding
        } else {
            &self.gos_inelastic
        };
        let db = crate::element::database::ElementDatabase::instance();
        let mut probs: HashMap<String, f64> = HashMap::new();
        let mut running_sum = 0.0_f64;
        for (name, shells) in map {
            let z = db.get(name).map(|e| e.atomic_number()).unwrap_or(0);
            let [_, num_inner] = Self::num_valence_electrons_subshells_z(z);
            let outer = num_inner as usize;
            let frac = if tot_inel > 0.0 && outer < shells.len() {
                shells[outer][3] / tot_inel
            } else {
                0.0
            };
            running_sum += frac;
            probs.insert(name.clone(), running_sum);
        }
        probs
    }

    /// Returns raw GOSinelastic state as flat probability map.
    /// Matches Java `getGOSVariable` (returns the map itself).
    pub fn calc_gos_variable(&self, surrounding: bool) -> HashMap<String, Vec<f64>> {
        let map = if surrounding {
            &self.gos_inelastic_surrounding
        } else {
            &self.gos_inelastic
        };
        map.iter()
            .map(|(k, v)| {
                let flat: Vec<f64> = v.iter().flat_map(|sh| sh.iter().copied()).collect();
                (k.clone(), flat)
            })
            .collect()
    }

    /// Plasmon (conduction-band) total cross-section (nm⁻¹).
    /// Matches Java `getPlasmonVariable` → returns cbInel[3].
    pub fn calc_plasmon_variable(&self, surrounding: bool) -> f64 {
        if surrounding {
            self.cb_inel_surrounding[3]
        } else {
            self.cb_inel[3]
        }
    }
}
