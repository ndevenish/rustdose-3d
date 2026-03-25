use std::collections::HashMap;
use std::sync::OnceLock;

/// Atomic mass unit in grams.
const ATOMIC_MASS_UNIT: f64 = 1.66e-24;

/// LJ_1 correction factor for light elements (from Fortran code).
const LJ_1: f64 = 1.160;

/// LJ_2 correction factor for light elements (from Fortran code).
const LJ_2: f64 = 1.41;

/// Absorption edge tolerance in keV (reserved for phase-6 edge-handling).
#[allow(dead_code)]
const ABSORPTION_EDGE_TOLERANCE: f64 = 0.001;

/// Number of polynomial expansion terms.
const POLYNOMIAL_EXPANSION: usize = 4;

/// Light atom threshold (Z <= 29 treated as light).
pub const LIGHT_ATOM_MAX_NUM: i32 = 29;

/// Database field indices into the tab-separated MuCalcConstants.txt.
/// Each variant holds its column index (0-based).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DatabaseField {
    EdgeK,                   // 2
    EdgeL,                   // 3
    EdgeM,                   // 4
    KCoeff0,                 // 5
    KCoeff1,                 // 6
    KCoeff2,                 // 7
    KCoeff3,                 // 8
    LCoeff0,                 // 9
    LCoeff1,                 // 10
    LCoeff2,                 // 11
    LCoeff3,                 // 12
    MCoeff0,                 // 13
    MCoeff1,                 // 14
    MCoeff2,                 // 15
    MCoeff3,                 // 16
    NCoeff0,                 // 17
    NCoeff1,                 // 18
    NCoeff2,                 // 19
    NCoeff3,                 // 20
    AtomicWeight,            // 23
    CoherentCoeff0,          // 24
    CoherentCoeff1,          // 25
    CoherentCoeff2,          // 26
    CoherentCoeff3,          // 27
    IncoherentCoeff0,        // 28
    IncoherentCoeff1,        // 29
    IncoherentCoeff2,        // 30
    IncoherentCoeff3,        // 31
    L2,                      // 36
    L3,                      // 37
    FluorescenceYieldK,      // 39
    FluorescenceYieldL1,     // 40
    FluorescenceYieldL2,     // 41
    FluorescenceYieldL3,     // 42
    KEdgeRatio,              // 43
    L1EdgeRatio,             // 44
    L2EdgeRatio,             // 45
    L3EdgeRatio,             // 46
    KFlAverage,              // 47
    LFlAverage,              // 48
    EdgeM2,                  // 49
    EdgeM3,                  // 50
    EdgeM4,                  // 51
    EdgeM5,                  // 52
    M1EdgeRatio,             // 53
    M2EdgeRatio,             // 54
    M3EdgeRatio,             // 55
    M4EdgeRatio,             // 56
    M5EdgeRatio,             // 57
    EminLow,                 // 58
    EmaxLow,                 // 59
    BKlow,                   // 60
    CKlow,                   // 61
    EminHigh,                // 62
    BKhigh,                  // 63
    CKhigh,                  // 64
    El50,                    // 65
    El100,                   // 66
    El150,                   // 67
    El200,                   // 68
    El250,                   // 69
    El300,                   // 70
    MeanIonisationPotential, // 71
    L1,                      // 72
    M1,                      // 73
    N1,                      // 74
    FluorescenceYieldM1,     // 75
    FluorescenceYieldM2,     // 76
    FluorescenceYieldM3,     // 77
    FluorescenceYieldM4,     // 78
    FluorescenceYieldM5,     // 79
}

impl DatabaseField {
    fn column(self) -> usize {
        match self {
            Self::EdgeK => 2,
            Self::EdgeL => 3,
            Self::EdgeM => 4,
            Self::KCoeff0 => 5,
            Self::KCoeff1 => 6,
            Self::KCoeff2 => 7,
            Self::KCoeff3 => 8,
            Self::LCoeff0 => 9,
            Self::LCoeff1 => 10,
            Self::LCoeff2 => 11,
            Self::LCoeff3 => 12,
            Self::MCoeff0 => 13,
            Self::MCoeff1 => 14,
            Self::MCoeff2 => 15,
            Self::MCoeff3 => 16,
            Self::NCoeff0 => 17,
            Self::NCoeff1 => 18,
            Self::NCoeff2 => 19,
            Self::NCoeff3 => 20,
            Self::AtomicWeight => 23,
            Self::CoherentCoeff0 => 24,
            Self::CoherentCoeff1 => 25,
            Self::CoherentCoeff2 => 26,
            Self::CoherentCoeff3 => 27,
            Self::IncoherentCoeff0 => 28,
            Self::IncoherentCoeff1 => 29,
            Self::IncoherentCoeff2 => 30,
            Self::IncoherentCoeff3 => 31,
            Self::L2 => 36,
            Self::L3 => 37,
            Self::FluorescenceYieldK => 39,
            Self::FluorescenceYieldL1 => 40,
            Self::FluorescenceYieldL2 => 41,
            Self::FluorescenceYieldL3 => 42,
            Self::KEdgeRatio => 43,
            Self::L1EdgeRatio => 44,
            Self::L2EdgeRatio => 45,
            Self::L3EdgeRatio => 46,
            Self::KFlAverage => 47,
            Self::LFlAverage => 48,
            Self::EdgeM2 => 49,
            Self::EdgeM3 => 50,
            Self::EdgeM4 => 51,
            Self::EdgeM5 => 52,
            Self::M1EdgeRatio => 53,
            Self::M2EdgeRatio => 54,
            Self::M3EdgeRatio => 55,
            Self::M4EdgeRatio => 56,
            Self::M5EdgeRatio => 57,
            Self::EminLow => 58,
            Self::EmaxLow => 59,
            Self::BKlow => 60,
            Self::CKlow => 61,
            Self::EminHigh => 62,
            Self::BKhigh => 63,
            Self::CKhigh => 64,
            Self::El50 => 65,
            Self::El100 => 66,
            Self::El150 => 67,
            Self::El200 => 68,
            Self::El250 => 69,
            Self::El300 => 70,
            Self::MeanIonisationPotential => 71,
            Self::L1 => 72,
            Self::M1 => 73,
            Self::N1 => 74,
            Self::FluorescenceYieldM1 => 75,
            Self::FluorescenceYieldM2 => 76,
            Self::FluorescenceYieldM3 => 77,
            Self::FluorescenceYieldM4 => 78,
            Self::FluorescenceYieldM5 => 79,
        }
    }

    fn all() -> &'static [DatabaseField] {
        use DatabaseField::*;
        &[
            EdgeK,
            EdgeL,
            EdgeM,
            KCoeff0,
            KCoeff1,
            KCoeff2,
            KCoeff3,
            LCoeff0,
            LCoeff1,
            LCoeff2,
            LCoeff3,
            MCoeff0,
            MCoeff1,
            MCoeff2,
            MCoeff3,
            NCoeff0,
            NCoeff1,
            NCoeff2,
            NCoeff3,
            AtomicWeight,
            CoherentCoeff0,
            CoherentCoeff1,
            CoherentCoeff2,
            CoherentCoeff3,
            IncoherentCoeff0,
            IncoherentCoeff1,
            IncoherentCoeff2,
            IncoherentCoeff3,
            L2,
            L3,
            FluorescenceYieldK,
            FluorescenceYieldL1,
            FluorescenceYieldL2,
            FluorescenceYieldL3,
            KEdgeRatio,
            L1EdgeRatio,
            L2EdgeRatio,
            L3EdgeRatio,
            KFlAverage,
            LFlAverage,
            EdgeM2,
            EdgeM3,
            EdgeM4,
            EdgeM5,
            M1EdgeRatio,
            M2EdgeRatio,
            M3EdgeRatio,
            M4EdgeRatio,
            M5EdgeRatio,
            EminLow,
            EmaxLow,
            BKlow,
            CKlow,
            EminHigh,
            BKhigh,
            CKhigh,
            El50,
            El100,
            El150,
            El200,
            El250,
            El300,
            MeanIonisationPotential,
            L1,
            M1,
            N1,
            FluorescenceYieldM1,
            FluorescenceYieldM2,
            FluorescenceYieldM3,
            FluorescenceYieldM4,
            FluorescenceYieldM5,
        ]
    }
}

/// Absorption edge type for polynomial coefficient lookup.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum AbsorptionEdge {
    K,
    L,
    M,
    N,
    /// Coherent scattering polynomial
    C,
    /// Incoherent scattering polynomial
    I,
}

/// Cross-section types returned by absorption coefficient calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CrossSection {
    Photoelectric,
    Coherent,
    Compton,
    Total,
}

/// An element with X-ray cross-section data from MuCalcConstants.txt.
#[derive(Debug, Clone)]
pub struct Element {
    name: String,
    atomic_number: i32,
    data: HashMap<DatabaseField, Option<f64>>,
    coefficients: HashMap<AbsorptionEdge, [f64; POLYNOMIAL_EXPANSION]>,
    // Ionisation probabilities (computed from edge ratios)
    prob_k_ionisation: f64,
    prob_l1_ionisation: f64,
    prob_l2_ionisation: f64,
    prob_l3_ionisation: f64,
    prob_m1_ionisation: f64,
    prob_m2_ionisation: f64,
    prob_m3_ionisation: f64,
    prob_m4_ionisation: f64,
    prob_m5_ionisation: f64,
}

impl Element {
    fn new(name: String, atomic_number: i32, data: HashMap<DatabaseField, Option<f64>>) -> Self {
        let coefficients = Self::build_coefficients(&data);

        let mut elem = Element {
            name,
            atomic_number,
            data,
            coefficients,
            prob_k_ionisation: 0.0,
            prob_l1_ionisation: 0.0,
            prob_l2_ionisation: 0.0,
            prob_l3_ionisation: 0.0,
            prob_m1_ionisation: 0.0,
            prob_m2_ionisation: 0.0,
            prob_m3_ionisation: 0.0,
            prob_m4_ionisation: 0.0,
            prob_m5_ionisation: 0.0,
        };
        elem.compute_edge_ratios();
        elem
    }

    fn build_coefficients(
        data: &HashMap<DatabaseField, Option<f64>>,
    ) -> HashMap<AbsorptionEdge, [f64; POLYNOMIAL_EXPANSION]> {
        let mut map = HashMap::new();

        let get = |f: DatabaseField| -> f64 { data.get(&f).and_then(|v| *v).unwrap_or(0.0) };

        map.insert(
            AbsorptionEdge::K,
            [
                get(DatabaseField::KCoeff0),
                get(DatabaseField::KCoeff1),
                get(DatabaseField::KCoeff2),
                get(DatabaseField::KCoeff3),
            ],
        );
        map.insert(
            AbsorptionEdge::L,
            [
                get(DatabaseField::LCoeff0),
                get(DatabaseField::LCoeff1),
                get(DatabaseField::LCoeff2),
                get(DatabaseField::LCoeff3),
            ],
        );
        map.insert(
            AbsorptionEdge::M,
            [
                get(DatabaseField::MCoeff0),
                get(DatabaseField::MCoeff1),
                get(DatabaseField::MCoeff2),
                get(DatabaseField::MCoeff3),
            ],
        );
        map.insert(
            AbsorptionEdge::N,
            [
                get(DatabaseField::NCoeff0),
                get(DatabaseField::NCoeff1),
                get(DatabaseField::NCoeff2),
                get(DatabaseField::NCoeff3),
            ],
        );
        map.insert(
            AbsorptionEdge::C,
            [
                get(DatabaseField::CoherentCoeff0),
                get(DatabaseField::CoherentCoeff1),
                get(DatabaseField::CoherentCoeff2),
                get(DatabaseField::CoherentCoeff3),
            ],
        );
        map.insert(
            AbsorptionEdge::I,
            [
                get(DatabaseField::IncoherentCoeff0),
                get(DatabaseField::IncoherentCoeff1),
                get(DatabaseField::IncoherentCoeff2),
                get(DatabaseField::IncoherentCoeff3),
            ],
        );

        map
    }

    /// Compute ionisation probabilities from edge ratios.
    /// Java: `EdgeRatio()` method.
    fn compute_edge_ratios(&mut self) {
        let compute = |ratio: Option<f64>| -> f64 {
            match ratio {
                Some(r) if r != 0.0 => {
                    let p = 1.0 - 1.0 / r;
                    if p == f64::NEG_INFINITY {
                        1.0
                    } else {
                        p
                    }
                }
                _ => 1.0,
            }
        };

        self.prob_k_ionisation = compute(self.get_field(DatabaseField::KEdgeRatio));
        self.prob_l1_ionisation = compute(self.get_field(DatabaseField::L1EdgeRatio));
        self.prob_l2_ionisation = compute(self.get_field(DatabaseField::L2EdgeRatio));
        self.prob_l3_ionisation = compute(self.get_field(DatabaseField::L3EdgeRatio));
        self.prob_m1_ionisation = compute(self.get_field(DatabaseField::M1EdgeRatio));
        self.prob_m2_ionisation = compute(self.get_field(DatabaseField::M2EdgeRatio));
        self.prob_m3_ionisation = compute(self.get_field(DatabaseField::M3EdgeRatio));
        self.prob_m4_ionisation = compute(self.get_field(DatabaseField::M4EdgeRatio));
        self.prob_m5_ionisation = compute(self.get_field(DatabaseField::M5EdgeRatio));
    }

    /// McMaster polynomial expansion: exp(sum of c_i * ln(E)^i).
    fn bax_for_edge(&self, energy: f64, edge: AbsorptionEdge) -> f64 {
        let coeffs = &self.coefficients[&edge];
        let mut sum = 0.0;

        for (i, &c) in coeffs[..POLYNOMIAL_EXPANSION].iter().enumerate() {
            if energy == 1.0 {
                // Match Java behavior: at E=1keV, ln(1)=0, so pow(0,i)
                // would be 0 for i>0. Java adds coeffs[i] directly.
                sum += c;
            } else {
                sum += c * energy.ln().powi(i as i32);
            }
        }

        sum.exp()
    }

    /// Get absorption coefficients (barns/atom) for a given X-ray energy (keV).
    /// Returns photoelectric, coherent, compton, and total cross-sections.
    pub fn get_abs_coefficients(&self, energy: f64) -> HashMap<CrossSection, f64> {
        let photoelectric = self.get_photoelectric_xs(energy);

        let elastic = if self.get_field(DatabaseField::CoherentCoeff0).unwrap_or(0.0) != 0.0 {
            self.bax_for_edge(energy, AbsorptionEdge::C)
        } else {
            0.0
        };

        let compton = if self
            .get_field(DatabaseField::IncoherentCoeff0)
            .unwrap_or(0.0)
            != 0.0
        {
            self.bax_for_edge(energy, AbsorptionEdge::I)
        } else {
            0.0
        };

        let total = photoelectric + elastic + compton;

        let mut results = HashMap::new();
        results.insert(CrossSection::Photoelectric, photoelectric);
        results.insert(CrossSection::Coherent, elastic);
        results.insert(CrossSection::Compton, compton);
        results.insert(CrossSection::Total, total);
        results
    }

    /// Determine photoelectric cross-section for a given energy.
    fn get_photoelectric_xs(&self, energy: f64) -> f64 {
        let edge_k = self
            .get_field(DatabaseField::EdgeK)
            .expect("K edge undefined");
        let edge_l = self.get_field(DatabaseField::EdgeL);
        let _edge_l2 = self.get_field(DatabaseField::L2);
        let edge_l3 = self.get_field(DatabaseField::L3);
        let edge_m = self.get_field(DatabaseField::EdgeM);

        // Select the appropriate edge polynomial
        let mut photoelectric = if energy > edge_k || edge_l.is_none() {
            self.bax_for_edge(energy, AbsorptionEdge::K)
        } else if energy > edge_l3.unwrap_or(0.0) || edge_m.is_none() {
            self.bax_for_edge(energy, AbsorptionEdge::L)
        } else if energy > edge_m.unwrap_or(0.0) {
            self.bax_for_edge(energy, AbsorptionEdge::M)
        } else {
            self.bax_for_edge(energy, AbsorptionEdge::N)
        };

        // L-edge correction (from Fortran): correct McMaster L1 edge
        // using edge jumps for correct cross-sections.
        let l3 = self.get_field(DatabaseField::L3).unwrap_or(0.0);
        let l2 = self.get_field(DatabaseField::L2).unwrap_or(0.0);
        let edge_l_val = edge_l.unwrap_or(0.0);

        if energy > l3 && energy < l2 {
            photoelectric /= LJ_1 * LJ_2;
        }
        if energy > l2 && energy < edge_l_val {
            photoelectric /= LJ_1;
        }

        photoelectric
    }

    // ---- Public getters ----

    fn get_field(&self, field: DatabaseField) -> Option<f64> {
        self.data.get(&field).and_then(|v| *v)
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn atomic_number(&self) -> i32 {
        self.atomic_number
    }

    pub fn atomic_weight(&self) -> f64 {
        self.get_field(DatabaseField::AtomicWeight).unwrap_or(0.0)
    }

    pub fn atomic_weight_in_grams(&self) -> f64 {
        self.atomic_weight() * ATOMIC_MASS_UNIT
    }

    pub fn k_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::EdgeK)
    }
    pub fn l1_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::EdgeL)
    }
    pub fn l2_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::L2)
    }
    pub fn l3_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::L3)
    }
    pub fn l1_binding(&self) -> Option<f64> {
        self.get_field(DatabaseField::L1)
    }
    pub fn m1_binding(&self) -> Option<f64> {
        self.get_field(DatabaseField::M1)
    }
    pub fn n1_binding(&self) -> Option<f64> {
        self.get_field(DatabaseField::N1)
    }

    pub fn m1_edge(&self) -> f64 {
        self.get_field(DatabaseField::EdgeM).unwrap_or(0.0)
    }
    pub fn m2_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::EdgeM2)
    }
    pub fn m3_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::EdgeM3)
    }
    pub fn m4_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::EdgeM4)
    }
    pub fn m5_edge(&self) -> Option<f64> {
        self.get_field(DatabaseField::EdgeM5)
    }

    pub fn k_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldK)
    }
    pub fn l1_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldL1)
    }
    pub fn l2_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldL2)
    }
    pub fn l3_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldL3)
    }
    pub fn m1_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldM1)
    }
    pub fn m2_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldM2)
    }
    pub fn m3_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldM3)
    }
    pub fn m4_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldM4)
    }
    pub fn m5_fluorescence_yield(&self) -> Option<f64> {
        self.get_field(DatabaseField::FluorescenceYieldM5)
    }

    pub fn k_ionisation_prob(&self) -> f64 {
        self.prob_k_ionisation
    }
    pub fn l1_ionisation_prob(&self) -> f64 {
        self.prob_l1_ionisation
    }
    pub fn l2_ionisation_prob(&self) -> f64 {
        self.prob_l2_ionisation
    }
    pub fn l3_ionisation_prob(&self) -> f64 {
        self.prob_l3_ionisation
    }
    pub fn m1_ionisation_prob(&self) -> f64 {
        self.prob_m1_ionisation
    }
    pub fn m2_ionisation_prob(&self) -> f64 {
        self.prob_m2_ionisation
    }
    pub fn m3_ionisation_prob(&self) -> f64 {
        self.prob_m3_ionisation
    }
    pub fn m4_ionisation_prob(&self) -> f64 {
        self.prob_m4_ionisation
    }
    pub fn m5_ionisation_prob(&self) -> f64 {
        self.prob_m5_ionisation
    }

    pub fn k_fl_average(&self) -> Option<f64> {
        self.get_field(DatabaseField::KFlAverage)
    }
    pub fn l_fl_average(&self) -> Option<f64> {
        self.get_field(DatabaseField::LFlAverage)
    }

    /// ELSEPA elastic cross sections at 50, 100, 150, 200, 250, 300 keV.
    pub fn elsepa_coefficients(&self) -> [f64; 6] {
        [
            self.get_field(DatabaseField::El50).unwrap_or(0.0),
            self.get_field(DatabaseField::El100).unwrap_or(0.0),
            self.get_field(DatabaseField::El150).unwrap_or(0.0),
            self.get_field(DatabaseField::El200).unwrap_or(0.0),
            self.get_field(DatabaseField::El250).unwrap_or(0.0),
            self.get_field(DatabaseField::El300).unwrap_or(0.0),
        ]
    }

    pub fn emin_low(&self) -> Option<f64> {
        self.get_field(DatabaseField::EminLow)
    }
    pub fn emax_low(&self) -> Option<f64> {
        self.get_field(DatabaseField::EmaxLow)
    }
    pub fn bk_low(&self) -> Option<f64> {
        self.get_field(DatabaseField::BKlow)
    }
    pub fn ck_low(&self) -> Option<f64> {
        self.get_field(DatabaseField::CKlow)
    }
    pub fn emin_high(&self) -> Option<f64> {
        self.get_field(DatabaseField::EminHigh)
    }
    pub fn bk_high(&self) -> Option<f64> {
        self.get_field(DatabaseField::BKhigh)
    }
    pub fn ck_high(&self) -> Option<f64> {
        self.get_field(DatabaseField::CKhigh)
    }
    pub fn mean_ionisation_potential(&self) -> f64 {
        self.get_field(DatabaseField::MeanIonisationPotential)
            .unwrap_or(0.0)
    }
}

/// Singleton X-ray element database loaded from MuCalcConstants.txt.
pub struct ElementDatabase {
    by_name: HashMap<String, Element>,
    by_number: HashMap<i32, Element>,
}

static INSTANCE: OnceLock<ElementDatabase> = OnceLock::new();

/// Embedded MuCalcConstants.txt data.
const MUCALC_DATA: &str = include_str!("../constants_data/MuCalcConstants.txt");

impl ElementDatabase {
    fn load() -> Self {
        let mut by_name = HashMap::new();
        let mut by_number = HashMap::new();

        for line in MUCALC_DATA.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let components: Vec<&str> = line.split('\t').collect();
            if components.len() < 2 {
                continue;
            }

            let atomic_number: i32 = match components[0].parse() {
                Ok(n) => n,
                Err(_) => continue,
            };
            let name = components[1].to_string();

            let mut data = HashMap::new();
            for &field in DatabaseField::all() {
                let col = field.column();
                let value = if col < components.len() && !components[col].is_empty() {
                    components[col].parse::<f64>().ok()
                } else {
                    None
                };
                data.insert(field, value);
            }

            let elem = Element::new(name.clone(), atomic_number, data);
            by_name.insert(name.to_lowercase(), elem.clone());
            by_number.insert(atomic_number, elem);
        }

        ElementDatabase { by_name, by_number }
    }

    /// Get the singleton instance of the element database.
    pub fn instance() -> &'static ElementDatabase {
        INSTANCE.get_or_init(Self::load)
    }

    /// Look up element by atomic number.
    pub fn get_by_number(&self, z: i32) -> Option<&Element> {
        self.by_number.get(&z)
    }

    /// Look up element by name (case-insensitive).
    pub fn get(&self, name: &str) -> Option<&Element> {
        self.by_name.get(&name.to_lowercase())
    }

    /// Number of elements in the database.
    pub fn len(&self) -> usize {
        self.by_number.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_number.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_database_loads() {
        let db = ElementDatabase::instance();
        assert!(
            db.len() >= 80,
            "Should have at least 80 elements, got {}",
            db.len()
        );
    }

    #[test]
    fn test_lookup_by_name_and_number() {
        let db = ElementDatabase::instance();
        let fe = db.get("Fe").unwrap();
        assert_eq!(fe.atomic_number(), 26);
        assert_eq!(fe.name(), "FE");

        let fe2 = db.get_by_number(26).unwrap();
        assert_eq!(fe2.name(), "FE");

        // Case insensitive
        let c = db.get("c").unwrap();
        assert_eq!(c.atomic_number(), 6);
    }

    #[test]
    fn test_hydrogen() {
        let db = ElementDatabase::instance();
        let h = db.get_by_number(1).unwrap();
        assert_eq!(h.name(), "H");
        assert!((h.atomic_weight() - 1.008).abs() < 0.01);
        assert!(h.k_edge().unwrap() > 0.0);
    }

    #[test]
    fn test_carbon_cross_sections() {
        let db = ElementDatabase::instance();
        let c = db.get_by_number(6).unwrap();
        assert_eq!(c.name(), "C");
        assert!((c.atomic_weight() - 12.01).abs() < 0.01);

        // Test cross-section calculation at 12.4 keV
        let xs = c.get_abs_coefficients(12.4);
        let pe = xs[&CrossSection::Photoelectric];
        let total = xs[&CrossSection::Total];

        // Photoelectric should be positive and less than total
        assert!(pe > 0.0, "Photoelectric XS should be positive: {}", pe);
        assert!(total > pe, "Total should exceed photoelectric");
    }

    #[test]
    fn test_iron_cross_sections() {
        let db = ElementDatabase::instance();
        let fe = db.get_by_number(26).unwrap();

        // Iron K-edge is ~7.112 keV
        let k_edge = fe.k_edge().unwrap();
        assert!((k_edge - 7.112).abs() < 0.01);

        // Cross sections above and below K-edge should differ significantly
        let xs_above = fe.get_abs_coefficients(7.2);
        let xs_below = fe.get_abs_coefficients(7.0);
        let pe_above = xs_above[&CrossSection::Photoelectric];
        let pe_below = xs_below[&CrossSection::Photoelectric];

        // Above K-edge should have much higher photoelectric XS
        assert!(
            pe_above > pe_below * 2.0,
            "PE above K-edge ({}) should be >> below ({})",
            pe_above,
            pe_below
        );
    }

    #[test]
    fn test_sulfur_properties() {
        let db = ElementDatabase::instance();
        let s = db.get_by_number(16).unwrap();
        assert_eq!(s.name(), "S");
        assert!((s.atomic_weight() - 32.06).abs() < 0.1);

        // Sulfur K-edge ~2.472 keV
        let k_edge = s.k_edge().unwrap();
        assert!((k_edge - 2.472).abs() < 0.01);
    }

    #[test]
    fn test_elsepa_coefficients() {
        let db = ElementDatabase::instance();
        let c = db.get_by_number(6).unwrap();
        let elsepa = c.elsepa_coefficients();
        // Carbon should have ELSEPA data at 50-300 keV
        assert!(elsepa[0] > 0.0, "ELSEPA at 50 keV should be positive");
    }
}
