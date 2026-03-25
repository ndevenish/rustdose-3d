use std::collections::{BTreeMap, HashMap};
use std::sync::OnceLock;

/// An element with electron elastic scattering cross-section data.
#[derive(Debug, Clone)]
pub struct ElementEM {
    name: String,
    atomic_number: i32,
    atomic_weight: f64,
    /// Energy (keV) → elastic cross-section. BTreeMap for ordered keys.
    data: BTreeMap<OrderedF64, f64>,
}

/// Wrapper for f64 that implements Ord (for BTreeMap keys).
/// Only used for positive energies, so NaN/negative not a concern.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderedF64(f64);

impl Eq for OrderedF64 {}

impl PartialOrd for OrderedF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl ElementEM {
    /// Get elastic scattering cross-section by linear interpolation.
    /// Valid range: 0.05–300 keV.
    pub fn elastic_coefficient(&self, electron_energy: f64) -> f64 {
        if !(0.05..=300.0).contains(&electron_energy) {
            return 0.0;
        }

        let key = OrderedF64(electron_energy);

        // Find floor and ceiling keys
        let before = self.data.range(..=key).next_back();
        let after = self.data.range(key..).next();

        let (before_key, before_val) = match before {
            Some((&k, &v)) => (k.0, v),
            None => match after {
                Some((&k, &v)) => (k.0, v),
                None => return 0.0,
            },
        };

        let (after_key, after_val) = match after {
            Some((&k, &v)) => (k.0, v),
            None => (before_key, before_val),
        };

        if (after_key - before_key).abs() < 1e-15 {
            return before_val;
        }

        let fraction = (electron_energy - before_key) / (after_key - before_key);
        before_val + (after_val - before_val) * fraction
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn atomic_number(&self) -> i32 {
        self.atomic_number
    }
    pub fn atomic_weight(&self) -> f64 {
        self.atomic_weight
    }
}

/// Singleton electron element database loaded from fullelsepa.csv.
pub struct ElementDatabaseEM {
    by_name: HashMap<String, ElementEM>,
    by_number: HashMap<i32, ElementEM>,
}

static INSTANCE: OnceLock<ElementDatabaseEM> = OnceLock::new();

const ELSEPA_DATA: &str = include_str!("../constants_data/fullelsepa.csv");

impl ElementDatabaseEM {
    fn load() -> Self {
        let mut by_name = HashMap::new();
        let mut by_number = HashMap::new();

        for line in ELSEPA_DATA.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let components: Vec<&str> = line.split(',').collect();
            if components.len() < 4 {
                continue;
            }

            let atomic_number: i32 = match components[0].trim().parse() {
                Ok(n) => n,
                Err(_) => continue,
            };
            let name = components[1].trim().to_string();
            let atomic_weight: f64 = match components[2].trim().parse() {
                Ok(w) => w,
                Err(_) => continue,
            };

            let mut data = BTreeMap::new();

            // Columns 3..=321: energy-dependent elastic cross sections
            // Columns 3-22 (i=3..22): energy = (i-2) * 0.05 keV
            // Columns 23-321 (i>=23): energy = (i-21) * 1.0 keV
            let max_col = components.len().min(322);
            for (i, comp) in components.iter().enumerate().take(max_col).skip(3) {
                let (minus_factor, multiply_factor) = if i >= 23 { (21, 1.0) } else { (2, 0.05) };
                let energy = (i as f64 - minus_factor as f64) * multiply_factor;
                let val = comp.trim();
                if !val.is_empty() {
                    if let Ok(xs) = val.parse::<f64>() {
                        data.insert(OrderedF64(energy), xs);
                    }
                }
            }

            let elem = ElementEM {
                name: name.clone(),
                atomic_number,
                atomic_weight,
                data,
            };
            by_name.insert(name.trim().to_lowercase(), elem.clone());
            by_number.insert(atomic_number, elem);
        }

        ElementDatabaseEM { by_name, by_number }
    }

    /// Get the singleton instance.
    pub fn instance() -> &'static ElementDatabaseEM {
        INSTANCE.get_or_init(Self::load)
    }

    pub fn get_by_number(&self, z: i32) -> Option<&ElementEM> {
        self.by_number.get(&z)
    }

    pub fn get(&self, name: &str) -> Option<&ElementEM> {
        self.by_name.get(&name.to_lowercase())
    }

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
    fn test_em_database_loads() {
        let db = ElementDatabaseEM::instance();
        assert!(
            db.len() >= 80,
            "Should have at least 80 elements, got {}",
            db.len()
        );
    }

    #[test]
    fn test_em_lookup() {
        let db = ElementDatabaseEM::instance();
        let h = db.get_by_number(1).unwrap();
        assert_eq!(h.atomic_number(), 1);
        assert!((h.atomic_weight() - 1.01).abs() < 0.01);

        let c = db.get("C").unwrap();
        assert_eq!(c.atomic_number(), 6);
    }

    #[test]
    fn test_em_interpolation() {
        let db = ElementDatabaseEM::instance();
        let c = db.get_by_number(6).unwrap();

        // At exact energy point (0.05 keV = first point)
        let xs_exact = c.elastic_coefficient(0.05);
        assert!(
            xs_exact > 0.0,
            "Cross section at 0.05 keV should be positive"
        );

        // Between energy points (0.075 keV = between 0.05 and 0.1)
        let xs_interp = c.elastic_coefficient(0.075);
        let xs_low = c.elastic_coefficient(0.05);
        let xs_high = c.elastic_coefficient(0.1);
        // Interpolated value should be between the two
        let (lo, hi) = if xs_low < xs_high {
            (xs_low, xs_high)
        } else {
            (xs_high, xs_low)
        };
        assert!(
            xs_interp >= lo && xs_interp <= hi,
            "Interpolated {} should be between {} and {}",
            xs_interp,
            lo,
            hi
        );

        // Out of range
        assert_eq!(c.elastic_coefficient(0.01), 0.0);
        assert_eq!(c.elastic_coefficient(400.0), 0.0);
    }

    #[test]
    fn test_em_hydrogen_values() {
        let db = ElementDatabaseEM::instance();
        let h = db.get_by_number(1).unwrap();

        // First data point from fullelsepa.csv line 1:
        // Column 3 = 7.350749E-03 at energy (3-2)*0.05 = 0.05 keV
        let xs = h.elastic_coefficient(0.05);
        assert!(
            (xs - 7.350749e-03).abs() < 1e-8,
            "H at 0.05 keV: expected ~7.35e-3, got {}",
            xs
        );

        // Last low-energy point: column 22 at energy (22-2)*0.05 = 1.0 keV
        // Also first high-energy point: column 23 at energy (23-21)*1.0 = 2.0 keV
        let xs_2kev = h.elastic_coefficient(2.0);
        assert!(xs_2kev > 0.0);
    }

    #[test]
    fn test_em_iron() {
        let db = ElementDatabaseEM::instance();
        let fe = db.get_by_number(26).unwrap();
        assert_eq!(fe.name().trim(), "FE");

        // Iron should have larger cross sections than hydrogen
        let fe_xs = fe.elastic_coefficient(100.0);
        let h = db.get_by_number(1).unwrap();
        let h_xs = h.elastic_coefficient(100.0);
        assert!(
            fe_xs > h_xs,
            "Fe XS ({}) should exceed H XS ({})",
            fe_xs,
            h_xs
        );
    }
}
