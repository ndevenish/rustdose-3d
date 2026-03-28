use crate::element::{CrossSection, ElementDatabase};

/// Avogadro's number.
const AVOGADRO: f64 = 6.022_140_857e23;

/// Conversion factor: microns to centimeters.
const MICRONS_TO_CM: f64 = 1e-4;

/// Container trait: models the sample container that attenuates the beam.
pub trait Container: std::fmt::Debug + Send + Sync {
    /// Calculate how much the container attenuates the beam.
    fn calculate_attenuation(&mut self, beam_energy: f64);
    /// Return the attenuation fraction (0 = transparent, 1 = fully opaque).
    fn attenuation_fraction(&self) -> f64;
    /// Return the container material name, if any.
    fn material_name(&self) -> Option<&str>;
    /// Print container information.
    fn info(&self) -> String;
}

/// Compute mass attenuation coefficient (cm²/g) for a compound given elemental
/// composition and beam energy, using the McMaster polynomial database.
///
/// For each element i with count n_i, atomic weight A_i, and total cross-section σ_i (barns):
///   µ/ρ = (N_A / (M × 1e24)) × Σ(n_i × σ_i)
/// where M = Σ(n_i × A_i) is the molecular weight.
fn compute_mass_attenuation_coeff(elements: &[(String, f64)], beam_energy: f64) -> f64 {
    let db = ElementDatabase::instance();

    let mut total_mol_weight = 0.0;
    let mut weighted_xs_sum = 0.0;

    for (symbol, count) in elements {
        if let Some(elem) = db.get(symbol) {
            let xs = elem.get_abs_coefficients(beam_energy);
            let total_xs = xs[&CrossSection::Total]; // barns/atom
            let atomic_weight = elem.atomic_weight(); // g/mol

            total_mol_weight += count * atomic_weight;
            weighted_xs_sum += count * total_xs;
        }
    }

    if total_mol_weight <= 0.0 {
        return 0.0;
    }

    // Convert: barns → cm² (×1e-24), then per gram via Avogadro/mol_weight
    AVOGADRO / (total_mol_weight * 1e24) * weighted_xs_sum
}

/// Transparent container: no attenuation at all.
#[derive(Debug, Clone)]
pub struct ContainerTransparent;

impl Container for ContainerTransparent {
    fn calculate_attenuation(&mut self, _beam_energy: f64) {}

    fn attenuation_fraction(&self) -> f64 {
        0.0
    }

    fn material_name(&self) -> Option<&str> {
        None
    }

    fn info(&self) -> String {
        "No container has been specified.".to_string()
    }
}

/// Container with a known material mixture (attenuates beam based on mass attenuation coefficient).
/// Note: proper attenuation requires NIST lookup (not yet implemented). Returns 0 attenuation.
#[derive(Debug)]
pub struct ContainerMixture {
    thickness_um: f64,
    density_g_per_ml: f64,
    material: String,
}

impl ContainerMixture {
    pub fn new(thickness_um: f64, density: f64, material: String) -> Self {
        ContainerMixture {
            thickness_um,
            density_g_per_ml: density,
            material,
        }
    }
}

impl Container for ContainerMixture {
    fn calculate_attenuation(&mut self, _beam_energy: f64) {
        println!(
            "ContainerMixture: NIST attenuation lookup not yet implemented, using 0 attenuation"
        );
    }

    fn attenuation_fraction(&self) -> f64 {
        0.0
    }

    fn material_name(&self) -> Option<&str> {
        Some(&self.material)
    }

    fn info(&self) -> String {
        format!(
            "Container: {} (thickness: {:.1} µm, density: {:.2} g/ml) - attenuation not yet calculated",
            self.material, self.thickness_um, self.density_g_per_ml
        )
    }
}

/// Container with elemental composition (attenuates beam via McMaster coefficients).
#[derive(Debug)]
pub struct ContainerElemental {
    thickness_um: f64,
    density_g_per_ml: f64,
    elements: Vec<(String, f64)>, // (symbol, count)
    mass_attenuation_coeff: f64,  // cm²/g
    attenuation_fraction: f64,
}

impl ContainerElemental {
    pub fn new(thickness_um: f64, density: f64, elements: Vec<(String, f64)>) -> Self {
        ContainerElemental {
            thickness_um,
            density_g_per_ml: density,
            elements,
            mass_attenuation_coeff: 0.0,
            attenuation_fraction: 0.0,
        }
    }
}

impl Container for ContainerElemental {
    fn calculate_attenuation(&mut self, beam_energy: f64) {
        if self.thickness_um <= 0.0 || self.elements.is_empty() {
            return;
        }

        self.mass_attenuation_coeff = compute_mass_attenuation_coeff(&self.elements, beam_energy);

        // mass_thickness = density (g/cm³) × thickness (cm)
        let thickness_cm = self.thickness_um * MICRONS_TO_CM;
        let mass_thickness = self.density_g_per_ml * thickness_cm;

        // Beer-Lambert law
        self.attenuation_fraction = 1.0 - (-self.mass_attenuation_coeff * mass_thickness).exp();
    }

    fn attenuation_fraction(&self) -> f64 {
        self.attenuation_fraction
    }

    fn material_name(&self) -> Option<&str> {
        None
    }

    fn info(&self) -> String {
        // Format as chemical formula: "SiO2" style
        let formula: String = self
            .elements
            .iter()
            .map(|(s, c)| {
                if (*c - 1.0).abs() < 1e-9 {
                    s.clone()
                } else if (*c - c.round()).abs() < 1e-9 {
                    format!("{}{}", s, *c as i64)
                } else {
                    format!("{}{}", s, c)
                }
            })
            .collect();
        if self.mass_attenuation_coeff > 0.0 {
            format!(
                "The mass attenuation coefficient of the {} container is {:.2} centimetres^2 per gram.\n\
                 The attenuation fraction of the beam due to the sample container of thickness {:.2} microns is: {:.2}.",
                formula,
                self.mass_attenuation_coeff,
                self.thickness_um,
                self.attenuation_fraction,
            )
        } else {
            format!(
                "Container: elemental [{}] (thickness: {:.1} µm, density: {:.2} g/ml)",
                formula, self.thickness_um, self.density_g_per_ml
            )
        }
    }
}

/// Create a container from CrystalConfig.
pub fn create_container(config: &crate::parser::config::CrystalConfig) -> Box<dyn Container> {
    use crate::parser::config::ContainerMaterialType;
    match config.container_material {
        Some(ContainerMaterialType::Mixture) => {
            let thickness = config.container_thickness.unwrap_or(0.0);
            let density = config.container_density.unwrap_or(1.0);
            let material = config.container_mixture.clone().unwrap_or_default();
            Box::new(ContainerMixture::new(thickness, density, material))
        }
        Some(ContainerMaterialType::Elemental) => {
            let thickness = config.container_thickness.unwrap_or(0.0);
            let density = config.container_density.unwrap_or(1.0);
            let elements = config
                .container_elements
                .iter()
                .map(|ec| (ec.symbol.clone(), ec.count))
                .collect();
            Box::new(ContainerElemental::new(thickness, density, elements))
        }
        _ => Box::new(ContainerTransparent),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transparent_container_zero_attenuation() {
        let mut c = ContainerTransparent;
        c.calculate_attenuation(8.05);
        assert!(
            c.attenuation_fraction() < 1e-6,
            "Transparent container should have zero attenuation"
        );
    }

    #[test]
    fn zero_thickness_containers_no_attenuation() {
        // ContainerMixture with zero thickness
        let mut cm = ContainerMixture::new(0.0, 2.23, "pyrex".to_string());
        cm.calculate_attenuation(8.05);
        assert!(
            cm.attenuation_fraction() < 1e-6,
            "Zero-thickness mixture container should not attenuate"
        );

        // ContainerElemental with zero thickness
        let mut ce = ContainerElemental::new(
            0.0,
            2.65,
            vec![("S".to_string(), 1.0), ("O".to_string(), 2.0)],
        );
        ce.calculate_attenuation(8.05);
        assert!(
            ce.attenuation_fraction() < 1e-6,
            "Zero-thickness elemental container should not attenuate"
        );
    }
}
