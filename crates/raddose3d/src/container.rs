use crate::element::ElementDatabase;
use crate::nist_tables;

/// Conversion factor: microns to centimeters.
const MICRONS_TO_CM: f64 = 1e-4;

/// Conversion factor: MeV to keV.
const MEV_TO_KEV: f64 = 1e3;

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

/// Linearly interpolate µ/ρ from a sorted (energy_MeV, mu_rho) table
/// for the given beam energy in keV.  Returns None if the table is empty.
fn interpolate_mu_rho(table: &[(f64, f64)], beam_energy_kev: f64) -> Option<f64> {
    if table.is_empty() {
        return None;
    }
    let mut prev_energy_kev = 0.0_f64;
    let mut prev_mu_rho = 0.0_f64;
    for &(energy_mev, mu_rho) in table {
        let energy_kev = energy_mev * MEV_TO_KEV;
        if energy_kev > beam_energy_kev {
            if prev_energy_kev <= 0.0 {
                return Some(mu_rho);
            }
            let frac = (beam_energy_kev - prev_energy_kev) / (energy_kev - prev_energy_kev);
            return Some(prev_mu_rho + (mu_rho - prev_mu_rho) * frac);
        }
        prev_energy_kev = energy_kev;
        prev_mu_rho = mu_rho;
    }
    // Beam energy beyond last table entry — use last value.
    if prev_mu_rho > 0.0 {
        Some(prev_mu_rho)
    } else {
        None
    }
}

/// Look up mass attenuation coefficient for a NIST ComTab mixture by name.
fn mixture_mu_rho(material: &str, beam_energy_kev: f64) -> Option<f64> {
    let lower = material.to_lowercase();
    let idx = nist_tables::MIXTURE_NAMES
        .iter()
        .position(|&name| name == lower.as_str())?;
    let table = nist_tables::MIXTURE_MU_RHO[idx];
    interpolate_mu_rho(table, beam_energy_kev)
}

/// Look up mass attenuation coefficient for an element by atomic number Z.
fn elem_mu_rho_by_z(z: i32, beam_energy_kev: f64) -> Option<f64> {
    let idx = (z as usize).checked_sub(1)?;
    let table = nist_tables::ELEM_MU_RHO.get(idx)?;
    interpolate_mu_rho(table, beam_energy_kev)
}

/// Compute the weighted-average µ/ρ (cm²/g) for a compound from its elements.
/// Weights are by mass fraction: w_i = count_i × A_i / Σ(count_j × A_j).
fn compute_elemental_mu_rho(elements: &[(String, f64)], beam_energy_kev: f64) -> Option<f64> {
    let db = ElementDatabase::instance();
    let mut total_weight = 0.0_f64;
    let mut weighted_sum = 0.0_f64;
    for (symbol, count) in elements {
        let elem = db.get(symbol)?;
        let z = elem.atomic_number();
        let atomic_weight = elem.atomic_weight();
        let mu_rho = elem_mu_rho_by_z(z, beam_energy_kev)?;
        let mass = count * atomic_weight;
        weighted_sum += mass * mu_rho;
        total_weight += mass;
    }
    if total_weight > 0.0 {
        Some(weighted_sum / total_weight)
    } else {
        None
    }
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

/// Container with a known material mixture — attenuates beam via NIST ComTab lookup.
#[derive(Debug)]
pub struct ContainerMixture {
    thickness_um: f64,
    density_g_per_ml: f64,
    material: String,
    mass_attenuation_coeff: f64, // cm²/g
    attenuation_fraction: f64,
}

impl ContainerMixture {
    pub fn new(thickness_um: f64, density: f64, material: String) -> Self {
        ContainerMixture {
            thickness_um,
            density_g_per_ml: density,
            material,
            mass_attenuation_coeff: 0.0,
            attenuation_fraction: 0.0,
        }
    }
}

impl Container for ContainerMixture {
    fn calculate_attenuation(&mut self, beam_energy: f64) {
        if self.thickness_um <= 0.0 || self.material.is_empty() {
            return;
        }
        match mixture_mu_rho(&self.material, beam_energy) {
            Some(mu_rho) => {
                self.mass_attenuation_coeff = mu_rho;
                let thickness_cm = self.thickness_um * MICRONS_TO_CM;
                let mass_thickness = self.density_g_per_ml * thickness_cm;
                self.attenuation_fraction =
                    1.0 - (-self.mass_attenuation_coeff * mass_thickness).exp();
            }
            None => {
                eprintln!(
                    "Warning: NIST ComTab data not found for container material '{}', using 0 attenuation",
                    self.material
                );
            }
        }
    }

    fn attenuation_fraction(&self) -> f64 {
        self.attenuation_fraction
    }

    fn material_name(&self) -> Option<&str> {
        Some(&self.material)
    }

    fn info(&self) -> String {
        if self.mass_attenuation_coeff > 0.0 {
            format!(
                "The mass attenuation coefficient of the {} container is {:.2} centimetres^2 per gram.\n\
                 The attenuation fraction of the beam due to the sample container of thickness {:.2} microns is: {:.2}.",
                self.material,
                self.mass_attenuation_coeff,
                self.thickness_um,
                self.attenuation_fraction,
            )
        } else {
            format!(
                "Container: {} (thickness: {:.1} µm, density: {:.2} g/ml)",
                self.material, self.thickness_um, self.density_g_per_ml
            )
        }
    }
}

/// Container with elemental composition (attenuates beam via NIST mass attenuation data).
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

    /// Build chemical formula string (e.g. "SiO2").
    fn formula(&self) -> String {
        self.elements
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
            .collect()
    }
}

impl Container for ContainerElemental {
    fn calculate_attenuation(&mut self, beam_energy: f64) {
        if self.thickness_um <= 0.0 || self.elements.is_empty() {
            return;
        }
        match compute_elemental_mu_rho(&self.elements, beam_energy) {
            Some(mu_rho) => {
                self.mass_attenuation_coeff = mu_rho;
                let thickness_cm = self.thickness_um * MICRONS_TO_CM;
                let mass_thickness = self.density_g_per_ml * thickness_cm;
                self.attenuation_fraction =
                    1.0 - (-self.mass_attenuation_coeff * mass_thickness).exp();
            }
            None => {
                eprintln!(
                    "Warning: failed to look up NIST attenuation data for container {}, using 0 attenuation",
                    self.formula()
                );
            }
        }
    }

    fn attenuation_fraction(&self) -> f64 {
        self.attenuation_fraction
    }

    fn material_name(&self) -> Option<&str> {
        None
    }

    fn info(&self) -> String {
        let formula = self.formula();
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

    #[test]
    fn water_mixture_attenuates_correctly() {
        // Water at 17 keV: NIST gives µ/ρ ≈ 1.355 cm²/g (interpolated near that range)
        // Beer-Lambert: 1 - exp(-µ/ρ × ρ × t) with ρ=1.0 g/cm³, t=100µm → small attenuation
        let mut cm = ContainerMixture::new(100.0, 1.0, "water".to_string());
        cm.calculate_attenuation(17.434);
        let frac = cm.attenuation_fraction();
        // Just verify it's a reasonable positive value
        assert!(
            frac > 0.0 && frac < 1.0,
            "Expected attenuation fraction in (0,1), got {frac}"
        );
    }

    #[test]
    fn elemental_container_attenuates_correctly() {
        // Silicon at 8.05 keV
        let mut ce = ContainerElemental::new(100.0, 2.33, vec![("Si".to_string(), 1.0)]);
        ce.calculate_attenuation(8.05);
        let frac = ce.attenuation_fraction();
        assert!(
            frac > 0.0 && frac < 1.0,
            "Expected attenuation fraction in (0,1), got {frac}"
        );
    }
}
