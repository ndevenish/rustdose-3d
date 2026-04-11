use crate::element::ElementDatabase;

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

/// Fetch mass attenuation coefficient (µ/ρ in cm²/g) for a single element
/// from the NIST X-Ray Mass Attenuation Coefficients database.
///
/// Downloads the table from NIST, parses the energy vs µ/ρ data, and
/// linearly interpolates to the requested beam energy.
///
/// Returns None if the download or parsing fails.
/// Not available on wasm32 (no blocking HTTP).
#[cfg(not(target_arch = "wasm32"))]
fn fetch_nist_mass_attenuation(atomic_number: i32, beam_energy_kev: f64) -> Option<f64> {
    let url = format!(
        "https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z{:02}.html",
        atomic_number
    );

    let body = ureq::get(&url)
        .call()
        .ok()?
        .body_mut()
        .read_to_string()
        .ok()?;

    // Parse data between <PRE> and </PRE> tags, matching lines with scientific notation
    let mut in_pre = false;
    let mut prev_energy_kev = 0.0_f64;
    let mut prev_mu_rho = 0.0_f64;

    for line in body.lines() {
        if line.contains("<PRE>") {
            in_pre = true;
            continue;
        }
        if line.contains("</PRE>") {
            in_pre = false;
            continue;
        }

        if !in_pre {
            continue;
        }

        // Look for lines containing scientific notation (digits, E, +/-)
        if !line.contains('E') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        // Find the first column containing 'E' (scientific notation)
        let energy_col = parts.iter().position(|p| p.contains('E'));
        let energy_col = match energy_col {
            Some(i) => i,
            None => continue,
        };

        if parts.len() < energy_col + 2 {
            continue;
        }

        let energy_mev: f64 = match parts[energy_col].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let mu_rho: f64 = match parts[energy_col + 1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let energy_kev = energy_mev * MEV_TO_KEV;

        if energy_kev > beam_energy_kev {
            // Linearly interpolate between previous and current
            if prev_energy_kev <= 0.0 {
                return Some(mu_rho);
            }
            let frac = (beam_energy_kev - prev_energy_kev) / (energy_kev - prev_energy_kev);
            return Some(prev_mu_rho + (mu_rho - prev_mu_rho) * frac);
        }

        prev_energy_kev = energy_kev;
        prev_mu_rho = mu_rho;
    }

    // Beam energy beyond last table entry — use last value
    if prev_mu_rho > 0.0 {
        Some(prev_mu_rho)
    } else {
        None
    }
}

/// Compute the weighted-average mass attenuation coefficient (cm²/g) for a
/// compound, using NIST data for each element. This matches the Java
/// RADDOSE-3D approach exactly (ContainerElemental.extractMassAttenuationCoef).
#[cfg(not(target_arch = "wasm32"))]
fn compute_nist_mass_attenuation(elements: &[(String, f64)], beam_energy_kev: f64) -> Option<f64> {
    let db = ElementDatabase::instance();

    let mut element_mu_rhos = Vec::with_capacity(elements.len());
    let mut element_weights = Vec::with_capacity(elements.len());
    let mut total_weight = 0.0_f64;

    for (symbol, count) in elements {
        let elem = db.get(symbol)?;
        let z = elem.atomic_number();
        let atomic_weight = elem.atomic_weight();

        let mu_rho = fetch_nist_mass_attenuation(z, beam_energy_kev)?;

        element_mu_rhos.push(mu_rho);
        element_weights.push(count * atomic_weight);
        total_weight += count * atomic_weight;
    }

    if total_weight <= 0.0 {
        return None;
    }

    // Weighted average: µ/ρ = Σ(w_i × µ_i/ρ_i)
    // where w_i = (count_i × A_i) / total_weight
    let mut mass_atten = 0.0;
    for i in 0..elements.len() {
        let relative_weight = element_weights[i] / total_weight;
        mass_atten += relative_weight * element_mu_rhos[i];
    }

    Some(mass_atten)
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

/// Fetch mass attenuation coefficient (µ/ρ in cm²/g) for a compound mixture
/// from the NIST X-Ray Mass Attenuation Coefficients ComTab database.
///
/// Downloads the HTML table, parses TD cells with scientific notation
/// in groups of 3 (energy MeV, µ/ρ, µ_en/ρ), and linearly interpolates
/// to the requested beam energy. Matches Java's ContainerMixture algorithm.
#[cfg(not(target_arch = "wasm32"))]
fn fetch_nist_mixture_attenuation(material: &str, beam_energy_kev: f64) -> Option<f64> {
    let url = format!(
        "https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/{}.html",
        material
    );

    let body = ureq::get(&url)
        .call()
        .ok()?
        .body_mut()
        .read_to_string()
        .ok()?;

    // Extract all TD cell values containing scientific notation (e.g. 1.000E-03)
    // matching Java's Pattern: <TD>([0-9]\.[0-9]+E[+-][0-9]+)</TD>
    let mut scientific_values: Vec<f64> = Vec::new();
    let mut search = body.as_str();
    while let Some(td_start) = search.find("<TD>") {
        let after = &search[td_start + 4..];
        if let Some(td_end) = after.find("</TD>") {
            let cell = &after[..td_end];
            // Check if cell matches scientific notation: digit.digits E +/- digits
            if looks_like_scientific(cell) {
                if let Ok(v) = cell.trim().parse::<f64>() {
                    scientific_values.push(v);
                }
            }
            search = &after[td_end + 5..];
        } else {
            break;
        }
    }

    // Process in groups of 3: (energy_mev, mu_rho, mu_en_rho)
    let mut prev_energy_kev = 0.0_f64;
    let mut prev_mu_rho = 0.0_f64;
    let mut i = 0;
    while i + 2 < scientific_values.len() {
        let energy_mev = scientific_values[i];
        let mu_rho = scientific_values[i + 1];
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
        i += 3;
    }

    // Beam energy beyond last table entry — use last value
    if prev_mu_rho > 0.0 {
        Some(prev_mu_rho)
    } else {
        None
    }
}

/// Return true if `s` looks like NIST scientific notation (e.g. "1.000E-03").
fn looks_like_scientific(s: &str) -> bool {
    let s = s.trim();
    if let Some(e_pos) = s.find('E') {
        let before = &s[..e_pos];
        let after = &s[e_pos + 1..];
        before.contains('.')
            && before.chars().next().is_some_and(|c| c.is_ascii_digit())
            && !after.is_empty()
    } else {
        false
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

        #[cfg(not(target_arch = "wasm32"))]
        match fetch_nist_mixture_attenuation(&self.material, beam_energy) {
            Some(mu_rho) => {
                self.mass_attenuation_coeff = mu_rho;
                let thickness_cm = self.thickness_um * MICRONS_TO_CM;
                let mass_thickness = self.density_g_per_ml * thickness_cm;
                self.attenuation_fraction =
                    1.0 - (-self.mass_attenuation_coeff * mass_thickness).exp();
            }
            None => {
                eprintln!(
                    "Warning: failed to fetch NIST ComTab data for container '{}', using 0 attenuation",
                    self.material
                );
            }
        }

        #[cfg(target_arch = "wasm32")]
        let _ = beam_energy;
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

        #[cfg(not(target_arch = "wasm32"))]
        match compute_nist_mass_attenuation(&self.elements, beam_energy) {
            Some(mu_rho) => {
                self.mass_attenuation_coeff = mu_rho;
            }
            None => {
                eprintln!(
                    "Warning: failed to fetch NIST data for container {}, using 0 attenuation",
                    self.formula()
                );
                return;
            }
        }

        #[cfg(not(target_arch = "wasm32"))]
        {
            // mass_thickness = density (g/cm³) × thickness (cm)
            let thickness_cm = self.thickness_um * MICRONS_TO_CM;
            let mass_thickness = self.density_g_per_ml * thickness_cm;

            // Beer-Lambert law
            self.attenuation_fraction = 1.0 - (-self.mass_attenuation_coeff * mass_thickness).exp();
        }

        #[cfg(target_arch = "wasm32")]
        let _ = beam_energy;
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
}
