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
        eprintln!(
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

/// Container with elemental composition (attenuates beam). Stub implementation.
#[derive(Debug)]
pub struct ContainerElemental {
    thickness_um: f64,
    density_g_per_ml: f64,
    elements: Vec<(String, f64)>, // (symbol, count)
}

impl ContainerElemental {
    pub fn new(thickness_um: f64, density: f64, elements: Vec<(String, f64)>) -> Self {
        ContainerElemental {
            thickness_um,
            density_g_per_ml: density,
            elements,
        }
    }
}

impl Container for ContainerElemental {
    fn calculate_attenuation(&mut self, _beam_energy: f64) {
        eprintln!(
            "ContainerElemental: attenuation calculation not yet implemented, using 0 attenuation"
        );
    }

    fn attenuation_fraction(&self) -> f64 {
        0.0
    }

    fn material_name(&self) -> Option<&str> {
        None
    }

    fn info(&self) -> String {
        let elem_str: Vec<String> = self
            .elements
            .iter()
            .map(|(s, c)| format!("{} {}", s, c))
            .collect();
        format!(
            "Container: elemental [{}] (thickness: {:.1} µm, density: {:.2} g/ml) - attenuation not yet calculated",
            elem_str.join(", "),
            self.thickness_um,
            self.density_g_per_ml
        )
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
