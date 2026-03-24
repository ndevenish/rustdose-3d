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
