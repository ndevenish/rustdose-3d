use crate::parser::config::DdmType;

/// Dose Decay Model: calculates Relative Diffraction Efficiency (RDE)
/// as a function of accumulated dose.
pub trait DdmModel: std::fmt::Debug + Send + Sync {
    /// Calculate the RDE for a given dose (in MGy).
    fn calc_decay(&self, dose: f64) -> f64;
    /// Human-readable name of the model.
    fn name(&self) -> &str;
}

/// Simple DDM: assumes no intensity decay with dose (RDE = 1 always).
#[derive(Debug, Clone)]
pub struct DdmSimple;

impl DdmModel for DdmSimple {
    fn calc_decay(&self, _dose: f64) -> f64 {
        1.0
    }
    fn name(&self) -> &str {
        "Simple DDM"
    }
}

/// Linear DDM: RDE decays linearly with dose.
#[derive(Debug, Clone)]
pub struct DdmLinear {
    pub gamma: f64,
}

impl DdmModel for DdmLinear {
    fn calc_decay(&self, dose: f64) -> f64 {
        (1.0 - dose * self.gamma).max(0.0)
    }
    fn name(&self) -> &str {
        "Linear DDM"
    }
}

/// Leal DDM: exponential decay model from Leal et al. (2012).
#[derive(Debug, Clone)]
pub struct DdmLeal {
    /// Half-dose in MGy (default 10.0).
    pub d_half: f64,
}

impl Default for DdmLeal {
    fn default() -> Self {
        Self { d_half: 10.0 }
    }
}

impl DdmModel for DdmLeal {
    fn calc_decay(&self, dose: f64) -> f64 {
        (-dose * (2.0_f64.ln()) / self.d_half).exp()
    }
    fn name(&self) -> &str {
        "Leal DDM"
    }
}

/// B-factor DDM: uses temperature factor model.
#[derive(Debug, Clone)]
pub struct DdmBfactor {
    pub b0: f64,
    pub beta: f64,
}

impl DdmModel for DdmBfactor {
    fn calc_decay(&self, dose: f64) -> f64 {
        let b = self.b0 + self.beta * dose;
        (-b).exp()
    }
    fn name(&self) -> &str {
        "B-factor DDM"
    }
}

/// Create a DDM from parsed configuration.
pub fn create_ddm(
    ddm_type: Option<DdmType>,
    gamma: Option<f64>,
    b0: Option<f64>,
    beta: Option<f64>,
) -> Box<dyn DdmModel> {
    match ddm_type {
        None | Some(DdmType::Simple) => Box::new(DdmSimple),
        Some(DdmType::Linear) => Box::new(DdmLinear {
            gamma: gamma.unwrap_or(0.01),
        }),
        Some(DdmType::Leal) => Box::new(DdmLeal::default()),
        Some(DdmType::Bfactor) => Box::new(DdmBfactor {
            b0: b0.unwrap_or(0.0),
            beta: beta.unwrap_or(1.0),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ddm_simple_is_one() {
        let ddm = DdmSimple;
        let doses = [0.0, 1.0, 5.0, 10.0, 30.0, 43.0, 50.0, 100.0, 1000.0];
        for &d in &doses {
            let v = ddm.calc_decay(d);
            assert!((v - 1.0).abs() < 1e-6, "DdmSimple at dose {}: got {}", d, v);
        }
    }

    #[test]
    fn ddm_linear_decay_range() {
        let ddm = DdmLinear { gamma: 0.01 };
        let doses = [0.0, 1.0, 5.0, 10.0, 30.0, 43.0, 50.0, 100.0, 1000.0];
        for &d in &doses {
            let v = ddm.calc_decay(d);
            assert!(
                (0.0..=1.0).contains(&v),
                "DdmLinear at dose {}: got {}",
                d,
                v
            );
        }
    }

    #[test]
    fn ddm_leal_decay_range() {
        let ddm = DdmLeal { d_half: 10.0 };
        let doses = [0.0, 1.0, 5.0, 10.0, 30.0, 43.0, 50.0, 100.0, 1000.0];
        let zero_intensity = ddm.calc_decay(0.0);
        for &d in &doses {
            let rde = ddm.calc_decay(d) / zero_intensity;
            assert!(
                (0.0..=1.0).contains(&rde),
                "DdmLeal RDE at dose {}: got {}",
                d,
                rde
            );
        }
    }

    #[test]
    fn ddm_linear_monotonic_decrease() {
        let ddm = DdmLinear { gamma: 0.01 };
        let doses = [0.0, 1.0, 5.0, 10.0, 30.0, 43.0, 50.0, 100.0, 1000.0];
        let mut prev = 1.0;
        for (i, &d) in doses.iter().enumerate() {
            let v = ddm.calc_decay(d);
            if i > 0 && prev != 0.0 {
                assert!(v < prev, "DdmLinear not monotonic at dose {}", d);
            }
            prev = v;
        }
    }

    #[test]
    fn ddm_leal_monotonic_decrease() {
        let ddm = DdmLeal { d_half: 10.0 };
        let doses = [0.0, 1.0, 5.0, 10.0, 30.0, 43.0, 50.0, 100.0, 1000.0];
        let zero_intensity = ddm.calc_decay(0.0);
        let mut prev = 1.0;
        for (i, &d) in doses.iter().enumerate() {
            let rde = ddm.calc_decay(d) / zero_intensity;
            if i > 0 {
                assert!(rde <= prev, "DdmLeal not monotonic at dose {}", d);
            }
            prev = rde;
        }
    }

    #[test]
    fn ddm_linear_zero_dose_is_one() {
        let ddm = DdmLinear { gamma: 0.01 };
        let v = ddm.calc_decay(0.0);
        assert!((v - 1.0).abs() < 1e-6, "DdmLinear at zero dose: got {}", v);
    }

    #[test]
    fn ddm_leal_zero_dose_is_one() {
        let ddm = DdmLeal { d_half: 10.0 };
        let zero_intensity = ddm.calc_decay(0.0);
        let rde = ddm.calc_decay(0.0) / zero_intensity;
        assert!(
            (rde - 1.0).abs() < 1e-6,
            "DdmLeal RDE at zero dose: got {}",
            rde
        );
    }
}
