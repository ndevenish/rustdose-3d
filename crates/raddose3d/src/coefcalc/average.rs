/// CoefCalcAverage: hardcoded average crystal composition (Holton 2010, 50% solvent).
///
/// Coefficients at 12.4 keV for an average protein crystal (Berman 2002 survey).
/// Absorption: 0.237/mm, Attenuation: 0.281/mm, Elastic: 0.01799/mm, Density: 1.2 g/ml.
#[derive(Debug)]
pub struct CoefCalcAverage;

/// Absorption coefficient in µm⁻¹ (Holton 2010).
const ABSORPTION_COEFFICIENT: f64 = 0.000237;
/// Attenuation coefficient in µm⁻¹.
const ATTENUATION_COEFFICIENT: f64 = 0.000281;
/// Elastic (coherent) scattering coefficient in µm⁻¹.
const ELASTIC_COEFFICIENT: f64 = 0.00001799;
/// Crystal density in g/ml (Fischer & Polikarpov 2004 + 50% solvent).
const DENSITY: f64 = 1.2;

impl super::CoefCalc for CoefCalcAverage {
    fn update_coefficients(&mut self, _photon_energy: f64) {
        // No-op: hardcoded values are energy-independent
    }

    fn absorption_coefficient(&self) -> f64 {
        ABSORPTION_COEFFICIENT
    }

    fn attenuation_coefficient(&self) -> f64 {
        ATTENUATION_COEFFICIENT
    }

    fn elastic_coefficient(&self) -> f64 {
        ELASTIC_COEFFICIENT
    }

    fn inelastic_coefficient(&self) -> f64 {
        // Compton already included in absorption for Average mode
        0.0
    }

    fn density(&self) -> f64 {
        DENSITY
    }

    fn fluorescent_escape_factors(&self, _beam_energy: f64) -> Vec<Vec<f64>> {
        println!(
            "WARNING: No X-ray fluorescent escape correction is implemented \
             for the Average/Dummy crystal composition. No correction applied."
        );
        // Return single row of zeros (27 factors)
        vec![vec![0.0; 27]]
    }
}
