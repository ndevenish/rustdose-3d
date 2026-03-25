use super::compute::CoefCalcCompute;
use super::CoefCalc;

/// CoefCalcFromCIF: compute coefficients from a CIF file's chemical formula.
///
/// Reads `_chemical_formula_sum` and `_cell_volume` from the CIF file.
#[derive(Debug)]
pub struct CoefCalcFromCIF {
    pub compute: CoefCalcCompute,
}

impl CoefCalcFromCIF {
    /// Parse a CIF file at `path` and build the coefficient calculator.
    pub fn from_file(path: &str) -> Result<Self, String> {
        let content =
            std::fs::read_to_string(path).map_err(|e| format!("Cannot read CIF file: {e}"))?;

        let mut compute = CoefCalcCompute::new();
        let mut found_formula = false;

        for line in content.lines() {
            let line = line.trim();

            if line.starts_with("_chemical_formula_sum") {
                parse_chemical_formula(line, &mut compute)?;
                found_formula = true;
            }

            if line.starts_with("_cell_volume") {
                if let Some(vol) = parse_cell_volume(line) {
                    compute.cell_volume = vol;
                }
            }
        }

        if !found_formula {
            return Err("CIF file must contain _chemical_formula_sum".to_string());
        }

        compute.calculate_density();
        Ok(CoefCalcFromCIF { compute })
    }
}

/// Parse `_chemical_formula_sum 'C3 H7 N O2'` into element occurrences.
fn parse_chemical_formula(line: &str, compute: &mut CoefCalcCompute) -> Result<(), String> {
    // Extract the quoted content
    let first = line.find('\'').ok_or("No quote in formula line")?;
    let last = line.rfind('\'').ok_or("No closing quote in formula line")?;
    if first == last {
        return Err("Malformed formula line".to_string());
    }
    let formula = &line[first + 1..last];

    // Parse space-separated tokens like "C3 H7 N O2"
    for token in formula.split_whitespace() {
        if token.is_empty() {
            continue;
        }
        // Find where digits start (if any)
        let split_pos = token.find(|c: char| c.is_ascii_digit() || c == '.');
        let (sym, count_str) = match split_pos {
            Some(p) => (&token[..p], &token[p..]),
            None => (token, "1"),
        };
        let count: f64 = count_str.parse().unwrap_or(1.0);
        let sym_upper = sym.to_uppercase();
        compute.increment_macro(&sym_upper, count);
    }
    Ok(())
}

/// Parse `_cell_volume 123.4(5)` → 123.4.
fn parse_cell_volume(line: &str) -> Option<f64> {
    // Find first digit in the line
    let start = line.find(|c: char| c.is_ascii_digit())?;
    let rest = &line[start..];
    // Stop at '(' for uncertainty or whitespace
    let end = rest
        .find(|c: char| c == '(' || c.is_whitespace())
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}

impl CoefCalc for CoefCalcFromCIF {
    fn update_coefficients(&mut self, photon_energy: f64) {
        let (photo, coherent, compton, total) =
            self.compute.calculate_coefficients_all(photon_energy);
        self.compute.abs_coeff_photo = photo;
        self.compute.elas_coeff = coherent;
        self.compute.abs_coeff_comp = compton;
        self.compute.att_coeff = total;
        let (_p, coh_m, _c, _t) = self.compute.calculate_coefficients_macro(photon_energy);
        self.compute.elas_coeff_macro = coh_m;
    }

    fn absorption_coefficient(&self) -> f64 {
        self.compute.abs_coeff_photo
    }
    fn attenuation_coefficient(&self) -> f64 {
        self.compute.att_coeff
    }
    fn elastic_coefficient(&self) -> f64 {
        self.compute.elas_coeff
    }
    fn inelastic_coefficient(&self) -> f64 {
        self.compute.abs_coeff_comp
    }
    fn density(&self) -> f64 {
        self.compute.crystal_density
    }
    fn fluorescent_escape_factors(&self, beam_energy: f64) -> Vec<Vec<f64>> {
        self.compute.calc_fluorescent_escape_factors(beam_energy)
    }
    fn solvent_fraction(&self) -> f64 {
        self.compute.sol_fraction
    }
}
