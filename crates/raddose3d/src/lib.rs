pub mod beam;
pub mod coefcalc;
pub mod container;
pub mod crystal;
pub mod ddm;
pub mod element;
pub mod energy_distribution;
pub mod experiment;
pub mod output;
pub mod residue;
pub mod simulation;
pub mod wedge;
pub mod writer;

pub use raddose3d_parser as parser;

/// Physical constants used throughout RADDOSE-3D.
pub mod constants {
    /// Elementary charge in coulombs.
    pub const ELEMENTARY_CHARGE: f64 = 1.602_176_565e-19;
    /// Conversion factor from keV to Joules.
    pub const KEV_TO_JOULES: f64 = 1000.0 * ELEMENTARY_CHARGE;
    /// Atomic mass unit in grams.
    pub const ATOMIC_MASS_UNIT: f64 = 1.66e-24;
    /// Avogadro's number.
    pub const AVOGADRO_NUM: f64 = 6.022e23;
    /// Water concentration in mM.
    pub const WATER_CONCENTRATION: f64 = 51666.0;
    /// Angstroms cubed to mL conversion.
    pub const ANGSTROMS_TO_ML: f64 = 1e-24;
    /// Number of static exposure steps when crystal is exposed without rotation.
    pub const STATIC_EXPOSURE: usize = 100;
    /// Conversion factor from Gy to MGy.
    pub const GY_TO_MGY: f64 = 1e-6;
    /// Unit conversion to get voxel mass in kg (1e-15 = 1e-18/1e-3).
    pub const UNIT_CONVERSION: f64 = 1e-15;
}

// ── High-level public API ────────────────────────────────────────────────────

/// Parse a RADDOSE-3D text input string into a `Config`.
///
/// This is the standard `.txt` format used by the Java RADDOSE-3D tool.
/// For JSON input, use [`parse_input_json`].
///
/// Relative file paths (SeqFile, PDB, CIF, ModelFile, beam File) are left
/// as-is. If the input came from a file, prefer [`parse_input_file`] which
/// resolves them relative to the file's directory.
///
/// # Errors
/// Returns an error string if parsing fails.
pub fn parse_input(input: &str) -> Result<parser::Config, String> {
    raddose3d_parser::parse(input).map_err(|e| e.to_string())
}

/// Parse a RADDOSE-3D input file, resolving relative paths against
/// the file's parent directory.
///
/// # Errors
/// Returns an error string if the file cannot be read or parsing fails.
pub fn parse_input_file(path: &std::path::Path) -> Result<parser::Config, String> {
    let input = std::fs::read_to_string(path)
        .map_err(|e| format!("Cannot read input file {}: {}", path.display(), e))?;
    let mut config = parse_input(&input)?;
    if let Some(base) = path.parent() {
        config.resolve_paths(base);
    }
    Ok(config)
}

/// Parse a RADDOSE-3D config from JSON.
///
/// The JSON must match the structure of [`parser::Config`] (the same struct
/// that `parse_input` produces). This is useful for programmatic input
/// generation.
///
/// # Errors
/// Returns an error string if JSON deserialization fails.
pub fn parse_input_json(json: &str) -> Result<parser::Config, String> {
    serde_json::from_str(json).map_err(|e| e.to_string())
}

/// Summary results from a completed RADDOSE-3D simulation.
///
/// When a config contains multiple wedges, these values reflect the
/// last wedge's cumulative exposure summary.
#[derive(Debug, Clone)]
pub struct RunResults {
    /// Average diffraction-weighted dose (MGy).
    pub average_dwd: f64,
    /// Last diffraction-weighted dose (MGy).
    pub last_dwd: f64,
    /// Maximum voxel dose (MGy).
    pub max_dose: f64,
    /// Average dose across whole crystal volume (MGy).
    pub avg_dose_whole_crystal: f64,
    /// Average dose across the exposed region only (MGy).
    pub avg_dose_exposed_region: f64,
    /// Total elastic (diffracted) photon count for the wedge.
    pub elastic_yield: f64,
    /// Total absorbed energy for the wedge (Joules).
    pub absorbed_energy: f64,
    /// Fraction of crystal volume that received any dose (0.0–100.0).
    pub used_volume_fraction: f64,
}

/// Run a RADDOSE-3D simulation from a parsed config, returning summary results.
///
/// The returned [`RunResults`] carries key numeric results from the last wedge.
///
/// For full output (CSV files, voxel dose, etc.) use [`experiment::Experiment`]
/// directly and attach the desired [`output::Output`] observers.
///
/// # Errors
/// Returns an error string if crystal/beam construction fails.
pub fn run(config: &parser::Config) -> Result<RunResults, String> {
    use std::sync::{Arc, Mutex};

    let results: Arc<Mutex<Option<RunResults>>> = Arc::new(Mutex::new(None));
    let collector = ResultsCollector(Arc::clone(&results));

    let mut exp = experiment::Experiment::new();
    exp.add_observer(Box::new(collector));

    use parser::config::ConfigItem;
    for item in &config.items {
        match item {
            ConfigItem::Crystal(c) => {
                let crystal = crystal::create_crystal(c)?;
                exp.set_crystal(crystal);
            }
            ConfigItem::Beam(b) => {
                let beam = beam::create_beam(b)?;
                exp.set_beam(beam);
            }
            ConfigItem::Wedge(w) => {
                let wdg = wedge::Wedge::from_config(w);
                exp.expose_wedge(&wdg);
            }
        }
    }
    exp.close();

    let guard = results.lock().unwrap();
    guard
        .clone()
        .ok_or_else(|| "No wedge was exposed".to_string())
}

// ── Internal: results-collecting Output observer ─────────────────────────────

use std::sync::{Arc, Mutex};

#[derive(Debug)]
struct ResultsCollector(Arc<Mutex<Option<RunResults>>>);

impl output::Output for ResultsCollector {
    fn publish_crystal(&mut self, _crystal: &dyn crystal::Crystal) {}
    fn publish_beam(&mut self, _beam: &dyn beam::Beam) {}
    fn publish_wedge(
        &mut self,
        _wedge: &wedge::Wedge,
        summary: &output::ExposureSummary,
        _crystal: Option<&dyn crystal::Crystal>,
    ) {
        *self.0.lock().unwrap() = Some(RunResults {
            average_dwd: summary.avg_diffracted_dose(),
            last_dwd: summary.last_dwd(),
            max_dose: summary.max_dose(),
            avg_dose_whole_crystal: summary.avg_dose_whole_crystal(),
            avg_dose_exposed_region: summary.avg_dose_exposed_region(),
            elastic_yield: summary.wedge_elastic(),
            absorbed_energy: summary.abs_energy_total(),
            used_volume_fraction: summary.used_volume_fraction(),
        });
    }
}
