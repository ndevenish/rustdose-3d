//! WebAssembly bindings for RADDOSE-3D.
//!
//! Build for the browser:
//! ```sh
//! wasm-pack build --target web
//! ```
//!
//! Build for Node.js:
//! ```sh
//! wasm-pack build --target nodejs
//! ```
//!
//! The resulting npm package can be published as `@raddose3d/wasm`.

use wasm_bindgen::prelude::*;

/// Run a RADDOSE-3D simulation from a text input string.
///
/// Accepts the standard RADDOSE-3D `.txt` input format.
/// Returns the summary text that would normally be written to stdout.
///
/// On error, throws a JavaScript exception with the error message.
#[wasm_bindgen]
pub fn run(input: &str) -> Result<String, JsValue> {
    use raddose3d::beam;
    use raddose3d::crystal;
    use raddose3d::experiment::Experiment;
    use raddose3d::output::OutputSummaryText;
    use raddose3d::wedge::Wedge;
    use raddose3d::writer::shared_string_writer;

    let config = raddose3d::parse_input(input).map_err(|e| JsValue::from_str(&e))?;

    // Capture summary text via a shared buffer that lives for 'static.
    let (writer, buf_arc) = shared_string_writer();
    let mut exp = Experiment::new();
    exp.add_observer(Box::new(OutputSummaryText::new(writer)));

    for crystal_config in &config.crystals {
        let crystal = crystal::create_crystal(crystal_config).map_err(|e| JsValue::from_str(&e))?;
        exp.set_crystal(crystal);
    }
    for beam_config in &config.beams {
        let beam = beam::create_beam(beam_config).map_err(|e| JsValue::from_str(&e))?;
        exp.set_beam(beam);
    }
    for wedge_config in &config.wedges {
        let wedge = Wedge::from_config(wedge_config);
        exp.expose_wedge(&wedge);
    }
    exp.close();

    let bytes = buf_arc.lock().unwrap().clone();
    String::from_utf8(bytes).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Run a RADDOSE-3D simulation and return structured results as a JS object.
///
/// The returned object has the following numeric fields (all in MGy unless noted):
/// - `averageDwd` — average diffraction-weighted dose
/// - `lastDwd` — last diffraction-weighted dose
/// - `maxDose` — maximum voxel dose
/// - `avgDoseWholeCrystal` — average dose, whole crystal
/// - `avgDoseExposedRegion` — average dose, exposed region
/// - `elasticYield` — total elastic photon count (dimensionless)
/// - `absorbedEnergy` — total absorbed energy (Joules)
/// - `usedVolumeFraction` — fraction of crystal voxels that received dose (0–100)
///
/// On error, throws a JavaScript exception.
#[wasm_bindgen(js_name = runStructured)]
pub fn run_structured(input: &str) -> Result<JsValue, JsValue> {
    let config = raddose3d::parse_input(input).map_err(|e| JsValue::from_str(&e))?;
    let results = raddose3d::run(&config).map_err(|e| JsValue::from_str(&e))?;

    let obj = serde_wasm_bindgen::to_value(&WasmResults {
        average_dwd: results.average_dwd,
        last_dwd: results.last_dwd,
        max_dose: results.max_dose,
        avg_dose_whole_crystal: results.avg_dose_whole_crystal,
        avg_dose_exposed_region: results.avg_dose_exposed_region,
        elastic_yield: results.elastic_yield,
        absorbed_energy: results.absorbed_energy,
        used_volume_fraction: results.used_volume_fraction,
    })
    .map_err(|e| JsValue::from_str(&e.to_string()))?;

    Ok(obj)
}

/// Run from a JSON-encoded config rather than the `.txt` input format.
///
/// The JSON must match the structure produced by `parse_input` serialised to
/// JSON — i.e. the same schema as [`raddose3d::parser::Config`].
#[wasm_bindgen(js_name = runFromJson)]
pub fn run_from_json(json: &str) -> Result<JsValue, JsValue> {
    let config = raddose3d::parse_input_json(json).map_err(|e| JsValue::from_str(&e))?;
    let results = raddose3d::run(&config).map_err(|e| JsValue::from_str(&e))?;

    let obj = serde_wasm_bindgen::to_value(&WasmResults {
        average_dwd: results.average_dwd,
        last_dwd: results.last_dwd,
        max_dose: results.max_dose,
        avg_dose_whole_crystal: results.avg_dose_whole_crystal,
        avg_dose_exposed_region: results.avg_dose_exposed_region,
        elastic_yield: results.elastic_yield,
        absorbed_energy: results.absorbed_energy,
        used_volume_fraction: results.used_volume_fraction,
    })
    .map_err(|e| JsValue::from_str(&e.to_string()))?;

    Ok(obj)
}

/// Run a full simulation and return scalar results, voxel dose grid, and dose histogram.
///
/// Returns a JS object combining all fields from `runStructured()` plus:
/// - `crystalSizeUm` — `[x, y, z]` crystal dimensions in µm
/// - `voxelX`, `voxelY`, `voxelZ` — flat arrays of voxel coordinates in µm
/// - `voxelDose` — flat array of voxel doses in MGy (0 for non-crystal voxels)
/// - `doseHist` — 11-element array of normalised voxel fractions by dose bin
///
/// The grid is downsampled to at most 91 voxels per axis (same as the HTML report).
///
/// On error, throws a JavaScript exception.
#[wasm_bindgen(js_name = runFull)]
pub fn run_full(input: &str) -> Result<JsValue, JsValue> {
    use raddose3d::beam;
    use raddose3d::crystal;
    use raddose3d::experiment::Experiment;
    use raddose3d::wedge::Wedge;
    use std::sync::{Arc, Mutex};

    let config = raddose3d::parse_input(input).map_err(|e| JsValue::from_str(&e))?;

    let shared = Arc::new(Mutex::new(FullData::default()));
    let mut exp = Experiment::new();
    exp.add_observer(Box::new(FullCollector(Arc::clone(&shared))));

    for crystal_config in &config.crystals {
        let crystal = crystal::create_crystal(crystal_config).map_err(|e| JsValue::from_str(&e))?;
        exp.set_crystal(crystal);
    }
    for beam_config in &config.beams {
        let beam = beam::create_beam(beam_config).map_err(|e| JsValue::from_str(&e))?;
        exp.set_beam(beam);
    }
    for wedge_config in &config.wedges {
        let wedge = Wedge::from_config(wedge_config);
        exp.expose_wedge(&wedge);
    }
    exp.close();

    let data = Arc::try_unwrap(shared)
        .map_err(|_| JsValue::from_str("internal error: collector Arc still shared"))?
        .into_inner()
        .map_err(|_| JsValue::from_str("internal error: collector mutex poisoned"))?;

    let results = data
        .results
        .ok_or_else(|| JsValue::from_str("No wedge was exposed"))?;

    let obj = serde_wasm_bindgen::to_value(&WasmFullResults {
        average_dwd: results.average_dwd,
        last_dwd: results.last_dwd,
        max_dose: results.max_dose,
        avg_dose_whole_crystal: results.avg_dose_whole_crystal,
        avg_dose_exposed_region: results.avg_dose_exposed_region,
        elastic_yield: results.elastic_yield,
        absorbed_energy: results.absorbed_energy,
        used_volume_fraction: results.used_volume_fraction,
        crystal_size_um: data.crystal_size_um,
        voxel_x: data.voxel_x,
        voxel_y: data.voxel_y,
        voxel_z: data.voxel_z,
        voxel_dose: data.voxel_dose,
        dose_hist: data.dose_hist.to_vec(),
    })
    .map_err(|e| JsValue::from_str(&e.to_string()))?;

    Ok(obj)
}

// ── Internal: full-results collector ─────────────────────────────────────────

/// Maximum voxels per axis before downsampling (matches html_report.rs).
const MAX_VOXELS_PER_AXIS: usize = 91;

#[derive(Debug, Default)]
struct FullData {
    results: Option<WasmResults>,
    dose_hist: [f64; 11],
    crystal_size_um: [f64; 3],
    voxel_x: Vec<f32>,
    voxel_y: Vec<f32>,
    voxel_z: Vec<f32>,
    voxel_dose: Vec<f32>,
}

#[derive(Debug)]
struct FullCollector(std::sync::Arc<std::sync::Mutex<FullData>>);

impl raddose3d::output::Output for FullCollector {
    fn publish_crystal(&mut self, _: &dyn raddose3d::crystal::Crystal) {}
    fn publish_beam(&mut self, _: &dyn raddose3d::beam::Beam) {}

    fn publish_wedge(
        &mut self,
        _wedge: &raddose3d::wedge::Wedge,
        summary: &raddose3d::output::ExposureSummary,
        _crystal: Option<&dyn raddose3d::crystal::Crystal>,
    ) {
        let mut data = self.0.lock().unwrap();
        data.results = Some(WasmResults {
            average_dwd: summary.avg_diffracted_dose(),
            last_dwd: summary.last_dwd(),
            max_dose: summary.max_dose(),
            avg_dose_whole_crystal: summary.avg_dose_whole_crystal(),
            avg_dose_exposed_region: summary.avg_dose_exposed_region(),
            elastic_yield: summary.wedge_elastic(),
            absorbed_energy: summary.abs_energy_total(),
            used_volume_fraction: summary.used_volume_fraction(),
        });
        data.dose_hist = summary.dose_hist_normalised();
    }

    fn close(&mut self, crystal: Option<&dyn raddose3d::crystal::Crystal>) {
        let Some(c) = crystal else { return };

        let size = c.cryst_size_voxels();
        let size_um = c.cryst_size_um();

        let step = [
            size[0].div_ceil(MAX_VOXELS_PER_AXIS),
            size[1].div_ceil(MAX_VOXELS_PER_AXIS),
            size[2].div_ceil(MAX_VOXELS_PER_AXIS),
        ];
        let nx = size[0].div_ceil(step[0]);
        let ny = size[1].div_ceil(step[1]);
        let nz = size[2].div_ceil(step[2]);

        let total = nx * ny * nz;
        let mut xs = Vec::with_capacity(total);
        let mut ys = Vec::with_capacity(total);
        let mut zs = Vec::with_capacity(total);
        let mut vals = Vec::with_capacity(total);

        for ii in 0..nx {
            let i = (ii * step[0]).min(size[0] - 1);
            let x = if nx > 1 {
                size_um[0] * ii as f64 / (nx - 1) as f64
            } else {
                0.0
            };
            for jj in 0..ny {
                let j = (jj * step[1]).min(size[1] - 1);
                let y = if ny > 1 {
                    size_um[1] * jj as f64 / (ny - 1) as f64
                } else {
                    0.0
                };
                for kk in 0..nz {
                    let k = (kk * step[2]).min(size[2] - 1);
                    let z = if nz > 1 {
                        size_um[2] * kk as f64 / (nz - 1) as f64
                    } else {
                        0.0
                    };
                    xs.push(x as f32);
                    ys.push(y as f32);
                    zs.push(z as f32);
                    let dose = if c.is_crystal_at(i, j, k) {
                        c.get_dose(i, j, k)
                    } else {
                        0.0
                    };
                    vals.push(dose as f32);
                }
            }
        }

        let mut data = self.0.lock().unwrap();
        data.crystal_size_um = size_um;
        data.voxel_x = xs;
        data.voxel_y = ys;
        data.voxel_z = zs;
        data.voxel_dose = vals;
    }
}

// ── Serialisable result structs ───────────────────────────────────────────────

#[derive(Debug, serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct WasmResults {
    average_dwd: f64,
    last_dwd: f64,
    max_dose: f64,
    avg_dose_whole_crystal: f64,
    avg_dose_exposed_region: f64,
    elastic_yield: f64,
    absorbed_energy: f64,
    used_volume_fraction: f64,
}

#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct WasmFullResults {
    average_dwd: f64,
    last_dwd: f64,
    max_dose: f64,
    avg_dose_whole_crystal: f64,
    avg_dose_exposed_region: f64,
    elastic_yield: f64,
    absorbed_energy: f64,
    used_volume_fraction: f64,
    crystal_size_um: [f64; 3],
    voxel_x: Vec<f32>,
    voxel_y: Vec<f32>,
    voxel_z: Vec<f32>,
    voxel_dose: Vec<f32>,
    dose_hist: Vec<f64>,
}
