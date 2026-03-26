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

// ── Internal serialisable results struct ─────────────────────────────────────

#[derive(serde::Serialize)]
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
