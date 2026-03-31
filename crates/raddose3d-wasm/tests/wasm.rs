use wasm_bindgen_test::*;

wasm_bindgen_test_configure!(run_in_node_experimental);

const INSULIN_INPUT: &str = "
Crystal
Type Cuboid
Dimensions 100 100 100
PixelsPerMicron 0.5
AbsCoefCalc RD3D
UnitCell 78.0 78.0 78.0
NumMonomers 24
NumResidues 51

Beam
Type Gaussian
Flux 2e12
FWHM 100 100
Energy 12.1
Collimation Rectangular 100 100

Wedge 0 360
ExposureTime 100
AngularResolution 2
";

/// run() returns non-empty summary text and contains expected keywords.
#[wasm_bindgen_test]
fn test_run_returns_summary_text() {
    let result = raddose3d_wasm::run(INSULIN_INPUT);
    assert!(result.is_ok(), "run() failed: {:?}", result.err());
    let text = result.unwrap();
    assert!(!text.is_empty(), "summary text is empty");
    assert!(
        text.contains("Max Dose") || text.contains("DWD") || text.contains("dose"),
        "summary text does not mention dose: {text}"
    );
}

/// runStructured() returns a JS object with the expected numeric fields.
#[wasm_bindgen_test]
fn test_run_structured_fields() {
    use js_sys::Reflect;
    use wasm_bindgen::JsValue;

    let result = raddose3d_wasm::run_structured(INSULIN_INPUT);
    assert!(result.is_ok(), "runStructured() failed: {:?}", result.err());
    let obj = result.unwrap();

    let fields = [
        "averageDwd",
        "lastDwd",
        "maxDose",
        "avgDoseWholeCrystal",
        "avgDoseExposedRegion",
        "elasticYield",
        "absorbedEnergy",
        "usedVolumeFraction",
    ];
    for field in &fields {
        let val = Reflect::get(&obj, &JsValue::from_str(field)).unwrap_or(JsValue::UNDEFINED);
        assert!(
            val.as_f64().is_some(),
            "field '{field}' is missing or not a number"
        );
    }
}

/// runStructured() values are physically reasonable for the insulin test case.
#[wasm_bindgen_test]
fn test_run_structured_values_reasonable() {
    use js_sys::Reflect;
    use wasm_bindgen::JsValue;

    let obj = raddose3d_wasm::run_structured(INSULIN_INPUT).unwrap();

    let get = |field: &str| -> f64 {
        Reflect::get(&obj, &JsValue::from_str(field))
            .unwrap()
            .as_f64()
            .unwrap()
    };

    let max_dose = get("maxDose");
    let avg_dwd = get("averageDwd");
    let elastic_yield = get("elasticYield");
    let used_vol = get("usedVolumeFraction");

    // Validated against Java: max ~12 MGy, avg DWD ~4 MGy, elastic ~3.6e11
    assert!(
        max_dose > 5.0 && max_dose < 25.0,
        "maxDose out of range: {max_dose}"
    );
    assert!(
        avg_dwd > 1.0 && avg_dwd < 15.0,
        "averageDwd out of range: {avg_dwd}"
    );
    assert!(
        elastic_yield > 1e10,
        "elasticYield too low: {elastic_yield}"
    );
    assert!(
        used_vol > 50.0 && used_vol <= 100.0,
        "usedVolumeFraction out of range: {used_vol}"
    );
}

/// Malformed input returns a JS error, not a panic.
#[wasm_bindgen_test]
fn test_run_bad_input_returns_error() {
    let result = raddose3d_wasm::run("this is not valid raddose input @@@@");
    assert!(result.is_err(), "expected error for bad input, got Ok");
}

/// runStructured() and run() agree that the simulation ran (cross-check).
#[wasm_bindgen_test]
fn test_run_and_run_structured_consistent() {
    use js_sys::Reflect;
    use wasm_bindgen::JsValue;

    let text = raddose3d_wasm::run(INSULIN_INPUT).unwrap();
    let obj = raddose3d_wasm::run_structured(INSULIN_INPUT).unwrap();

    let max_dose = Reflect::get(&obj, &JsValue::from_str("maxDose"))
        .unwrap()
        .as_f64()
        .unwrap();

    // Both succeeded and max dose is positive
    assert!(!text.is_empty());
    assert!(max_dose > 0.0);
}
