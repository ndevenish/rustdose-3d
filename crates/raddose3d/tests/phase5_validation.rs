//! Phase 5 validation tests.
//!
//! Runs the full pipeline with all output types enabled and compares
//! against Java golden files captured with:
//!   java -jar raddose3d.jar -i insulin_test.txt -p golden/phase5/insulin_

use raddose3d::{
    experiment::Experiment,
    output::{
        OutputDWDs, OutputFinalDoseStateCSV, OutputFinalDoseStateR, OutputRDECSV, OutputSummaryCSV,
        OutputSummaryText, OutputVoxelDose, OutputVoxelFluences,
    },
    parser,
    writer::shared_string_writer,
};

const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures");
const GOLDEN_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/phase5");

fn run_insulin() -> std::collections::HashMap<&'static str, String> {
    let path = format!("{}/insulin_test.txt", FIXTURES_DIR);
    let input = std::fs::read_to_string(&path).expect("insulin fixture missing");
    let config = parser::parse(&input).expect("parse failed");

    // Create shared string writers for each output type
    let (w_sumtxt, b_sumtxt) = shared_string_writer();
    let (w_sumcsv, b_sumcsv) = shared_string_writer();
    let (w_dwds, b_dwds) = shared_string_writer();
    let (w_rde, b_rde) = shared_string_writer();
    let (w_dosecsv, b_dosecsv) = shared_string_writer();
    let (w_dose_r, b_dose_r) = shared_string_writer();
    let (w_voxdose, b_voxdose) = shared_string_writer();
    let (w_voxflu, b_voxflu) = shared_string_writer();

    let mut exp = Experiment::new();
    exp.add_observer(Box::new(OutputSummaryText::new(w_sumtxt)));
    exp.add_observer(Box::new(OutputSummaryCSV::new(w_sumcsv)));
    exp.add_observer(Box::new(OutputDWDs::new(w_dwds)));
    exp.add_observer(Box::new(OutputRDECSV::new(w_rde)));
    exp.add_observer(Box::new(OutputFinalDoseStateCSV::new(w_dosecsv)));
    exp.add_observer(Box::new(OutputFinalDoseStateR::new(w_dose_r)));
    exp.add_observer(Box::new(OutputVoxelDose::new(w_voxdose)));
    exp.add_observer(Box::new(OutputVoxelFluences::new(w_voxflu)));

    // Wire up crystals/beams/wedges
    for cc in &config.crystals {
        let c = raddose3d::crystal::create_crystal(cc).unwrap();
        exp.set_crystal(c);
    }
    for bc in &config.beams {
        let b = raddose3d::beam::create_beam(bc).unwrap();
        exp.set_beam(b);
    }
    for wc in &config.wedges {
        let w = raddose3d::wedge::Wedge::from_config(wc);
        exp.expose_wedge(&w);
    }
    exp.close();

    let to_string = |buf: std::sync::Arc<std::sync::Mutex<Vec<u8>>>| {
        String::from_utf8(buf.lock().unwrap().clone()).unwrap()
    };

    let mut out = std::collections::HashMap::new();
    out.insert("Summary.txt", to_string(b_sumtxt));
    out.insert("Summary.csv", to_string(b_sumcsv));
    out.insert("DWDs.csv", to_string(b_dwds));
    out.insert("RDE.csv", to_string(b_rde));
    out.insert("DoseState.csv", to_string(b_dosecsv));
    out.insert("DoseState.R", to_string(b_dose_r));
    out.insert("VoxelDose.csv", to_string(b_voxdose));
    out.insert("VoxelFluences.csv", to_string(b_voxflu));
    out
}

fn golden(name: &str) -> String {
    let path = format!("{}/insulin_{}", GOLDEN_DIR, name);
    std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("golden file missing: {}", path))
}

// ── Summary.txt ───────────────────────────────────────────────────────────────

#[test]
fn phase5_summary_txt_avg_dwd() {
    let out = run_insulin();
    let rust = &out["Summary.txt"];
    let java = golden("Summary.txt");

    // Extract Average DWD from both
    let rust_dwd = extract_metric(rust, "Average Diffraction Weighted Dose");
    let java_dwd = extract_metric(&java, "Average Diffraction Weighted Dose");
    assert!(
        (rust_dwd - java_dwd).abs() / java_dwd < 0.001,
        "Avg DWD: rust={rust_dwd} java={java_dwd} (>0.1% diff)"
    );
}

#[test]
fn phase5_summary_txt_max_dose() {
    let out = run_insulin();
    let rust = &out["Summary.txt"];
    let java = golden("Summary.txt");

    let rust_val = extract_metric(rust, "Max Dose");
    let java_val = extract_metric(&java, "Max Dose");
    assert!(
        (rust_val - java_val).abs() / java_val < 0.001,
        "Max Dose: rust={rust_val} java={java_val}"
    );
}

#[test]
fn phase5_summary_txt_used_volume() {
    let out = run_insulin();
    let rust = &out["Summary.txt"];
    let java = golden("Summary.txt");

    assert!(
        rust.contains("Used Volume") && rust.contains("100.0%"),
        "Used Volume 100% not found in Rust output"
    );
    assert!(
        java.contains("100.0%"),
        "Used Volume 100% not found in Java golden"
    );
}

#[test]
fn phase5_summary_txt_histogram_present() {
    let out = run_insulin();
    let rust = &out["Summary.txt"];
    assert!(
        rust.contains("Final Dose Histogram:"),
        "Dose histogram missing from Summary.txt"
    );
    // Should have 11 bins (Bin 1 through Bin 11)
    assert!(
        rust.contains("Bin 11,"),
        "Bin 11 missing — histogram incomplete"
    );
}

#[test]
fn phase5_summary_txt_dose_inefficiency_pe() {
    let out = run_insulin();
    let rust = &out["Summary.txt"];
    assert!(
        rust.contains("Dose Inefficiency PE"),
        "Dose Inefficiency PE line missing"
    );
}

#[test]
fn phase5_summary_txt_crystal_coefcalc_info() {
    let out = run_insulin();
    let rust = &out["Summary.txt"];
    assert!(
        rust.contains("Crystal coefficients calculated with RADDOSE-3D"),
        "Coefficient info missing from Summary.txt"
    );
    assert!(
        rust.contains("Photelectric Coefficient:"),
        "Photelectric Coefficient line missing"
    );
}

// ── Summary.csv ───────────────────────────────────────────────────────────────

#[test]
fn phase5_summary_csv_header() {
    let out = run_insulin();
    let rust = &out["Summary.csv"];
    let java = golden("Summary.csv");

    let rust_header = rust.lines().next().unwrap_or("");
    let java_header = java.lines().next().unwrap_or("");
    assert_eq!(
        rust_header, java_header,
        "Summary.csv header mismatch\nrust: {rust_header}\njava: {java_header}"
    );
}

#[test]
fn phase5_summary_csv_row_count() {
    let out = run_insulin();
    let rust = &out["Summary.csv"];
    let java = golden("Summary.csv");
    // Header + 1 wedge row
    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "Summary.csv row count mismatch: rust={rust_rows} java={java_rows}"
    );
}

// ── DWDs.csv ─────────────────────────────────────────────────────────────────

#[test]
fn phase5_dwds_csv_header() {
    let out = run_insulin();
    let rust = &out["DWDs.csv"];
    let java = golden("DWDs.csv");

    let rust_hdr = rust.lines().next().unwrap_or("").trim();
    let java_hdr = java.lines().next().unwrap_or("").trim();
    assert_eq!(
        rust_hdr, java_hdr,
        "DWDs.csv header mismatch\nrust: {rust_hdr}\njava: {java_hdr}"
    );
}

#[test]
fn phase5_dwds_csv_row_count() {
    let out = run_insulin();
    let rust = &out["DWDs.csv"];
    let java = golden("DWDs.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "DWDs.csv row count mismatch: rust={rust_rows} java={java_rows}"
    );
}

#[test]
fn phase5_dwds_csv_first_row_values() {
    let out = run_insulin();
    let rust = &out["DWDs.csv"];
    let java = golden("DWDs.csv");

    let rust_first = rust.lines().nth(1).unwrap_or("");
    let java_first = java.lines().nth(1).unwrap_or("");

    let rust_vals: Vec<f64> = rust_first
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    let java_vals: Vec<f64> = java_first
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    assert_eq!(
        rust_vals.len(),
        java_vals.len(),
        "DWDs first row column count mismatch"
    );
    for (i, (r, j)) in rust_vals.iter().zip(java_vals.iter()).enumerate() {
        let tol = if *j == 0.0 { 1e-10 } else { (*j).abs() * 0.01 };
        assert!(
            (r - j).abs() <= tol,
            "DWDs.csv first row col {i}: rust={r} java={j} (>1% diff)"
        );
    }
}

// ── RDE.csv ──────────────────────────────────────────────────────────────────

#[test]
fn phase5_rde_csv_header() {
    let out = run_insulin();
    let rust = &out["RDE.csv"];
    let java = golden("RDE.csv");

    let rust_hdr = rust.lines().next().unwrap_or("").trim();
    let java_hdr = java.lines().next().unwrap_or("").trim();
    assert_eq!(
        rust_hdr, java_hdr,
        "RDE.csv header mismatch\nrust: {rust_hdr}\njava: {java_hdr}"
    );
}

#[test]
fn phase5_rde_csv_row_count() {
    let out = run_insulin();
    let rust = &out["RDE.csv"];
    let java = golden("RDE.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "RDE.csv row count mismatch: rust={rust_rows} java={java_rows}"
    );
}

// ── DoseState.csv ─────────────────────────────────────────────────────────────

#[test]
fn phase5_dose_state_csv_row_count() {
    let out = run_insulin();
    let rust = &out["DoseState.csv"];
    let java = golden("DoseState.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "DoseState.csv row count: rust={rust_rows} java={java_rows}"
    );
}

#[test]
fn phase5_dose_state_csv_columns_per_row() {
    let out = run_insulin();
    let rust = &out["DoseState.csv"];
    // Each row should have 6 comma-separated values: x,y,z,dose,fluence,elastic
    for (i, line) in rust.lines().enumerate() {
        let cols = line.split(',').count();
        assert_eq!(
            cols, 6,
            "DoseState.csv row {i}: expected 6 columns, got {cols}"
        );
    }
}

// ── DoseState.R ───────────────────────────────────────────────────────────────

#[test]
fn phase5_dose_state_r_structure() {
    let out = run_insulin();
    let rust = &out["DoseState.R"];
    assert!(
        rust.contains("require(\"rgl\")"),
        "DoseState.R missing rgl require"
    );
    assert!(
        rust.contains("dose <- array("),
        "DoseState.R missing array declaration"
    );
    assert!(
        rust.contains("contour3d("),
        "DoseState.R missing contour3d call"
    );
}

// ── VoxelDose.csv / VoxelFluences.csv ────────────────────────────────────────

#[test]
fn phase5_voxel_dose_row_count() {
    let out = run_insulin();
    let rust = &out["VoxelDose.csv"];
    let java = golden("VoxelDose.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "VoxelDose.csv row count: rust={rust_rows} java={java_rows}"
    );
}

#[test]
fn phase5_voxel_fluences_row_count() {
    let out = run_insulin();
    let rust = &out["VoxelFluences.csv"];
    let java = golden("VoxelFluences.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "VoxelFluences.csv row count: rust={rust_rows} java={java_rows}"
    );
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Extract the numeric value after `: ` on the first line matching `label`.
fn extract_metric(text: &str, label: &str) -> f64 {
    for line in text.lines() {
        if line.contains(label) {
            // Find the last numeric token on the line
            if let Some(val) = line
                .split_whitespace()
                .filter_map(|s| s.trim_matches(',').parse::<f64>().ok())
                .next_back()
            {
                return val;
            }
        }
    }
    panic!("label '{label}' not found in output")
}
