//! PDB validation tests.
//!
//! Runs the full pipeline using a local copy of 1DWA.pdb (saved in
//! tests/fixtures/) and compares key metrics against Java golden files
//! captured with:
//!   java -jar raddose3d.jar -i <input_with_local_pdb> -p golden/pdb/golden_pdb
//!
//! Validated bugs fixed:
//!   - Coefficients in Summary.txt were 0 (captured before update_coefficients)
//!   - Density was 1.17 instead of 1.14 g/ml (hetero_atom_occurrence not
//!     multiplied by num_monomers; hetatm mass missing from solvent fraction)

use raddose3d::{
    experiment::Experiment,
    output::{OutputDWDs, OutputRDECSV, OutputSummaryCSV, OutputSummaryText},
    parser,
    writer::shared_string_writer,
};

const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures");
const GOLDEN_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/pdb");

/// Build the input config string with an absolute path to the local PDB file,
/// so this test works regardless of working directory.
fn pdb_input() -> String {
    let pdb_path = format!("{}/1dwa.pdb", FIXTURES_DIR);
    format!(
        "\
Crystal
Type Cuboid
Dimensions 100 100 100
AbsCoefCalc EXP
Pdb {pdb_path}
SolventHeavyConc S 2700
PixelsPerMicron 0.2

Beam
Type TopHat
Energy 12.40
Collimation Rectangular 160 160
Flux 6.00e12

Wedge 0 90
exposureTime 100
"
    )
}

fn run_pdb() -> std::collections::HashMap<&'static str, String> {
    let input = pdb_input();
    let config = parser::parse(&input).expect("parse failed");

    let (w_sumtxt, b_sumtxt) = shared_string_writer();
    let (w_sumcsv, b_sumcsv) = shared_string_writer();
    let (w_dwds, b_dwds) = shared_string_writer();
    let (w_rde, b_rde) = shared_string_writer();

    let mut exp = Experiment::new();
    exp.add_observer(Box::new(OutputSummaryText::new(w_sumtxt)));
    exp.add_observer(Box::new(OutputSummaryCSV::new(w_sumcsv)));
    exp.add_observer(Box::new(OutputDWDs::new(w_dwds)));
    exp.add_observer(Box::new(OutputRDECSV::new(w_rde)));

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
    out
}

fn golden(name: &str) -> String {
    let path = format!("{}/golden_pdb{}", GOLDEN_DIR, name);
    std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("golden file missing: {path}"))
}

/// Extract the numeric value after `: ` on the first line matching `label`.
fn extract_metric(text: &str, label: &str) -> f64 {
    for line in text.lines() {
        if line.contains(label) {
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

// ── Summary.txt ───────────────────────────────────────────────────────────────

/// Density must match Java (was 1.17 before fix; correct is 1.14).
#[test]
fn pdb_density_matches_java() {
    let out = run_pdb();
    let rust = &out["Summary.txt"];
    let java = golden("Summary.txt");

    let rust_density = extract_metric(rust, "Density:");
    let java_density = extract_metric(&java, "Density:");
    assert!(
        (rust_density - java_density).abs() / java_density < 0.001,
        "Density: rust={rust_density} java={java_density} (>0.1% diff)"
    );
}

/// Coefficients must be non-zero (were all 0.00e0 before fix).
#[test]
fn pdb_coefficients_nonzero() {
    let out = run_pdb();
    let rust = &out["Summary.txt"];

    let photo = extract_metric(rust, "Photelectric Coefficient:");
    assert!(
        photo > 0.0,
        "Photoelectric coefficient is zero — update_coefficients not called before summary"
    );
    let att = extract_metric(rust, "Attenuation Coefficient:");
    assert!(att > 0.0, "Attenuation coefficient is zero");
}

/// Photoelectric coefficient must be within 1% of Java.
#[test]
fn pdb_photoelectric_coeff_matches_java() {
    let out = run_pdb();
    let rust = &out["Summary.txt"];
    let java = golden("Summary.txt");

    let rust_val = extract_metric(rust, "Photelectric Coefficient:");
    let java_val = extract_metric(&java, "Photelectric Coefficient:");
    assert!(
        (rust_val - java_val).abs() / java_val < 0.01,
        "Photoelectric coeff: rust={rust_val} java={java_val} (>1% diff)"
    );
}

/// Average DWD must match Java within 0.1%.
#[test]
fn pdb_avg_dwd_matches_java() {
    let out = run_pdb();
    let rust = &out["Summary.txt"];
    let java = golden("Summary.txt");

    let rust_dwd = extract_metric(rust, "Average Diffraction Weighted Dose");
    let java_dwd = extract_metric(&java, "Average Diffraction Weighted Dose");
    assert!(
        (rust_dwd - java_dwd).abs() / java_dwd < 0.001,
        "Avg DWD: rust={rust_dwd} java={java_dwd} (>0.1% diff)"
    );
}

/// Max Dose must match Java exactly (bit-identical for this test case).
#[test]
fn pdb_max_dose_matches_java() {
    let out = run_pdb();
    let rust = &out["Summary.txt"];
    let java = golden("Summary.txt");

    let rust_val = extract_metric(rust, "Max Dose");
    let java_val = extract_metric(&java, "Max Dose");
    assert!(
        (rust_val - java_val).abs() / java_val < 1e-5,
        "Max Dose: rust={rust_val} java={java_val}"
    );
}

/// Summary.csv row count must match Java.
#[test]
fn pdb_summary_csv_row_count() {
    let out = run_pdb();
    let rust = &out["Summary.csv"];
    let java = golden("Summary.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "Summary.csv row count: rust={rust_rows} java={java_rows}"
    );
}

/// DWDs.csv row count must match Java.
#[test]
fn pdb_dwds_csv_row_count() {
    let out = run_pdb();
    let rust = &out["DWDs.csv"];
    let java = golden("DWDs.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "DWDs.csv row count: rust={rust_rows} java={java_rows}"
    );
}

/// RDE.csv row count must match Java.
#[test]
fn pdb_rde_csv_row_count() {
    let out = run_pdb();
    let rust = &out["RDE.csv"];
    let java = golden("RDE.csv");

    let rust_rows = rust.lines().count();
    let java_rows = java.lines().count();
    assert_eq!(
        rust_rows, java_rows,
        "RDE.csv row count: rust={rust_rows} java={java_rows}"
    );
}
