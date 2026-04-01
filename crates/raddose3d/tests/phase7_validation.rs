//! Phase 7 validation tests.
//!
//! Covers the public library API (`parse_input`, `parse_input_json`, `run`),
//! CLI binary behaviour (spawned as a subprocess), and WASM-crate logic
//! exercised at the native Rust level.

use raddose3d::{parse_input, parse_input_json, run};

const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures");
const GOLDEN_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/insulin");

// ── Helpers ───────────────────────────────────────────────────────────────────

fn insulin_input() -> String {
    std::fs::read_to_string(format!("{}/insulin_test.txt", FIXTURES_DIR))
        .expect("insulin fixture missing")
}

fn golden(name: &str) -> String {
    std::fs::read_to_string(format!("{}/insulin_{}", GOLDEN_DIR, name))
        .unwrap_or_else(|_| panic!("golden file missing: {}", name))
}

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

fn cli_bin() -> std::path::PathBuf {
    // Use the debug binary built by cargo.
    let mut p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("../../../target/debug/raddose3d");
    p
}

// ── parse_input ───────────────────────────────────────────────────────────────

#[test]
fn parse_input_returns_same_as_parser_crate() {
    let text = insulin_input();
    let via_lib = parse_input(&text).expect("parse_input failed");
    let via_parser = raddose3d::parser::parse(&text).expect("parser::parse failed");
    assert_eq!(via_lib.crystals.len(), via_parser.crystals.len());
    assert_eq!(via_lib.beams.len(), via_parser.beams.len());
    assert_eq!(via_lib.wedges.len(), via_parser.wedges.len());
}

#[test]
fn parse_input_error_on_garbage() {
    assert!(parse_input("this is not valid input @@@").is_err());
}

// ── parse_input_json ──────────────────────────────────────────────────────────

#[test]
fn parse_input_json_roundtrip() {
    // Parse the text format → serialise to JSON → deserialise back.
    let text = insulin_input();
    let config = parse_input(&text).expect("parse_input failed");
    let json = serde_json::to_string(&config).expect("serialise to JSON failed");
    let config2 = parse_input_json(&json).expect("parse_input_json failed");

    assert_eq!(config.crystals.len(), config2.crystals.len());
    assert_eq!(config.beams.len(), config2.beams.len());
    assert_eq!(config.wedges.len(), config2.wedges.len());

    // Key numeric fields survive the round-trip.
    let c1 = &config.crystals[0];
    let c2 = &config2.crystals[0];
    assert_eq!(c1.dim_x, c2.dim_x);
    assert_eq!(c1.dim_y, c2.dim_y);
    assert_eq!(c1.dim_z, c2.dim_z);
    assert_eq!(c1.pixels_per_micron, c2.pixels_per_micron);

    let b1 = &config.beams[0];
    let b2 = &config2.beams[0];
    assert_eq!(b1.flux, b2.flux);
    assert_eq!(b1.energy, b2.energy);

    let w1 = &config.wedges[0];
    let w2 = &config2.wedges[0];
    assert_eq!(w1.start_ang, w2.start_ang);
    assert_eq!(w1.end_ang, w2.end_ang);
}

#[test]
fn parse_input_json_error_on_garbage() {
    assert!(parse_input_json("{not: valid json}").is_err());
}

// ── run() / RunResults ────────────────────────────────────────────────────────

#[test]
fn run_returns_results_for_insulin() {
    let config = parse_input(&insulin_input()).unwrap();
    let results = run(&config).expect("run failed");

    // Values must be physically reasonable.
    assert!(results.average_dwd > 0.0, "average_dwd must be positive");
    assert!(results.max_dose > 0.0, "max_dose must be positive");
    assert!(results.max_dose > results.average_dwd, "max > avg");
    assert!(
        results.used_volume_fraction > 90.0,
        "insulin should fill > 90% of crystal"
    );
    assert!(results.elastic_yield > 0.0, "elastic yield must be > 0");
    assert!(results.absorbed_energy > 0.0, "absorbed energy must be > 0");
}

#[test]
fn run_avg_dwd_matches_golden() {
    let config = parse_input(&insulin_input()).unwrap();
    let results = run(&config).expect("run failed");

    let java_summary = golden("Summary.txt");
    let java_dwd = extract_metric(&java_summary, "Average Diffraction Weighted Dose");

    assert!(
        (results.average_dwd - java_dwd).abs() / java_dwd < 0.001,
        "run() average_dwd {:.6} differs from Java golden {:.6} by > 0.1%",
        results.average_dwd,
        java_dwd
    );
}

#[test]
fn run_max_dose_matches_golden() {
    let config = parse_input(&insulin_input()).unwrap();
    let results = run(&config).expect("run failed");

    let java_summary = golden("Summary.txt");
    let java_max = extract_metric(&java_summary, "Max Dose");

    assert!(
        (results.max_dose - java_max).abs() / java_max < 0.001,
        "run() max_dose {:.6} differs from Java golden {:.6} by > 0.1%",
        results.max_dose,
        java_max
    );
}

#[test]
fn run_no_wedge_returns_error() {
    use raddose3d::parser::{BeamConfig, Config, CrystalConfig};
    let config = Config {
        crystals: vec![CrystalConfig::default()],
        beams: vec![BeamConfig::default()],
        wedges: vec![],
    };
    // run() with no wedges should return Err because nothing is exposed.
    assert!(run(&config).is_err(), "expected error with no wedges");
}

// ── CLI binary tests ──────────────────────────────────────────────────────────

#[test]
fn cli_version_flag() {
    let bin = cli_bin();
    if !bin.exists() {
        eprintln!("Skipping CLI test: binary not found at {:?}", bin);
        return;
    }
    let out = std::process::Command::new(&bin)
        .arg("-V")
        .output()
        .expect("failed to run binary");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("RADDOSE-3D version"),
        "Expected version string in output, got: {stdout}"
    );
    assert!(out.status.success(), "-V should exit 0");
}

#[test]
fn cli_test_flag() {
    let bin = cli_bin();
    if !bin.exists() {
        return;
    }
    let out = std::process::Command::new(&bin)
        .arg("-t")
        .output()
        .expect("failed to run binary");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("Test run"),
        "Expected 'Test run' in stdout, got: {stdout}"
    );
    assert!(out.status.success(), "-t should exit 0");
}

#[test]
fn cli_no_args_exits_nonzero() {
    let bin = cli_bin();
    if !bin.exists() {
        return;
    }
    let out = std::process::Command::new(&bin)
        .output()
        .expect("failed to run binary");
    assert!(
        !out.status.success(),
        "Expected non-zero exit with no arguments"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("No input file"),
        "Expected 'No input file' in stderr, got: {stderr}"
    );
}

#[test]
fn cli_stdin_produces_output_files() {
    let bin = cli_bin();
    if !bin.exists() {
        return;
    }
    let tmp = tempdir();
    let prefix = format!("{}/rd3d-stdin-", tmp.path().display());
    let input = insulin_input();

    let out = std::process::Command::new(&bin)
        .args(["-i", "-", "-p", &prefix])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("failed to spawn")
        .wait_with_output_from_stdin(input.as_bytes());

    match out {
        Ok(o) => {
            assert!(o.status.success(), "CLI exited non-zero: {:?}", o.status);
            for suffix in &[
                "Summary.csv",
                "Summary.txt",
                "DoseState.csv",
                "DoseState.R",
                "RDE.csv",
                "DWDs.csv",
                "VoxelDose.csv",
                "VoxelFluences.csv",
            ] {
                let path = format!("{}{}", prefix, suffix);
                assert!(
                    std::path::Path::new(&path).exists(),
                    "Expected output file: {}",
                    path
                );
            }
        }
        Err(e) => panic!("CLI stdin run failed: {}", e),
    }
}

#[test]
fn cli_prefix_output_files_use_prefix() {
    let bin = cli_bin();
    if !bin.exists() {
        return;
    }
    let tmp = tempdir();
    let prefix = format!("{}/myprefix-", tmp.path().display());
    let fixture = format!("{}/insulin_test.txt", FIXTURES_DIR);

    let out = std::process::Command::new(&bin)
        .args(["-i", &fixture, "-p", &prefix])
        .output()
        .expect("failed to run binary");

    assert!(out.status.success(), "CLI exited non-zero");

    // Check that files with the given prefix exist.
    let summary = format!("{}Summary.txt", prefix);
    assert!(
        std::path::Path::new(&summary).exists(),
        "Expected {}",
        summary
    );
    let csv = format!("{}Summary.csv", prefix);
    assert!(std::path::Path::new(&csv).exists(), "Expected {}", csv);
}

#[test]
fn cli_summary_txt_avg_dwd_matches_golden() {
    let bin = cli_bin();
    if !bin.exists() {
        return;
    }
    let tmp = tempdir();
    let prefix = format!("{}/val-", tmp.path().display());
    let fixture = format!("{}/insulin_test.txt", FIXTURES_DIR);

    let out = std::process::Command::new(&bin)
        .args(["-i", &fixture, "-p", &prefix])
        .output()
        .expect("failed to run binary");
    assert!(out.status.success(), "CLI exited non-zero");

    let summary_path = format!("{}Summary.txt", prefix);
    let rust_summary = std::fs::read_to_string(&summary_path).expect("Summary.txt not written");
    let java_summary = golden("Summary.txt");

    let rust_dwd = extract_metric(&rust_summary, "Average Diffraction Weighted Dose");
    let java_dwd = extract_metric(&java_summary, "Average Diffraction Weighted Dose");

    assert!(
        (rust_dwd - java_dwd).abs() / java_dwd < 0.001,
        "CLI Avg DWD {:.6} differs from Java golden {:.6} by > 0.1%",
        rust_dwd,
        java_dwd
    );
}

// ── WASM crate logic (native Rust level) ─────────────────────────────────────

/// Exercise the WASM crate entry-points via the underlying Rust functions.
/// (wasm-bindgen machinery is skipped; we call the inner logic directly.)
#[test]
fn wasm_run_produces_summary_text() {
    // Replicate what raddose3d_wasm::run() does, but without wasm-bindgen.
    use raddose3d::{
        beam, crystal, experiment::Experiment, output::OutputSummaryText, wedge::Wedge,
        writer::shared_string_writer,
    };

    let input = insulin_input();
    let config = parse_input(&input).expect("parse failed");

    let (writer, buf_arc) = shared_string_writer();
    let mut exp = Experiment::new();
    exp.add_observer(Box::new(OutputSummaryText::new(writer)));

    for cc in &config.crystals {
        exp.set_crystal(crystal::create_crystal(cc).unwrap());
    }
    for bc in &config.beams {
        exp.set_beam(beam::create_beam(bc).unwrap());
    }
    for wc in &config.wedges {
        exp.expose_wedge(&Wedge::from_config(wc));
    }
    exp.close();

    let bytes = buf_arc.lock().unwrap().clone();
    let text = String::from_utf8(bytes).unwrap();

    assert!(
        text.contains("Average Diffraction Weighted Dose"),
        "WASM run() should produce summary text"
    );
    assert!(
        text.contains("Max Dose"),
        "WASM run() output missing Max Dose"
    );
}

#[test]
fn wasm_run_structured_matches_run_results() {
    let input = insulin_input();
    let config = parse_input(&input).expect("parse failed");
    let results = run(&config).expect("run failed");

    // runStructured() delegates to raddose3d::run() — verify the underlying
    // values are consistent.
    assert!(results.average_dwd > 0.0);
    assert!(results.max_dose > results.average_dwd);
    let results2 = run(&config).unwrap();
    let rel_diff = (results.average_dwd - results2.average_dwd).abs() / results.average_dwd;
    assert!(
        rel_diff < 1e-10,
        "run() must be deterministic: {:.16} vs {:.16} (rel diff {:.2e})",
        results.average_dwd,
        results2.average_dwd,
        rel_diff,
    );
}

#[test]
fn wasm_run_from_json_matches_text_input() {
    let text = insulin_input();
    let config_from_text = parse_input(&text).expect("parse_input failed");
    let json = serde_json::to_string(&config_from_text).unwrap();
    let config_from_json = parse_input_json(&json).expect("parse_input_json failed");

    let r1 = run(&config_from_text).expect("run from text failed");
    let r2 = run(&config_from_json).expect("run from json failed");

    // serde_json serialisation does not guarantee bit-exact f64 round-trips;
    // allow up to 1 ULP of difference (< 1e-12 relative for these values).
    let tol = 1e-12;
    assert!(
        (r1.average_dwd - r2.average_dwd).abs() < tol,
        "JSON round-trip average_dwd mismatch: {} vs {}",
        r1.average_dwd,
        r2.average_dwd
    );
    assert!(
        (r1.max_dose - r2.max_dose).abs() < tol,
        "JSON round-trip max_dose mismatch: {} vs {}",
        r1.max_dose,
        r2.max_dose
    );
}

// ── Temp-directory helper ─────────────────────────────────────────────────────

/// Very small temp-dir helper so we don't need a `tempfile` crate dependency.
struct TempDir(std::path::PathBuf);

impl TempDir {
    fn path(&self) -> &std::path::Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn tempdir() -> TempDir {
    use std::time::{SystemTime, UNIX_EPOCH};
    let ns = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .subsec_nanos();
    let path = std::env::temp_dir().join(format!("rd3d_test_{}", ns));
    std::fs::create_dir_all(&path).unwrap();
    TempDir(path)
}

// ── Stdin helper extension ────────────────────────────────────────────────────

trait WriteStdinExt {
    fn wait_with_output_from_stdin(self, data: &[u8]) -> std::io::Result<std::process::Output>;
}

impl WriteStdinExt for std::process::Child {
    fn wait_with_output_from_stdin(mut self, data: &[u8]) -> std::io::Result<std::process::Output> {
        use std::io::Write;
        if let Some(ref mut stdin) = self.stdin.take() {
            stdin.write_all(data)?;
        }
        self.wait_with_output()
    }
}
