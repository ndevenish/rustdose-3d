/// Build script: pre-parse ELSEPA and transition CSV data into compact binary blobs.
///
/// Binary formats
/// ─────────────
/// ELSEPA (`elsepa_{above|below}_{N}.bin`):
///   [n_rows: u32 LE][n_cols: u32 LE]
///   [n_rows × u64 LE]            energy keys (f64::to_bits)
///   [n_rows × n_cols × f32 LE]   probability values, row-major
///
/// TransitionData (`{auger|fl}_{stem}.bin`):
///   [n: u32 LE]
///   [n × f32 LE] × 6             linewidths, probs, energies,
///                                 cumulative_probs, exit_index, drop_index
use std::{
    env, fs,
    io::Write as _,
    path::{Path, PathBuf},
};

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let data = manifest.join("src/constants_data");
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());

    println!("cargo:rerun-if-changed=src/constants_data");

    process_elsepa(&data.join("above_20000"), &out, "above", 52);
    process_elsepa(&data.join("below_20000"), &out, "below", 94);

    let auger_stems = process_transitions(&data.join("auger_linewidths"), &out, "auger", true);
    let fl_stems = process_transitions(&data.join("fl_linewidths"), &out, "fl", false);

    generate_elsepa_rs(&out, 52, 94);
    generate_transition_rs(&out, &auger_stems, &fl_stems);
}

// ── ELSEPA ─────────────────────────────────────────────────────────────────────

fn process_elsepa(dir: &Path, out: &Path, prefix: &str, max_n: usize) {
    for n in 1..=max_n {
        let csv = dir.join(format!("{n}.csv"));
        let content = fs::read_to_string(&csv)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", csv.display()));

        let mut keys: Vec<u64> = Vec::new();
        let mut probs: Vec<f32> = Vec::new();
        let mut n_cols: Option<usize> = None;
        let mut first = true;

        for line in content.lines() {
            if first {
                first = false;
                continue; // skip header
            }
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 2 {
                continue;
            }
            let energy: f64 = match parts[0].trim().parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            keys.push(energy.to_bits());
            let row: Vec<f32> = parts[1..]
                .iter()
                .filter_map(|s| s.trim().parse::<f32>().ok())
                .collect();
            match n_cols {
                None => n_cols = Some(row.len()),
                Some(nc) => assert_eq!(nc, row.len(), "col mismatch in {}", csv.display()),
            }
            probs.extend_from_slice(&row);
        }

        let nc = n_cols.unwrap_or(0);
        let nr = keys.len();
        let out_path = out.join(format!("elsepa_{prefix}_{n}.bin"));
        let mut f = fs::File::create(&out_path).unwrap();
        f.write_all(&(nr as u32).to_le_bytes()).unwrap();
        f.write_all(&(nc as u32).to_le_bytes()).unwrap();
        for k in &keys {
            f.write_all(&k.to_le_bytes()).unwrap();
        }
        for p in &probs {
            f.write_all(&p.to_le_bytes()).unwrap();
        }
    }
}

fn generate_elsepa_rs(out: &Path, above_max: usize, below_max: usize) {
    let out_str = out.to_string_lossy().replace('\\', "/");
    let rs_path = out.join("elsepa_arrays.rs");
    let mut f = fs::File::create(&rs_path).unwrap();

    // ABOVE
    writeln!(
        f,
        "pub(super) static ELSEPA_ABOVE: [&[u8]; {}] = [",
        above_max + 1
    )
    .unwrap();
    writeln!(f, "    b\"\",  // index 0 unused").unwrap();
    for n in 1..=above_max {
        writeln!(f, "    include_bytes!(\"{out_str}/elsepa_above_{n}.bin\"),").unwrap();
    }
    writeln!(f, "];").unwrap();
    writeln!(f).unwrap();

    // BELOW
    writeln!(
        f,
        "pub(super) static ELSEPA_BELOW: [&[u8]; {}] = [",
        below_max + 1
    )
    .unwrap();
    writeln!(f, "    b\"\",  // index 0 unused").unwrap();
    for n in 1..=below_max {
        writeln!(f, "    include_bytes!(\"{out_str}/elsepa_below_{n}.bin\"),").unwrap();
    }
    writeln!(f, "];").unwrap();
}

// ── Transition data ─────────────────────────────────────────────────────────────

/// Returns sorted list of stems (filename without .csv) that were processed.
fn process_transitions(dir: &Path, out: &Path, prefix: &str, is_auger: bool) -> Vec<String> {
    let mut entries: Vec<_> = fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("cannot read dir {}: {e}", dir.display()))
        .map(|e| e.unwrap())
        .collect();
    entries.sort_by_key(|e| e.file_name());

    let mut stems = Vec::new();

    for entry in entries {
        let fname = entry.file_name();
        let fname = fname.to_string_lossy();
        if !fname.ends_with(".csv") {
            continue;
        }
        let stem = fname.trim_end_matches(".csv").to_string();

        let content = fs::read_to_string(entry.path()).unwrap();
        let mut linewidths: Vec<f32> = Vec::new();
        let mut probs_raw: Vec<f32> = Vec::new();
        let mut energies: Vec<f32> = Vec::new();
        let mut exit_idx: Vec<f32> = Vec::new();
        let mut drop_idx: Vec<f32> = Vec::new();
        let mut sum_prob = 0.0f64;

        for line in content.lines() {
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 4 {
                continue;
            }
            let lw: f64 = match parts[1].trim().parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            let prob: f64 = match parts[2].trim().parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            let energy: f64 = match parts[3].trim().parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            linewidths.push(lw as f32);
            probs_raw.push(prob as f32);
            energies.push(energy as f32);
            sum_prob += prob;

            if is_auger && parts.len() >= 6 {
                let exit: f64 = parts[4].trim().parse().unwrap_or(0.0);
                let drop: f64 = parts[5].trim().parse().unwrap_or(0.0);
                exit_idx.push(exit as f32);
                drop_idx.push(drop as f32);
            } else if !is_auger && parts.len() >= 5 {
                let drop: f64 = parts[4].trim().parse().unwrap_or(0.0);
                exit_idx.push(0.0);
                drop_idx.push(drop as f32);
            } else {
                exit_idx.push(0.0);
                drop_idx.push(0.0);
            }
        }

        if linewidths.is_empty() {
            continue;
        }

        // Compute normalised cumulative probs (same logic as parse_transition_csv)
        let mut cum_probs: Vec<f32> = Vec::with_capacity(probs_raw.len());
        let mut cum = 0.0f64;
        for &p in &probs_raw {
            cum += p as f64;
            cum_probs.push(if sum_prob > 0.0 {
                (cum / sum_prob) as f32
            } else {
                0.0
            });
        }

        let n = linewidths.len();
        let out_path = out.join(format!("{prefix}_{stem}.bin"));
        let mut file = fs::File::create(&out_path).unwrap();
        file.write_all(&(n as u32).to_le_bytes()).unwrap();
        for slice in [
            &linewidths,
            &probs_raw,
            &energies,
            &cum_probs,
            &exit_idx,
            &drop_idx,
        ] {
            for v in slice {
                file.write_all(&v.to_le_bytes()).unwrap();
            }
        }

        stems.push(stem);
    }

    stems
}

fn generate_transition_rs(out: &Path, auger_stems: &[String], fl_stems: &[String]) {
    let out_str = out.to_string_lossy().replace('\\', "/");
    let rs_path = out.join("transition_data.rs");
    let mut f = fs::File::create(&rs_path).unwrap();

    writeln!(
        f,
        "pub(super) fn get_transition_bin(set: &str, stem: &str) -> Option<&'static [u8]> {{"
    )
    .unwrap();
    writeln!(f, "    match (set, stem) {{").unwrap();
    for stem in auger_stems {
        writeln!(
            f,
            "        (\"auger\", \"{stem}\") => Some(include_bytes!(\"{out_str}/auger_{stem}.bin\")),"
        )
        .unwrap();
    }
    for stem in fl_stems {
        writeln!(
            f,
            "        (\"fl\", \"{stem}\") => Some(include_bytes!(\"{out_str}/fl_{stem}.bin\")),"
        )
        .unwrap();
    }
    writeln!(f, "        _ => None,").unwrap();
    writeln!(f, "    }}").unwrap();
    writeln!(f, "}}").unwrap();
}
