mod blob;
mod compare;
mod csv;
mod iso;
mod mc;
mod simplify;

use std::fs;
use std::path::PathBuf;

use clap::Parser;

use blob::{write_blob, Surface};
use csv::parse_voxel_dose_csv;
use iso::compute_iso_levels;
use mc::marching_cubes;
use simplify::simplify;

#[derive(Parser, Debug)]
#[command(
    name = "raddose3d-isomesh",
    about = "Convert a VoxelDose.csv into a compact isosurface blob (.rdiso)"
)]
struct Args {
    /// Input VoxelDose CSV file
    voxel_dose_csv: PathBuf,

    /// Output .rdiso file
    output: PathBuf,

    /// Target triangles per surface after simplification
    #[arg(long, default_value_t = 2000)]
    target_tris: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // --- Parse CSV ---
    eprintln!("Reading {}...", args.voxel_dose_csv.display());
    let content = fs::read_to_string(&args.voxel_dose_csv)
        .map_err(|e| format!("cannot read {}: {e}", args.voxel_dose_csv.display()))?;

    let grid = parse_voxel_dose_csv(&content).map_err(|e| format!("CSV parse error: {e}"))?;

    let max_dose = grid.max_dose();
    eprintln!(
        "Grid: {}×{}×{}, spacing={:.3} µm, max_dose={:.4} MGy",
        grid.nx, grid.ny, grid.nz, grid.spacing, max_dose
    );

    // --- Iso levels ---
    let levels = compute_iso_levels(max_dose);
    if levels.is_empty() {
        eprintln!("Max dose is zero or negative — nothing to render.");
        return Ok(());
    }
    eprintln!(
        "Iso levels: {}",
        levels
            .iter()
            .map(|l| format!("{:.3}", l.value))
            .collect::<Vec<_>>()
            .join(", ")
    );

    // --- Bounds ---
    let min = grid.min;
    let max = [
        min[0] + (grid.nx as f32 - 1.0) * grid.spacing,
        min[1] + (grid.ny as f32 - 1.0) * grid.spacing,
        min[2] + (grid.nz as f32 - 1.0) * grid.spacing,
    ];

    // --- Marching cubes + simplification per surface ---
    let mut surfaces: Vec<Surface> = Vec::new();
    for level in levels {
        let iso_val = level.value;
        eprint!("  MC iso={:.4}...", iso_val);
        let raw = marching_cubes(&grid, iso_val);
        eprintln!(" {} verts, {} tris", raw.verts.len(), raw.tris.len());

        let simplified = if raw.tris.len() > args.target_tris {
            eprint!("  Simplifying to {} tris...", args.target_tris);
            let s = simplify(raw, args.target_tris);
            eprintln!(" → {} verts, {} tris", s.verts.len(), s.tris.len());
            s
        } else {
            raw
        };

        surfaces.push(Surface {
            level,
            mesh: simplified,
        });
    }

    // --- Size comparison ---
    compare::print_comparison(&grid, &surfaces);

    // --- Write blob ---
    eprintln!("Writing {}...", args.output.display());
    let out_file = fs::File::create(&args.output)
        .map_err(|e| format!("cannot create {}: {e}", args.output.display()))?;

    write_blob(out_file, &surfaces, max_dose, min, max)?;

    let meta = fs::metadata(&args.output)?;
    eprintln!("Done. Output: {} bytes", meta.len());

    Ok(())
}
