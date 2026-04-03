/// Size comparison between voxel grid and mesh blob representations.
use std::io::Write;

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::blob::Surface;
use crate::csv::Grid;

fn gzip_size(data: &[u8]) -> usize {
    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    gz.write_all(data).ok();
    gz.finish().map(|v| v.len()).unwrap_or(0)
}

pub fn print_comparison(grid: &Grid, surfaces: &[Surface]) {
    // Voxel grid: sparse (occupied only) — 12-byte header + 4 bytes per cell
    let n = grid.nx * grid.ny * grid.nz;
    let occupied: usize = grid.data.iter().filter(|&&v| v > 0.0).count();

    // Header: nx u32, ny u32, nz u32 (12 bytes) + dose values
    let voxel_raw = 12 + occupied * 4;

    let mut voxel_bytes = Vec::with_capacity(voxel_raw);
    voxel_bytes.extend_from_slice(&(grid.nx as u32).to_le_bytes());
    voxel_bytes.extend_from_slice(&(grid.ny as u32).to_le_bytes());
    voxel_bytes.extend_from_slice(&(grid.nz as u32).to_le_bytes());
    // Write occupied doses in grid order
    for v in &grid.data {
        if *v > 0.0 {
            voxel_bytes.extend_from_slice(&v.to_le_bytes());
        }
    }
    let voxel_gz = gzip_size(&voxel_bytes);

    // Mesh blob (raw, without gzip wrapper — just the per-surface data)
    let mesh_raw: usize = surfaces
        .iter()
        .map(|s| s.mesh.verts.len() * 12 + s.mesh.tris.len() * 12)
        .sum::<usize>()
        + 8   // magic + version + num_surfaces
        + 28; // bounds (24) + max_dose (4)

    // Approximate gzip of mesh payload
    let mut mesh_bytes: Vec<u8> = Vec::with_capacity(mesh_raw);
    mesh_bytes.extend_from_slice(b"RDISO\0");
    mesh_bytes.push(1u8);
    mesh_bytes.push(surfaces.len() as u8);
    for s in surfaces {
        for v in &s.mesh.verts {
            for c in v {
                mesh_bytes.extend_from_slice(&c.to_le_bytes());
            }
        }
        for t in &s.mesh.tris {
            for idx in t {
                mesh_bytes.extend_from_slice(&idx.to_le_bytes());
            }
        }
    }
    let mesh_gz = gzip_size(&mesh_bytes);

    let total_verts: usize = surfaces.iter().map(|s| s.mesh.verts.len()).sum();
    let total_tris: usize = surfaces.iter().map(|s| s.mesh.tris.len()).sum();

    println!();
    println!("┌────────────────────────────────┬──────────┬──────────┐");
    println!("│ Format                         │ Raw      │ Gzip     │");
    println!("├────────────────────────────────┼──────────┼──────────┤");
    println!(
        "│ Voxel grid (sparse, {}×{}×{}) │ {:>7}B │ {:>7}B │",
        grid.nx,
        grid.ny,
        grid.nz,
        format_bytes(voxel_raw),
        format_bytes(voxel_gz),
    );
    println!(
        "│ Mesh blob ({} surf, {} tri) │ {:>7}B │ {:>7}B │",
        surfaces.len(),
        total_tris,
        format_bytes(mesh_raw),
        format_bytes(mesh_gz),
    );
    println!("└────────────────────────────────┴──────────┴──────────┘");
    println!(
        "  Voxels: {} total, {} occupied ({:.1}%)",
        n,
        occupied,
        100.0 * occupied as f64 / n as f64
    );
    println!(
        "  Mesh: {} vertices, {} triangles across {} surfaces",
        total_verts,
        total_tris,
        surfaces.len()
    );

    let ratio = mesh_gz as f64 / voxel_gz as f64;
    if ratio < 1.0 {
        println!(
            "  Mesh blob is {:.1}× smaller than voxel grid (gzip)",
            1.0 / ratio
        );
    } else {
        println!("  Mesh blob is {:.1}× larger than voxel grid (gzip)", ratio);
    }
    println!();

    let _ = n; // suppress unused warning
}

fn format_bytes(n: usize) -> String {
    if n < 1024 {
        format!("{}", n)
    } else if n < 1024 * 1024 {
        format!("{:.1}K", n as f64 / 1024.0)
    } else {
        format!("{:.1}M", n as f64 / (1024.0 * 1024.0))
    }
}
