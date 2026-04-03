/// Binary blob writer for .rdiso files.
///
/// Format (little-endian, then gzip-compressed):
///   [6]   magic: "RDISO\0"
///   [1]   version: u8 = 1
///   [1]   num_surfaces: u8
///   [4*6] bounds: min_x,min_y,min_z,max_x,max_y,max_z (f32, µm)
///   [4]   max_dose: f32
///
///   Per surface:
///   [4]   iso_value: f32
///   [12]  color: r, g, b (f32 each, 0..1)
///   [4]   opacity: f32
///   [4]   num_verts: u32
///   [4]   num_tris: u32
///   [num_verts * 12]  verts: x,y,z interleaved (f32)
///   [num_tris * 12]   indices: i,j,k per tri (u32)
use std::io::{self, Write};

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::iso::IsoLevel;
use crate::mc::Mesh;

pub struct Surface {
    pub level: IsoLevel,
    pub mesh: Mesh,
}

pub fn write_blob<W: Write>(
    writer: W,
    surfaces: &[Surface],
    max_dose: f32,
    bounds_min: [f32; 3],
    bounds_max: [f32; 3],
) -> io::Result<()> {
    let mut gz = GzEncoder::new(writer, Compression::default());

    // Magic + version + num_surfaces
    gz.write_all(b"RDISO\0")?;
    gz.write_all(&[1u8])?; // version
    gz.write_all(&[surfaces.len() as u8])?;

    // Bounds
    for v in &bounds_min {
        gz.write_all(&v.to_le_bytes())?;
    }
    for v in &bounds_max {
        gz.write_all(&v.to_le_bytes())?;
    }

    // Max dose
    gz.write_all(&max_dose.to_le_bytes())?;

    // Surfaces
    for surf in surfaces {
        gz.write_all(&surf.level.value.to_le_bytes())?;
        for c in &surf.level.color {
            gz.write_all(&c.to_le_bytes())?;
        }
        gz.write_all(&surf.level.opacity.to_le_bytes())?;

        let nv = surf.mesh.verts.len() as u32;
        let nt = surf.mesh.tris.len() as u32;
        gz.write_all(&nv.to_le_bytes())?;
        gz.write_all(&nt.to_le_bytes())?;

        for v in &surf.mesh.verts {
            for c in v {
                gz.write_all(&c.to_le_bytes())?;
            }
        }
        for t in &surf.mesh.tris {
            for idx in t {
                gz.write_all(&idx.to_le_bytes())?;
            }
        }
    }

    gz.finish()?;
    Ok(())
}
