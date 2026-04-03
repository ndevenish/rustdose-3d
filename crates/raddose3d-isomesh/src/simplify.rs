/// Mesh simplification via meshopt (wraps meshoptimizer C++).
///
/// Steps:
///   1. Vertex welding via generate_vertex_remap
///   2. QEM simplification
use crate::mc::Mesh;
use meshopt::{SimplifyOptions, VertexDataAdapter};

pub fn simplify(mesh: Mesh, target_tris: usize) -> Mesh {
    if mesh.tris.is_empty() {
        return mesh;
    }

    let indices: Vec<u32> = mesh.tris.iter().flat_map(|t| t.iter().cloned()).collect();

    // Step 1: vertex welding
    let (unique_count, remap) =
        meshopt::generate_vertex_remap::<[f32; 3]>(&mesh.verts, Some(&indices));

    let welded_verts: Vec<[f32; 3]> =
        meshopt::remap_vertex_buffer::<[f32; 3]>(&mesh.verts, unique_count, &remap);
    let welded_indices = meshopt::remap_index_buffer(Some(&indices), mesh.verts.len(), &remap);

    // Step 2: QEM simplification
    // Build a VertexDataAdapter from the raw byte view of welded_verts
    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(
            welded_verts.as_ptr() as *const u8,
            welded_verts.len() * std::mem::size_of::<[f32; 3]>(),
        )
    };
    let adapter = VertexDataAdapter::new(byte_slice, std::mem::size_of::<[f32; 3]>(), 0)
        .expect("valid vertex adapter");

    let target_index_count = (target_tris * 3).min(welded_indices.len());
    let simplified_indices = meshopt::simplify(
        &welded_indices,
        &adapter,
        target_index_count,
        0.01,
        SimplifyOptions::None,
        None,
    );

    let out_tris = simplified_indices
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();

    Mesh {
        verts: welded_verts,
        tris: out_tris,
    }
}
