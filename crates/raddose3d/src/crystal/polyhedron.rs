use crate::beam::Beam;
use crate::coefcalc::{self, CoefCalc};
use crate::container::{self, Container, ContainerTransparent};
use crate::ddm::{self, DdmModel};
use crate::output::ExposureSummary;
use crate::parser::config::CrystalConfig;
use crate::wedge::Wedge;

/// General polyhedron crystal: arbitrary vertices and triangulated faces.
/// Used for OBJ-loaded crystals, cylinders (hardcoded geometry), and
/// icosphere approximations of spheres (SphericalNew).
#[derive(Debug)]
#[allow(dead_code)] // phase-6 fields (photo_electron_escape, fluorescent_escape)
pub struct CrystalPolyhedron {
    // Geometry (post P/L rotation, pre-wedge)
    vertices: Vec<[f64; 3]>,
    indices: Vec<[usize; 3]>,
    normals: Vec<[f64; 3]>,
    origin_distances: Vec<f64>,
    // Rotated geometry (updated per-angle in setup_depth_finding)
    rotated_vertices: Vec<[f64; 3]>,
    rotated_normals: Vec<[f64; 3]>,
    rotated_origin_distances: Vec<f64>,
    /// Pre-expanded rotated vertices for fast triangle lookup in find_depth.
    expanded_rotated_vertices: Vec<[[f64; 3]; 3]>,
    // Crystal size
    cryst_size_um: [f64; 3],
    cryst_size_voxels: [usize; 3],
    pix_per_um: f64,
    // Voxel data
    cryst_coord: Vec<[f64; 3]>,
    dose: Vec<f64>,
    fluence: Vec<f64>,
    elastic: Vec<f64>,
    // Components
    coefcalc: Box<dyn CoefCalc>,
    ddm: Box<dyn DdmModel>,
    container: Box<dyn Container>,
    exposure_summary: ExposureSummary,
    subprogram: String,
    photo_electron_escape: bool,
    fluorescent_escape: bool,
    calc_surrounding: bool,
    crystal_name: String,
    angle_p: f64,
    angle_l: f64,
}

impl CrystalPolyhedron {
    const DEFAULT_RESOLUTION: f64 = 0.5;

    /// Create from raw geometry (vertices already in final coordinate space,
    /// indices already 0-based) and a crystal config.
    pub fn from_geometry_and_config(
        mut vertices_raw: Vec<[f64; 3]>,
        indices_0based: Vec<[usize; 3]>,
        config: &CrystalConfig,
        name: &str,
    ) -> Result<Self, String> {
        // Apply P/L rotation to vertices (Java's sign convention)
        let angle_p = config.angle_p.unwrap_or(0.0).to_radians();
        let angle_l = config.angle_l.unwrap_or(0.0).to_radians();

        if angle_p != 0.0 || angle_l != 0.0 {
            let cos_p = angle_p.cos();
            let sin_p = angle_p.sin();
            let cos_l = angle_l.cos();
            let sin_l = angle_l.sin();
            for v in &mut vertices_raw {
                let x = v[0];
                let y = v[1];
                let z = v[2];
                // Rotate about Z (Java's convention: x2 = x*cos_p + y*sin_p)
                let x2 = x * cos_p + y * sin_p;
                let y2 = -x * sin_p + y * cos_p;
                // Rotate about X
                let y3 = y2 * cos_l + z * sin_l;
                let z3 = -y2 * sin_l + z * cos_l;
                v[0] = x2;
                v[1] = y3;
                v[2] = z3;
            }
        }

        // Bounding box of rotated vertices
        let mut min_coords = [f64::INFINITY; 3];
        let mut max_coords = [f64::NEG_INFINITY; 3];
        for v in &vertices_raw {
            for dim in 0..3 {
                min_coords[dim] = min_coords[dim].min(v[dim]);
                max_coords[dim] = max_coords[dim].max(v[dim]);
            }
        }

        let dim_x = max_coords[0] - min_coords[0];
        let dim_y = max_coords[1] - min_coords[1];
        let dim_z = max_coords[2] - min_coords[2];
        let max_dim = dim_x.max(dim_y).max(dim_z).max(1e-9);

        let pix_per_um = config
            .pixels_per_micron
            .unwrap_or_else(|| (10.0 / max_dim).min(Self::DEFAULT_RESOLUTION));

        let nx = (dim_x * pix_per_um).round() as usize + 1;
        let ny = (dim_y * pix_per_um).round() as usize + 1;
        let nz = (dim_z * pix_per_um).round() as usize + 1;
        let total_voxels = nx * ny * nz;

        // Voxel coordinates start at bounding box minimum
        let mut cryst_coord = vec![[0.0f64; 3]; total_voxels];
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = i * ny * nz + j * nz + k;
                    cryst_coord[idx] = [
                        min_coords[0] + i as f64 / pix_per_um,
                        min_coords[1] + j as f64 / pix_per_um,
                        min_coords[2] + k as f64 / pix_per_um,
                    ];
                }
            }
        }

        // Compute normals and origin distances
        let num_tris = indices_0based.len();
        let mut normals = vec![[0.0f64; 3]; num_tris];
        let mut origin_distances = vec![0.0f64; num_tris];
        for (t, tri) in indices_0based.iter().enumerate() {
            let (n, d) = calc_triangle_normal_and_dist(&vertices_raw, tri);
            normals[t] = n;
            origin_distances[t] = d;
        }

        // Build expanded rotated vertices (initially same as vertices)
        let expanded_rotated_vertices = build_expanded_vertices(&vertices_raw, &indices_0based);

        let coefcalc = coefcalc::create_coefcalc(config)?;
        let ddm = ddm::create_ddm(
            config.ddm,
            config.gamma_param,
            config.b0_param,
            config.beta_param,
        );

        let subprogram = config.program.as_deref().unwrap_or("RD3D").to_uppercase();
        match subprogram.as_str() {
            "RD3D" | "" | "MONTECARLO" | "GOS" | "XFEL" | "EMSP" | "EMED" => {}
            other => {
                return Err(format!("Unrecognised subprogram '{other}'"));
            }
        }

        let photo_electron_escape = config
            .calculate_pe_escape
            .as_deref()
            .is_some_and(|s| s.eq_ignore_ascii_case("true"));
        let fluorescent_escape = config
            .calculate_fl_escape
            .as_deref()
            .is_some_and(|s| s.eq_ignore_ascii_case("true"));

        let calc_surrounding = config
            .calc_surrounding
            .as_deref()
            .is_some_and(|s| s.eq_ignore_ascii_case("true"));

        let container = container::create_container(config);

        let rotated_normals = normals.clone();
        let rotated_origin_distances = origin_distances.clone();
        let rotated_vertices = vertices_raw.clone();

        Ok(CrystalPolyhedron {
            vertices: vertices_raw,
            indices: indices_0based,
            normals,
            origin_distances,
            rotated_vertices,
            rotated_normals,
            rotated_origin_distances,
            expanded_rotated_vertices,
            cryst_size_um: [dim_x, dim_y, dim_z],
            cryst_size_voxels: [nx, ny, nz],
            pix_per_um,
            cryst_coord,
            dose: vec![0.0; total_voxels],
            fluence: vec![0.0; total_voxels],
            elastic: vec![0.0; total_voxels],
            coefcalc,
            ddm,
            container,
            exposure_summary: ExposureSummary::new(),
            subprogram,
            photo_electron_escape,
            fluorescent_escape,
            calc_surrounding,
            crystal_name: name.to_string(),
            angle_p,
            angle_l,
        })
    }

    /// Create from config for OBJ-loaded polyhedron ("polyhedron" / "obj" type).
    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        // Dimensions must be specified explicitly (matches Java behaviour).
        if config.dim_x.is_none() {
            return Err("Polyhedron crystal requires Dimensions to be specified".to_string());
        }

        let model_file = config
            .model_file
            .as_deref()
            .ok_or("Polyhedron crystal requires ModelFile")?;

        let (vertices, indices) = load_obj(model_file)?;

        // OBJ vertices are in absolute µm coordinates (matching Java behaviour).
        // Dimensions from config define the bounding box; no auto-derivation
        // from vertex coordinates (Java does not do this either).
        Self::from_geometry_and_config(vertices, indices, config, "Polyhedron")
    }

    #[inline]
    fn voxel_idx(&self, i: usize, j: usize, k: usize) -> usize {
        let [_, ny, nz] = self.cryst_size_voxels;
        i * ny * nz + j * nz + k
    }

    /// Check if a point is inside the polyhedron using ray-casting (+Z direction).
    /// Matches Java's `calculateCrystalOccupancy` exactly.
    fn is_inside_polyhedron(&self, coord: &[f64; 3]) -> bool {
        let mut crossings = 0usize;

        for (t, tri) in self.indices.iter().enumerate() {
            let n = &self.normals[t];
            let d = self.origin_distances[t];

            let denom = n[2]; // n · [0,0,1]

            // Ray-plane: t = -(n·origin + d) / (n·dir), dir=[0,0,1], n·dir = n[2]
            let t_val = -(n[0] * coord[0] + n[1] * coord[1] + n[2] * coord[2] + d) / denom;
            if t_val < 0.0 || t_val.is_nan() || t_val.is_infinite() {
                continue;
            }

            let p = [coord[0], coord[1], coord[2] + t_val];

            if point_in_triangle(
                &p,
                &self.vertices[tri[0]],
                &self.vertices[tri[1]],
                &self.vertices[tri[2]],
            ) {
                crossings += 1;
            }
        }

        crossings % 2 == 1
    }
}

/// OBJ geometry: vertices and triangle indices.
pub type ObjGeometry = (Vec<[f64; 3]>, Vec<[usize; 3]>);

/// Load an OBJ file and return (vertices, 0-based indices).
pub fn load_obj(path: &str) -> Result<ObjGeometry, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read OBJ file '{}': {}", path, e))?;

    let mut vertices = Vec::new();
    let mut indices = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with("v ") {
            let parts: Vec<&str> = line
                .strip_prefix("v ")
                .unwrap_or("")
                .split_whitespace()
                .collect();
            if parts.len() < 3 {
                return Err(format!("Invalid vertex line in OBJ: '{}'", line));
            }
            let x: f64 = parts[0]
                .parse()
                .map_err(|_| format!("Bad float: {}", parts[0]))?;
            let y: f64 = parts[1]
                .parse()
                .map_err(|_| format!("Bad float: {}", parts[1]))?;
            let z: f64 = parts[2]
                .parse()
                .map_err(|_| format!("Bad float: {}", parts[2]))?;
            vertices.push([x, y, z]);
        } else if line.starts_with("f ") {
            let parts: Vec<&str> = line
                .strip_prefix("f ")
                .unwrap_or("")
                .split_whitespace()
                .collect();
            if parts.len() < 3 {
                return Err(format!("Invalid face line in OBJ: '{}'", line));
            }
            // Parse index: "a", "a/b", or "a//b" — take the first integer (1-based)
            let parse_idx = |s: &str| -> Result<usize, String> {
                let first = s.split('/').next().unwrap_or(s);
                let idx: usize = first
                    .parse()
                    .map_err(|_| format!("Bad face index: {}", s))?;
                Ok(idx - 1) // Convert to 0-based
            };
            // Triangulate face (fan triangulation for polygons)
            let first = parse_idx(parts[0])?;
            for i in 1..(parts.len() - 1) {
                let a = first;
                let b = parse_idx(parts[i])?;
                let c = parse_idx(parts[i + 1])?;
                indices.push([a, b, c]);
            }
        }
    }

    Ok((vertices, indices))
}

/// Build cylinder geometry matching Java's CrystalCylinder.
/// Returns (vertices, 0-based indices).
pub fn cylinder_geometry(radius: f64, height: f64) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let num_vertices: usize = 32;
    let mut vertices = Vec::with_capacity(num_vertices * 2);

    // Build interleaved base/top vertices with Java's angle formula,
    // then apply 90° rotation about Z using the actual rotation matrix
    // (matching Java's Math.cos/sin with floating-point imprecision).
    let angle_rad = 90.0_f64.to_radians();
    let cos_a = angle_rad.cos(); // ~6.12e-17, not exactly 0
    let sin_a = angle_rad.sin(); // 1.0
                                 // R_z = [[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]]

    for vertex in 0..num_vertices {
        let angle = -(2.0 * std::f64::consts::PI / num_vertices as f64) * vertex as f64;
        let y_coord = radius * angle.cos();
        let z_coord = radius * angle.sin();

        // Base vertex pre-rotation: (-height/2, y_coord, z_coord)
        let bx = -height / 2.0;
        vertices.push([
            cos_a * bx - sin_a * y_coord,
            sin_a * bx + cos_a * y_coord,
            z_coord,
        ]);
        // Top vertex pre-rotation: (+height/2, y_coord, z_coord)
        let tx = height / 2.0;
        vertices.push([
            cos_a * tx - sin_a * y_coord,
            sin_a * tx + cos_a * y_coord,
            z_coord,
        ]);
    }

    // Java indices are 1-based with 96 triangles, convert to 0-based
    let java_indices: &[[usize; 3]] = &[
        [2, 4, 3],
        [3, 4, 6],
        [6, 8, 7],
        [8, 10, 9],
        [10, 12, 11],
        [12, 14, 13],
        [14, 16, 15],
        [15, 16, 18],
        [18, 20, 19],
        [20, 22, 21],
        [21, 22, 24],
        [24, 26, 25],
        [26, 28, 27],
        [27, 28, 30],
        [30, 32, 31],
        [32, 34, 33],
        [34, 36, 35],
        [35, 36, 38],
        [38, 40, 39],
        [39, 40, 42],
        [42, 44, 43],
        [44, 46, 45],
        [45, 46, 48],
        [48, 50, 49],
        [49, 50, 52],
        [52, 54, 53],
        [53, 54, 56],
        [56, 58, 57],
        [58, 60, 59],
        [59, 60, 62],
        [62, 38, 22],
        [64, 2, 1],
        [62, 64, 63],
        [53, 55, 63],
        [1, 2, 3],
        [5, 3, 6],
        [5, 6, 7],
        [7, 8, 9],
        [9, 10, 11],
        [11, 12, 13],
        [13, 14, 15],
        [17, 15, 18],
        [17, 18, 19],
        [19, 20, 21],
        [23, 21, 24],
        [23, 24, 25],
        [25, 26, 27],
        [29, 27, 30],
        [29, 30, 31],
        [31, 32, 33],
        [33, 34, 35],
        [37, 35, 38],
        [37, 38, 39],
        [41, 39, 42],
        [41, 42, 43],
        [43, 44, 45],
        [47, 45, 48],
        [47, 48, 49],
        [51, 49, 52],
        [51, 52, 53],
        [55, 53, 56],
        [55, 56, 57],
        [57, 58, 59],
        [61, 59, 62],
        [2, 64, 4],
        [64, 62, 14],
        [58, 54, 60],
        [58, 56, 54],
        [54, 52, 50],
        [48, 46, 44],
        [42, 48, 44],
        [42, 40, 38],
        [34, 32, 36],
        [30, 22, 32],
        [26, 22, 28],
        [26, 24, 22],
        [22, 20, 18],
        [14, 22, 16],
        [10, 8, 12],
        [54, 62, 60],
        [64, 6, 4],
        [8, 14, 12],
        [22, 30, 28],
        [38, 48, 42],
        [38, 36, 22],
        [6, 64, 8],
        [61, 62, 63],
        [36, 32, 22],
        [54, 38, 62],
        [38, 50, 48],
        [38, 54, 50],
        [8, 64, 14],
        [63, 64, 1],
        [22, 18, 16],
        [62, 22, 14],
        [63, 1, 3],
        [3, 5, 7],
        [15, 9, 13],
        [15, 17, 19],
        [21, 23, 31],
        [23, 25, 27],
        [27, 29, 31],
        [31, 33, 35],
        [35, 37, 39],
        [39, 41, 43],
        [47, 51, 45],
        [47, 49, 51],
        [53, 39, 43],
        [59, 63, 57],
        [59, 61, 63],
        [7, 31, 3],
        [9, 11, 13],
        [63, 55, 57],
        [23, 27, 31],
        [39, 63, 35],
        [45, 51, 53],
        [15, 7, 9],
        [15, 19, 21],
        [7, 21, 31],
        [15, 21, 7],
        [39, 53, 63],
        [31, 63, 3],
        [63, 31, 35],
        [43, 45, 53],
    ];

    // Convert 1-based to 0-based
    let indices: Vec<[usize; 3]> = java_indices
        .iter()
        .map(|tri| [tri[0] - 1, tri[1] - 1, tri[2] - 1])
        .collect();

    (vertices, indices)
}

/// Build icosphere geometry matching Java's CrystalSphericalNew (42 vertices, 80 triangles).
/// Vertices are scaled by `diameter`.
pub fn icosphere_geometry(diameter: f64) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let raw_vertices: &[[f64; 3]] = &[
        [0.000000, -0.500000, 0.000000],
        [0.361804, -0.223610, 0.262863],
        [-0.138194, -0.223610, 0.425325],
        [-0.447213, -0.223608, 0.000000],
        [-0.138194, -0.223610, -0.425325],
        [0.361804, -0.223610, -0.262863],
        [0.138194, 0.223610, 0.425325],
        [-0.361804, 0.223610, 0.262863],
        [-0.361804, 0.223610, -0.262863],
        [0.138194, 0.223610, -0.425325],
        [0.447213, 0.223608, 0.000000],
        [0.000000, 0.500000, 0.000000],
        [-0.081228, -0.425327, 0.249998],
        [0.212661, -0.425327, 0.154506],
        [0.131434, -0.262869, 0.404506],
        [0.425324, -0.262868, 0.000000],
        [0.212661, -0.425327, -0.154506],
        [-0.262865, -0.425326, 0.000000],
        [-0.344095, -0.262868, 0.249998],
        [-0.081228, -0.425327, -0.249998],
        [-0.344095, -0.262868, -0.249998],
        [0.131434, -0.262869, -0.404506],
        [0.475529, 0.000000, 0.154506],
        [0.475529, 0.000000, -0.154506],
        [0.000000, 0.000000, 0.500000],
        [0.293893, 0.000000, 0.404508],
        [-0.475529, 0.000000, 0.154506],
        [-0.293893, 0.000000, 0.404508],
        [-0.293893, 0.000000, -0.404508],
        [-0.475529, 0.000000, -0.154506],
        [0.293893, 0.000000, -0.404508],
        [0.000000, 0.000000, -0.500000],
        [0.344095, 0.262868, 0.249998],
        [-0.131434, 0.262869, 0.404506],
        [-0.425324, 0.262868, 0.000000],
        [-0.131434, 0.262869, -0.404506],
        [0.344095, 0.262868, -0.249998],
        [0.081228, 0.425327, 0.249998],
        [0.262865, 0.425326, 0.000000],
        [-0.212661, 0.425327, 0.154506],
        [-0.212661, 0.425327, -0.154506],
        [0.081228, 0.425327, -0.249998],
    ];

    let vertices: Vec<[f64; 3]> = raw_vertices
        .iter()
        .map(|v| [v[0] * diameter, v[1] * diameter, v[2] * diameter])
        .collect();

    // 1-based Java indices, converted to 0-based
    let java_indices: &[[usize; 3]] = &[
        [1, 14, 13],
        [2, 14, 16],
        [1, 13, 18],
        [1, 18, 20],
        [1, 20, 17],
        [2, 16, 23],
        [3, 15, 25],
        [4, 19, 27],
        [5, 21, 29],
        [6, 22, 31],
        [2, 23, 26],
        [3, 25, 28],
        [4, 27, 30],
        [5, 29, 32],
        [6, 31, 24],
        [7, 33, 38],
        [8, 34, 40],
        [9, 35, 41],
        [10, 36, 42],
        [11, 37, 39],
        [39, 42, 12],
        [39, 37, 42],
        [37, 10, 42],
        [42, 41, 12],
        [42, 36, 41],
        [36, 9, 41],
        [41, 40, 12],
        [41, 35, 40],
        [35, 8, 40],
        [40, 38, 12],
        [40, 34, 38],
        [34, 7, 38],
        [38, 39, 12],
        [38, 33, 39],
        [33, 11, 39],
        [24, 37, 11],
        [24, 31, 37],
        [31, 10, 37],
        [32, 36, 10],
        [32, 29, 36],
        [29, 9, 36],
        [30, 35, 9],
        [30, 27, 35],
        [27, 8, 35],
        [28, 34, 8],
        [28, 25, 34],
        [25, 7, 34],
        [26, 33, 7],
        [26, 23, 33],
        [23, 11, 33],
        [31, 32, 10],
        [31, 22, 32],
        [22, 5, 32],
        [29, 30, 9],
        [29, 21, 30],
        [21, 4, 30],
        [27, 28, 8],
        [27, 19, 28],
        [19, 3, 28],
        [25, 26, 7],
        [25, 15, 26],
        [15, 2, 26],
        [23, 24, 11],
        [23, 16, 24],
        [16, 6, 24],
        [17, 22, 6],
        [17, 20, 22],
        [20, 5, 22],
        [20, 21, 5],
        [20, 18, 21],
        [18, 4, 21],
        [18, 19, 4],
        [18, 13, 19],
        [13, 3, 19],
        [16, 17, 6],
        [16, 14, 17],
        [14, 1, 17],
        [13, 15, 3],
        [13, 14, 15],
        [14, 2, 15],
    ];

    let indices: Vec<[usize; 3]> = java_indices
        .iter()
        .map(|tri| [tri[0] - 1, tri[1] - 1, tri[2] - 1])
        .collect();

    (vertices, indices)
}

/// Create a cylinder crystal from config.
pub fn crystal_cylinder_from_config(config: &CrystalConfig) -> Result<CrystalPolyhedron, String> {
    let diameter = config
        .dim_x
        .ok_or("Cylinder crystal requires DimX (diameter)")?;
    let height = config.dim_y.unwrap_or(diameter);
    let radius = diameter / 2.0;

    let (vertices, indices) = cylinder_geometry(radius, height);
    CrystalPolyhedron::from_geometry_and_config(vertices, indices, config, "Cylinder")
}

/// Create an icosphere (SphericalNew) crystal from config.
pub fn crystal_spherical_new_from_config(
    config: &CrystalConfig,
) -> Result<CrystalPolyhedron, String> {
    let diameter = config
        .dim_x
        .ok_or("SphericalNew crystal requires DimX (diameter)")?;
    let (vertices, indices) = icosphere_geometry(diameter);
    CrystalPolyhedron::from_geometry_and_config(vertices, indices, config, "SphericalNew")
}

impl super::Crystal for CrystalPolyhedron {
    fn setup_depth_finding(&mut self, ang_rad: f64, wedge: &Wedge) {
        let cos_a = ang_rad.cos();
        let sin_a = ang_rad.sin();

        let delta = ang_rad - wedge.start_ang;
        let tx = wedge.start_x + wedge.trans_x * delta;
        let ty = wedge.start_y + wedge.trans_y * delta;
        let tz = wedge.start_z + wedge.trans_z * delta;

        for (i, v) in self.vertices.iter().enumerate() {
            let x = v[0] + tx;
            let y = v[1] + ty;
            let z = v[2] + tz;

            // Rotate in X-Z plane (about Y axis)
            self.rotated_vertices[i] = [x * cos_a + z * sin_a, y, -x * sin_a + z * cos_a];
        }

        // Recompute normals for rotated vertices
        for (t, tri) in self.indices.iter().enumerate() {
            let (n, d) = calc_triangle_normal_and_dist(&self.rotated_vertices, tri);
            self.rotated_normals[t] = n;
            self.rotated_origin_distances[t] = d;
        }

        // Rebuild expanded rotated vertices
        self.expanded_rotated_vertices =
            build_expanded_vertices(&self.rotated_vertices, &self.indices);
    }

    fn find_depth(&self, vox_coord: &[f64; 3], _delta_phi: f64, _wedge: &Wedge) -> f64 {
        let dir_z = -1.0f64; // Ray direction [0,0,-1] (towards beam source)
        let mut distances = Vec::new();

        for (t, _tri) in self.indices.iter().enumerate() {
            let n = &self.rotated_normals[t];
            let d = self.rotated_origin_distances[t];

            let denom = n[2] * dir_z; // n · [0,0,-1]
            if denom.abs() < 1e-12 {
                continue;
            }

            let t_val =
                -(n[0] * vox_coord[0] + n[1] * vox_coord[1] + n[2] * vox_coord[2] + d) / denom;
            if t_val <= 0.0 || t_val.is_nan() || t_val.is_infinite() {
                continue;
            }

            let p = [vox_coord[0], vox_coord[1], vox_coord[2] + t_val * dir_z];

            let exp_verts = &self.expanded_rotated_vertices[t];
            if point_in_triangle(&p, &exp_verts[0], &exp_verts[1], &exp_verts[2]) {
                distances.push(t_val);
            }
        }

        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
        distances.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

        if distances.is_empty() || distances.len() % 2 == 0 {
            return 0.0;
        }

        let mut depth = distances[0];
        let mut idx = 1;
        while idx + 1 < distances.len() {
            depth += distances[idx + 1] - distances[idx];
            idx += 2;
        }

        depth
    }

    fn get_cryst_coord(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        self.cryst_coord[self.voxel_idx(i, j, k)]
    }

    fn num_images(&self, wedge: &Wedge) -> usize {
        if (wedge.start_ang - wedge.end_ang).abs() < wedge.ang_res {
            crate::constants::STATIC_EXPOSURE
        } else {
            let sign: f64 = if wedge.end_ang < wedge.start_ang {
                -1.0
            } else {
                1.0
            };
            (sign * (wedge.end_ang - wedge.start_ang) / wedge.ang_res + 1.0) as usize
        }
    }

    fn is_crystal_at(&self, i: usize, j: usize, k: usize) -> bool {
        let [nx, ny, nz] = self.cryst_size_voxels;
        if i >= nx || j >= ny || k >= nz {
            return false;
        }
        let coord = self.cryst_coord[self.voxel_idx(i, j, k)];
        self.is_inside_polyhedron(&coord)
    }

    fn add_dose(&mut self, i: usize, j: usize, k: usize, dose_increase: f64) {
        let idx = self.voxel_idx(i, j, k);
        self.dose[idx] += dose_increase;
    }

    fn add_fluence(&mut self, i: usize, j: usize, k: usize, fluence_increase: f64) {
        let idx = self.voxel_idx(i, j, k);
        self.fluence[idx] += fluence_increase;
    }

    fn add_elastic(&mut self, i: usize, j: usize, k: usize, elastic_increase: f64) {
        let idx = self.voxel_idx(i, j, k);
        self.elastic[idx] += elastic_increase;
    }

    fn crystal_info(&self) -> String {
        let base = if self.crystal_name == "Cylinder" {
            format!(
                "Cylinder (Polyhedron) crystal of diameter {:.2} mm and height {:.2} mm at a resolution of {:.2} microns per voxel edge.",
                self.cryst_size_um[0] / 1000.0,
                self.cryst_size_um[1] / 1000.0,
                1.0 / self.pix_per_um
            )
        } else {
            format!(
                "{} crystal of size [{:.0}, {:.0}, {:.0}] um [x, y, z] at a resolution of {:.2} microns per voxel edge.",
                self.crystal_name,
                self.cryst_size_um[0],
                self.cryst_size_um[1],
                self.cryst_size_um[2],
                1.0 / self.pix_per_um
            )
        };
        if self.angle_p != 0.0 || self.angle_l != 0.0 {
            format!(
                "{} P={:.2} deg, L={:.2} deg.",
                base,
                self.angle_p.to_degrees(),
                self.angle_l.to_degrees()
            )
        } else {
            base
        }
    }

    fn cryst_size_voxels(&self) -> [usize; 3] {
        self.cryst_size_voxels
    }

    fn cryst_size_um(&self) -> [f64; 3] {
        self.cryst_size_um
    }

    fn get_dose(&self, i: usize, j: usize, k: usize) -> f64 {
        self.dose[self.voxel_idx(i, j, k)]
    }

    fn get_fluence(&self, i: usize, j: usize, k: usize) -> f64 {
        self.fluence[self.voxel_idx(i, j, k)]
    }

    fn get_elastic(&self, i: usize, j: usize, k: usize) -> f64 {
        self.elastic[self.voxel_idx(i, j, k)]
    }

    fn crystal_pix_per_um(&self) -> f64 {
        self.pix_per_um
    }

    fn coefcalc(&self) -> &dyn CoefCalc {
        &*self.coefcalc
    }

    fn coefcalc_mut(&mut self) -> &mut dyn CoefCalc {
        &mut *self.coefcalc
    }

    fn ddm(&self) -> &dyn DdmModel {
        &*self.ddm
    }

    fn exposure_summary(&self) -> &ExposureSummary {
        &self.exposure_summary
    }

    fn exposure_summary_mut(&mut self) -> &mut ExposureSummary {
        &mut self.exposure_summary
    }

    fn subprogram(&self) -> &str {
        &self.subprogram
    }

    fn photo_electron_escape(&self) -> bool {
        self.photo_electron_escape
    }

    fn fluorescent_escape(&self) -> bool {
        self.fluorescent_escape
    }

    fn calc_surrounding(&self) -> bool {
        self.calc_surrounding
    }

    fn expose(&mut self, beam: &mut dyn Beam, wedge: &Wedge) {
        match self.subprogram.as_str() {
            "RD3D" | "" => {
                let mut container: Box<dyn Container> =
                    std::mem::replace(&mut self.container, Box::new(ContainerTransparent));
                super::expose_rd3d(self, beam, wedge, &mut *container);
                self.container = container;
            }
            other => {
                unreachable!(
                    "unrecognised subprogram '{}' should have been caught at construction",
                    other
                );
            }
        }
    }
}

// ── Shared geometry helpers ────────────────────────────────────────────────────

/// Build expanded vertex array for fast triangle lookup.
fn build_expanded_vertices(vertices: &[[f64; 3]], indices: &[[usize; 3]]) -> Vec<[[f64; 3]; 3]> {
    indices
        .iter()
        .map(|tri| [vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]])
        .collect()
}

/// Calculate triangle normal (unit) and signed distance from origin: d = -(n · v0).
fn calc_triangle_normal_and_dist(vertices: &[[f64; 3]], tri: &[usize; 3]) -> ([f64; 3], f64) {
    let v0 = vertices[tri[0]];
    let v1 = vertices[tri[1]];
    let v2 = vertices[tri[2]];

    let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];

    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];

    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    let normal = if len > 0.0 {
        [n[0] / len, n[1] / len, n[2] / len]
    } else {
        [0.0, 0.0, 0.0]
    };

    let d = -(normal[0] * v0[0] + normal[1] * v0[1] + normal[2] * v0[2]);
    (normal, d)
}

/// 2D polygon inclusion test (ray-casting in X-Y plane).
/// Matches Java's `Vector.polygonInclusionTest`.
fn point_in_triangle(p: &[f64; 3], a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> bool {
    let verts = [*a, *b, *c];
    let mut inside = false;
    let n = verts.len();
    let mut j = n - 1;
    for i in 0..n {
        if ((verts[i][1] > p[1]) != (verts[j][1] > p[1]))
            && (p[0]
                < (verts[j][0] - verts[i][0]) * (p[1] - verts[i][1]) / (verts[j][1] - verts[i][1])
                    + verts[i][0])
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::Crystal;
    use crate::parser::config::{CrystalConfig, CrystalType};

    fn make_cylinder_config(diameter: f64, height: f64) -> CrystalConfig {
        CrystalConfig {
            crystal_type: Some(CrystalType::Cylinder),
            dim_x: Some(diameter),
            dim_y: Some(height),
            // Minimal unit cell for CoefCalcFromParams
            cell_a: Some(78.0),
            cell_b: Some(78.0),
            cell_c: Some(37.0),
            num_residues: Some(60),
            ..Default::default()
        }
    }

    fn make_spherical_new_config(diameter: f64) -> CrystalConfig {
        CrystalConfig {
            crystal_type: Some(CrystalType::SphericalNew),
            dim_x: Some(diameter),
            // Minimal unit cell for CoefCalcFromParams
            cell_a: Some(78.0),
            cell_b: Some(78.0),
            cell_c: Some(37.0),
            num_residues: Some(60),
            ..Default::default()
        }
    }

    #[test]
    fn test_cylinder_geometry_is_valid() {
        let diameter = 100.0f64;
        let height = 150.0f64;
        let config = make_cylinder_config(diameter, height);
        let crystal =
            crystal_cylinder_from_config(&config).expect("Cylinder crystal should construct");

        // Bounding box should be approximately [diameter, height, diameter]
        let [sx, sy, sz] = crystal.cryst_size_um;
        assert!(
            (sx - diameter).abs() < 2.0,
            "x size {:.1} ≈ diameter {:.1}",
            sx,
            diameter
        );
        assert!(
            (sy - height).abs() < 2.0,
            "y size {:.1} ≈ height {:.1}",
            sy,
            height
        );
        assert!(
            (sz - diameter).abs() < 2.0,
            "z size {:.1} ≈ diameter {:.1}",
            sz,
            diameter
        );

        // Center voxel should be inside
        let [nx, ny, nz] = crystal.cryst_size_voxels;
        let ci = nx / 2;
        let cj = ny / 2;
        let ck = nz / 2;
        assert!(
            crystal.is_crystal_at(ci, cj, ck),
            "Center voxel ({},{},{}) should be inside cylinder",
            ci,
            cj,
            ck
        );

        // Corner voxel (0,0,0) should be outside (corner of bounding box)
        assert!(
            !crystal.is_crystal_at(0, 0, 0),
            "Corner voxel (0,0,0) should be outside cylinder"
        );
    }

    // ── Helpers for polyhedron OBJ tests ─────────────────────────────────────────

    fn make_polyhedron_config(obj_path: &str) -> CrystalConfig {
        CrystalConfig {
            crystal_type: Some(CrystalType::Polyhedron),
            model_file: Some(obj_path.to_string()),
            pixels_per_micron: Some(0.5),
            angle_p: Some(0.0),
            angle_l: Some(0.0),
            dim_x: Some(60.0),
            dim_y: Some(20.0),
            dim_z: Some(40.0),
            // Minimal unit cell for CoefCalcFromParams
            cell_a: Some(78.0),
            cell_b: Some(78.0),
            cell_c: Some(37.0),
            num_residues: Some(60),
            ..Default::default()
        }
    }

    fn zero_wedge() -> crate::wedge::Wedge {
        crate::wedge::Wedge {
            start_ang: 0.0,
            end_ang: 0.0,
            ang_res: 2.0_f64.to_radians(),
            exposure_time: 100.0,
            start_x: 0.0,
            start_y: 0.0,
            start_z: 0.0,
            trans_x: 0.0,
            trans_y: 0.0,
            trans_z: 0.0,
            off_axis_um: 0.0,
            max_resolution: 2.0,
        }
    }

    /// Port of Java CrystalPolyhedronTests.testFindDepthSimple:
    /// Verify get_cryst_coord and find_depth for a cuboid OBJ at angles 0° and 90°.
    #[test]
    fn test_polyhedron_cryst_coord_and_depth_simple() {
        let manifest = env!("CARGO_MANIFEST_DIR");
        let obj_path = format!(
            "{}/../../tests/fixtures/CrystalPolyhedron-cuboid-30-20-10.obj",
            manifest
        );
        let config = make_polyhedron_config(&obj_path);

        let mut crystal =
            CrystalPolyhedron::from_config(&config).expect("Cuboid OBJ polyhedron should build");

        let wedge = zero_wedge();

        // Crystal dims from OBJ (≈60×20×40 µm); resolution = pix_per_um = 0.5
        let xdim = 60.0f64;
        let ydim = 20.0f64;
        let zdim = 40.0f64;
        let resolution = 0.5f64; // pixels per µm

        // Inner voxel range: 1..(dim*resolution - 1), matching Java loop
        let x_max = (xdim * resolution) as usize - 1; // 29
        let y_max = (ydim * resolution) as usize - 1; // 9
        let z_max = (zdim * resolution) as usize - 1; // 19

        // --- Test get_cryst_coord ---
        for x in 1..x_max {
            for y in 1..y_max {
                for z in 1..z_max {
                    let coord = crystal.get_cryst_coord(x, y, z);
                    let ex = -(xdim / 2.0) + x as f64 / resolution;
                    let ey = -(ydim / 2.0) + y as f64 / resolution;
                    let ez = -(zdim / 2.0) + z as f64 / resolution;
                    assert!(
                        (coord[0] - ex).abs() < 0.1,
                        "cryst_coord[0] for ({},{},{}): got {:.4}, expected {:.4}",
                        x,
                        y,
                        z,
                        coord[0],
                        ex
                    );
                    assert!(
                        (coord[1] - ey).abs() < 0.1,
                        "cryst_coord[1] for ({},{},{}): got {:.4}, expected {:.4}",
                        x,
                        y,
                        z,
                        coord[1],
                        ey
                    );
                    assert!(
                        (coord[2] - ez).abs() < 0.1,
                        "cryst_coord[2] for ({},{},{}): got {:.4}, expected {:.4}",
                        x,
                        y,
                        z,
                        coord[2],
                        ez
                    );
                }
            }
        }

        // --- Test find_depth at angle = 0° ---
        crystal.setup_depth_finding(0.0, &wedge);
        for x in 1..x_max {
            for y in 1..y_max {
                for z in 1..z_max {
                    let coord = crystal.get_cryst_coord(x, y, z);
                    let depth = crystal.find_depth(&coord, 0.0, &wedge);
                    let true_depth = z as f64 / resolution;
                    assert!(
                        (depth - true_depth).abs() < 2.0,
                        "find_depth angle=0 for ({},{},{}): got {:.4}, expected {:.4}",
                        x,
                        y,
                        z,
                        depth,
                        true_depth
                    );
                }
            }
        }

        // --- Test find_depth at angle = 90° ---
        // Java swaps crystCoords[0] and crystCoords[2] before calling findDepth(90°).
        // This mirrors what exposeAngle does when rotating the crystal in the beam.
        crystal.setup_depth_finding(std::f64::consts::FRAC_PI_2, &wedge);
        for x in 1..x_max {
            for y in 1..y_max {
                for z in 1..z_max {
                    let orig = crystal.get_cryst_coord(x, y, z);
                    // Swap x and z as Java does
                    let rotated_coord = [orig[2], orig[1], orig[0]];
                    let depth = crystal.find_depth(&rotated_coord, 0.0, &wedge);
                    let true_depth = x as f64 / resolution;
                    assert!(
                        (depth - true_depth).abs() < 2.0,
                        "find_depth angle=90 for ({},{},{}): got {:.4}, expected {:.4}",
                        x,
                        y,
                        z,
                        depth,
                        true_depth
                    );
                }
            }
        }
    }

    /// Port of Java CrystalPolyhedronTests.testFindDepthConcave:
    /// Verify find_depth for a horseshoe-shaped concave crystal — thick region
    /// should return ≈60 µm and gap region ≈40 µm (skipping the hollow interior).
    #[test]
    fn test_polyhedron_find_depth_concave() {
        let manifest = env!("CARGO_MANIFEST_DIR");
        let obj_path = format!(
            "{}/../../tests/fixtures/CrystalPolyhedron-concave_cuboid-30-20-10.obj",
            manifest
        );

        let config = CrystalConfig {
            crystal_type: Some(CrystalType::Polyhedron),
            model_file: Some(obj_path),
            angle_p: Some(0.0),
            angle_l: Some(0.0),
            dim_x: Some(27.5),
            dim_y: Some(41.44),
            dim_z: Some(62.67),
            cell_a: Some(78.0),
            cell_b: Some(78.0),
            cell_c: Some(37.0),
            num_residues: Some(60),
            ..Default::default()
        };

        let mut crystal =
            CrystalPolyhedron::from_config(&config).expect("Concave OBJ polyhedron should build");

        let wedge = zero_wedge();
        crystal.setup_depth_finding(0.0, &wedge);

        // Thick region: beam traverses a solid arm (~60 µm)
        let thick_coord = [3.5f64, -8.65, 29.9];
        let thick_depth = crystal.find_depth(&thick_coord, 0.0, &wedge);
        assert!(
            (thick_depth - 60.0).abs() < 1.0,
            "Concave crystal thick depth: got {:.4}, expected 60.0 ± 1.0",
            thick_depth
        );

        // Gap region: beam skips the hollow interior, total crystal ~40 µm
        let thin_coord = [5.0f64, 6.0, 30.0];
        let thin_depth = crystal.find_depth(&thin_coord, 0.0, &wedge);
        assert!(
            (thin_depth - 40.0).abs() < 1.0,
            "Concave crystal thin depth: got {:.4}, expected 40.0 ± 1.0",
            thin_depth
        );
    }

    #[test]
    fn test_spherical_new_occupancy() {
        let diameter = 100.0f64;
        let config = make_spherical_new_config(diameter);
        let crystal =
            crystal_spherical_new_from_config(&config).expect("SphericalNew should construct");

        // Center voxel should be inside
        let [nx, ny, nz] = crystal.cryst_size_voxels;
        let ci = nx / 2;
        let cj = ny / 2;
        let ck = nz / 2;
        assert!(
            crystal.is_crystal_at(ci, cj, ck),
            "Center voxel should be inside SphericalNew"
        );

        // Extreme corner should be outside
        assert!(
            !crystal.is_crystal_at(0, 0, 0),
            "Corner voxel (0,0,0) should be outside SphericalNew"
        );
    }

    const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures");

    fn make_polyhedron_from_obj(
        obj_file: &str,
        dim_x: f64,
        dim_y: f64,
        dim_z: f64,
        resolution: f64,
    ) -> CrystalPolyhedron {
        let config = CrystalConfig {
            crystal_type: Some(CrystalType::Polyhedron),
            coefcalc: Some(crate::parser::config::CoefCalcType::Average),
            dim_x: Some(dim_x),
            dim_y: Some(dim_y),
            dim_z: Some(dim_z),
            pixels_per_micron: Some(resolution),
            angle_p: Some(0.0),
            angle_l: Some(0.0),
            model_file: Some(format!("{}/{}", FIXTURES_DIR, obj_file)),
            cell_a: Some(100.0),
            cell_b: Some(100.0),
            cell_c: Some(100.0),
            ..Default::default()
        };
        CrystalPolyhedron::from_config(&config).unwrap()
    }

    #[test]
    fn polyhedron_depth_matches_cuboid_for_box_obj() {
        use crate::crystal::cuboid::CrystalCuboid;

        let xdim = 60.0;
        let ydim = 20.0;
        let zdim = 40.0;
        let resolution = 0.5;

        // OBJ file has vertices at ±30,±10,±20 — already 60×20×40 µm.
        // dim_x is not used for OBJ scaling (matches Java behavior).
        let mut poly = make_polyhedron_from_obj(
            "CrystalPolyhedron-cuboid-30-20-10.obj",
            xdim,
            ydim,
            zdim,
            resolution,
        );

        let cub_config = CrystalConfig {
            crystal_type: Some(CrystalType::Cuboid),
            coefcalc: Some(crate::parser::config::CoefCalcType::Average),
            dim_x: Some(xdim),
            dim_y: Some(ydim),
            dim_z: Some(zdim),
            pixels_per_micron: Some(resolution),
            angle_p: Some(0.0),
            angle_l: Some(0.0),
            cell_a: Some(100.0),
            cell_b: Some(100.0),
            cell_c: Some(100.0),
            ..Default::default()
        };
        let mut cub = CrystalCuboid::from_config(&cub_config).unwrap();

        let w = crate::wedge::Wedge::from_config(&crate::parser::config::WedgeConfig {
            start_ang: 0.0,
            end_ang: 0.0,
            exposure_time: Some(100.0),
            angular_resolution: None,
            start_offset_x: None,
            start_offset_y: None,
            start_offset_z: None,
            translate_x: None,
            translate_y: None,
            translate_z: None,
            rot_ax_beam_offset: None,
            max_resolution: Some(2.0),
        });

        let [nx, ny, nz] = poly.cryst_size_voxels();

        // Check interior voxels (skip edges due to rounding)
        poly.setup_depth_finding(0.0, &w);
        cub.setup_depth_finding(0.0, &w);

        for x in 1..(nx - 1) {
            for y in 1..(ny - 1) {
                for z in 1..(nz - 1) {
                    let pc = poly.get_cryst_coord(x, y, z);
                    let cc = cub.get_cryst_coord(x, y, z);

                    let depth_poly = poly.find_depth(&pc, 0.0, &w);
                    let depth_cub = cub.find_depth(&cc, 0.0, &w);

                    let expected = z as f64 / resolution;
                    assert!(
                        (depth_poly - expected).abs() < 2.0,
                        "Polyhedron depth at z={}: got {}, expected {}",
                        z,
                        depth_poly,
                        expected
                    );
                    assert!(
                        (depth_cub - expected).abs() < 2.0,
                        "Cuboid depth at z={}: got {}, expected {}",
                        z,
                        depth_cub,
                        expected
                    );
                }
            }
        }
    }

    #[test]
    fn polyhedron_concave_depth() {
        let resolution = 0.5;
        let mut poly = make_polyhedron_from_obj(
            "CrystalPolyhedron-concave_cuboid-30-20-10.obj",
            1.0, // scale = 1.0 (OBJ coordinates used directly)
            1.0,
            1.0,
            resolution,
        );

        let w = crate::wedge::Wedge::from_config(&crate::parser::config::WedgeConfig {
            start_ang: 0.0,
            end_ang: 0.0,
            exposure_time: Some(100.0),
            angular_resolution: None,
            start_offset_x: None,
            start_offset_y: None,
            start_offset_z: None,
            translate_x: None,
            translate_y: None,
            translate_z: None,
            rot_ax_beam_offset: None,
            max_resolution: Some(2.0),
        });

        poly.setup_depth_finding(0.0, &w);

        // Thick part of crystal
        let thick_coord = [3.5, -8.65, 29.9];
        let thick_depth = poly.find_depth(&thick_coord, 0.0, &w);
        assert!(
            (thick_depth - 60.0).abs() < 2.0,
            "Thick part depth: got {}, expected ~60",
            thick_depth
        );

        // Thin part (through horseshoe gap)
        let thin_coord = [5.0, 6.0, 30.0];
        let thin_depth = poly.find_depth(&thin_coord, 0.0, &w);
        assert!(
            (thin_depth - 40.0).abs() < 2.0,
            "Thin part depth: got {}, expected ~40",
            thin_depth
        );
    }
}
