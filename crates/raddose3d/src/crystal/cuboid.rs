use crate::beam::Beam;
use crate::coefcalc::{self, CoefCalc};
use crate::container::{self, Container, ContainerTransparent};
use crate::ddm::{self, DdmModel};
use crate::output::ExposureSummary;
use crate::parser::config::CrystalConfig;
use crate::simulation::{MicroEdSimulation, MonteCarloSimulation, XfelSimulation};
use crate::wedge::Wedge;

/// Cuboid crystal with rectangular voxel grid.
///
/// In the Java code, CrystalCuboid extends CrystalPolyhedron which uses
/// ray-tracing for occupancy. For a cuboid, occupancy is trivially determined
/// by bounding box, so we implement it directly with a simpler approach.
#[derive(Debug)]
#[allow(dead_code)] // phase-6 fields (angle_p, angle_l, photo_electron_escape, fluorescent_escape, first_wedge)
pub struct CrystalCuboid {
    /// Crystal dimensions in µm [x, y, z].
    cryst_size_um: [f64; 3],
    /// Crystal dimensions in voxels [nx, ny, nz].
    cryst_size_voxels: [usize; 3],
    /// Resolution in voxels per µm.
    pix_per_um: f64,
    /// Initial rotation angles (radians).
    angle_p: f64,
    angle_l: f64,
    /// 3D dose array (MGy).
    dose: Vec<f64>,
    /// 3D fluence array.
    fluence: Vec<f64>,
    /// 3D elastic scattering array.
    elastic: Vec<f64>,
    /// Voxel coordinates: stored as [x, y, z] for each voxel.
    cryst_coord: Vec<[f64; 3]>,
    /// Coefficient calculator.
    coefcalc: Box<dyn CoefCalc>,
    /// Dose decay model.
    ddm: Box<dyn DdmModel>,
    /// Container.
    container: Box<dyn Container>,
    /// Exposure summary.
    exposure_summary: ExposureSummary,
    /// Subprogram mode.
    subprogram: String,
    /// Whether PE escape is enabled.
    photo_electron_escape: bool,
    /// Whether fluorescent escape is enabled.
    fluorescent_escape: bool,
    /// Is this the first wedge?
    first_wedge: bool,
    /// Vertices of the cuboid (8 vertices × 3 coords).
    vertices: [[f64; 3]; 8],
    /// Rotated vertices for current angle.
    rotated_vertices: [[f64; 3]; 8],
    /// Triangle indices (12 triangles × 3 vertex indices).
    indices: [[usize; 3]; 12],
    /// Triangle normals (12 triangles).
    normals: [[f64; 3]; 12],
    /// Rotated normals.
    rotated_normals: [[f64; 3]; 12],
    /// Origin distances for each triangle plane.
    origin_distances: [f64; 12],
    /// Rotated origin distances.
    rotated_origin_distances: [f64; 12],
    // ── MC/XFEL/MicroED simulation config ──────────────────────────────────
    /// Number of MC/XFEL runs (default 1).
    runs: usize,
    /// Vertical goniometer axis (true = 90°).
    vertical_goniometer: bool,
    /// Vertical polarisation (true = 90°).
    vertical_polarisation: bool,
    /// Surrounding cryoprotectant thickness [x, y, z] in µm.
    surrounding_thickness: [f64; 3],
    /// Crystal type string for MicroED ("CUBOID", etc.).
    crystal_type: String,
}

impl CrystalCuboid {
    /// Default resolution in voxels/µm.
    const DEFAULT_RESOLUTION: f64 = 0.5;
    /// Max voxels before reducing resolution.
    #[allow(dead_code)]
    const MAX_VOXELS: usize = 1_000_000;

    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let dim_x = config.dim_x.ok_or("Crystal requires DimX")?;
        let dim_y = config.dim_y.unwrap_or(dim_x);
        let dim_z = config.dim_z.unwrap_or(dim_y);

        // Resolution
        let pix_per_um = config.pixels_per_micron.unwrap_or_else(|| {
            let default_res = 10.0 / dim_x;
            default_res.min(Self::DEFAULT_RESOLUTION)
        });

        // Crystal size in voxels
        let nx = (dim_x * pix_per_um).round() as usize + 1;
        let ny = (dim_y * pix_per_um).round() as usize + 1;
        let nz = (dim_z * pix_per_um).round() as usize + 1;

        let _total_voxels = nx * ny * nz;

        // Angles
        let angle_p = config.angle_p.unwrap_or(0.0).to_radians();
        let angle_l = config.angle_l.unwrap_or(0.0).to_radians();

        // Build vertices for cuboid centered at origin
        let hx = dim_x / 2.0;
        let hy = dim_y / 2.0;
        let hz = dim_z / 2.0;

        let mut vertices = [
            [-hx, -hy, hz],
            [-hx, -hy, -hz],
            [-hx, hy, -hz],
            [-hx, hy, hz],
            [hx, -hy, hz],
            [hx, -hy, -hz],
            [hx, hy, -hz],
            [hx, hy, hz],
        ];

        // Apply initial rotations (P about Z, L about X)
        if angle_p != 0.0 || angle_l != 0.0 {
            let cos_p = angle_p.cos();
            let sin_p = angle_p.sin();
            let cos_l = angle_l.cos();
            let sin_l = angle_l.sin();
            for v in &mut vertices {
                // Rotate about Z
                let x = v[0] * cos_p - v[1] * sin_p;
                let y = v[0] * sin_p + v[1] * cos_p;
                v[0] = x;
                v[1] = y;
                // Rotate about X
                let y2 = v[1] * cos_l - v[2] * sin_l;
                let z2 = v[1] * sin_l + v[2] * cos_l;
                v[1] = y2;
                v[2] = z2;
            }
        }

        // Triangle indices (1-based in Java, 0-based here)
        let indices: [[usize; 3]; 12] = [
            [0, 2, 1],
            [3, 2, 0],
            [2, 5, 1],
            [6, 5, 2],
            [1, 4, 0],
            [1, 5, 4],
            [3, 7, 2],
            [7, 6, 2],
            [3, 0, 7],
            [0, 4, 7],
            [7, 4, 6],
            [6, 4, 5],
        ];

        // Calculate normals and origin distances
        let mut normals = [[0.0; 3]; 12];
        let mut origin_distances = [0.0; 12];
        for (t, tri) in indices.iter().enumerate() {
            let (n, d) = calc_triangle_normal_and_dist(&vertices, tri);
            normals[t] = n;
            origin_distances[t] = d;
        }

        // Initialize voxel coordinates
        // Bounding box from rotated vertices
        let mut min_coords = [f64::INFINITY; 3];
        let mut max_coords = [f64::NEG_INFINITY; 3];
        for v in &vertices {
            for dim in 0..3 {
                min_coords[dim] = min_coords[dim].min(v[dim]);
                max_coords[dim] = max_coords[dim].max(v[dim]);
            }
        }

        let actual_dim_x = max_coords[0] - min_coords[0];
        let actual_dim_y = max_coords[1] - min_coords[1];
        let actual_dim_z = max_coords[2] - min_coords[2];

        let nx = (actual_dim_x * pix_per_um).round() as usize + 1;
        let ny = (actual_dim_y * pix_per_um).round() as usize + 1;
        let nz = (actual_dim_z * pix_per_um).round() as usize + 1;
        let total_voxels = nx * ny * nz;

        let x_shift = actual_dim_x / 2.0;
        let y_shift = actual_dim_y / 2.0;
        let z_shift = actual_dim_z / 2.0;

        let mut cryst_coord = vec![[0.0; 3]; total_voxels];
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = i * ny * nz + j * nz + k;
                    cryst_coord[idx] = [
                        -x_shift + i as f64 / pix_per_um,
                        -y_shift + j as f64 / pix_per_um,
                        -z_shift + k as f64 / pix_per_um,
                    ];
                }
            }
        }

        // Create CoefCalc
        let coefcalc = coefcalc::create_coefcalc(config)?;

        // Create container
        let container = container::create_container(config);

        // Create DDM
        let ddm = ddm::create_ddm(
            config.ddm,
            config.gamma_param,
            config.b0_param,
            config.beta_param,
        );

        // Subprogram
        let subprogram = config.program.as_deref().unwrap_or("RD3D").to_uppercase();

        let photo_electron_escape = config
            .calculate_pe_escape
            .as_deref()
            .is_some_and(|s| s.eq_ignore_ascii_case("true"));
        let fluorescent_escape = config
            .calculate_fl_escape
            .as_deref()
            .is_some_and(|s| s.eq_ignore_ascii_case("true"));

        // MC/XFEL/MicroED config
        let runs = config.runs.unwrap_or(1).max(1) as usize;
        let vertical_goniometer = config
            .goniometer_axis
            .map(|v| (v - 90.0).abs() < 1e-6)
            .unwrap_or(false);
        let vertical_polarisation = config
            .polarisation_direction
            .map(|v| (v - 90.0).abs() < 1e-6)
            .unwrap_or(false);
        let surrounding_thickness = config.surrounding_thickness.unwrap_or([0.0; 3]);

        Ok(CrystalCuboid {
            cryst_size_um: [actual_dim_x, actual_dim_y, actual_dim_z],
            cryst_size_voxels: [nx, ny, nz],
            pix_per_um,
            angle_p,
            angle_l,
            dose: vec![0.0; total_voxels],
            fluence: vec![0.0; total_voxels],
            elastic: vec![0.0; total_voxels],
            cryst_coord,
            coefcalc,
            ddm,
            container,
            exposure_summary: ExposureSummary::new(),
            subprogram,
            photo_electron_escape,
            fluorescent_escape,
            first_wedge: true,
            vertices,
            rotated_vertices: vertices,
            indices,
            normals,
            rotated_normals: normals,
            origin_distances,
            rotated_origin_distances: origin_distances,
            runs,
            vertical_goniometer,
            vertical_polarisation,
            surrounding_thickness,
            crystal_type: "CUBOID".to_string(),
        })
    }

    /// Linear index into flat voxel arrays.
    #[inline]
    fn voxel_idx(&self, i: usize, j: usize, k: usize) -> usize {
        let [_, ny, nz] = self.cryst_size_voxels;
        i * ny * nz + j * nz + k
    }

    /// Check if voxel is inside the cuboid using ray-casting.
    fn is_inside_polyhedron(&self, coord: &[f64; 3]) -> bool {
        // Ray direction along Z axis
        let dir = [0.0, 0.0, 1.0];
        let mut crossings = 0;

        for (t, tri) in self.indices.iter().enumerate() {
            let n = &self.normals[t];
            let d = self.origin_distances[t];

            // Ray-plane intersection
            let denom = n[0] * dir[0] + n[1] * dir[1] + n[2] * dir[2];
            if denom.abs() < 1e-12 {
                continue;
            }

            let t_val = (n[0] * coord[0] + n[1] * coord[1] + n[2] * coord[2] + d) / denom;
            if t_val <= 0.0 || t_val.is_nan() || t_val.is_infinite() {
                continue;
            }

            // Intersection point
            let p = [
                coord[0] + t_val * dir[0],
                coord[1] + t_val * dir[1],
                coord[2] + t_val * dir[2],
            ];

            // Check if point is inside triangle
            if point_in_triangle(
                &p,
                &self.vertices[tri[0]],
                &self.vertices[tri[1]],
                &self.vertices[tri[2]],
            ) {
                crossings += 1;
            }
        }

        // Odd number of crossings = inside
        crossings % 2 == 1
    }
}

impl super::Crystal for CrystalCuboid {
    fn setup_depth_finding(&mut self, ang_rad: f64, wedge: &Wedge) {
        let cos_a = ang_rad.cos();
        let sin_a = ang_rad.sin();

        for (i, v) in self.vertices.iter().enumerate() {
            let y = v[1] + wedge.start_y + wedge.trans_y * (ang_rad - wedge.start_ang);
            let x = v[0] + wedge.start_x + wedge.trans_x * (ang_rad - wedge.start_ang);
            let z = v[2] + wedge.start_z + wedge.trans_z * (ang_rad - wedge.start_ang);

            // Rotate in X-Z plane
            self.rotated_vertices[i] = [x * cos_a + z * sin_a, y, -x * sin_a + z * cos_a];
        }

        // Recalculate normals for rotated vertices
        for (t, tri) in self.indices.iter().enumerate() {
            let (n, d) = calc_triangle_normal_and_dist(&self.rotated_vertices, tri);
            self.rotated_normals[t] = n;
            self.rotated_origin_distances[t] = d;
        }
    }

    fn find_depth(&self, vox_coord: &[f64; 3], _delta_phi: f64, _wedge: &Wedge) -> f64 {
        // Cast ray in Z direction, find intersections with polyhedron surface
        let dir = [0.0, 0.0, 1.0];
        let mut distances = Vec::new();

        for (t, tri) in self.indices.iter().enumerate() {
            let n = &self.rotated_normals[t];
            let d = self.rotated_origin_distances[t];

            let denom = n[0] * dir[0] + n[1] * dir[1] + n[2] * dir[2];
            if denom.abs() < 1e-12 {
                continue;
            }

            let t_val =
                (n[0] * vox_coord[0] + n[1] * vox_coord[1] + n[2] * vox_coord[2] + d) / denom;
            if t_val <= 0.0 || t_val.is_nan() || t_val.is_infinite() {
                continue;
            }

            let p = [
                vox_coord[0] + t_val * dir[0],
                vox_coord[1] + t_val * dir[1],
                vox_coord[2] + t_val * dir[2],
            ];

            if point_in_triangle(
                &p,
                &self.rotated_vertices[tri[0]],
                &self.rotated_vertices[tri[1]],
                &self.rotated_vertices[tri[2]],
            ) {
                distances.push(t_val);
            }
        }

        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // Remove duplicates
        distances.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

        // Even number of intersections (or none) = outside crystal, return 0
        if distances.is_empty() || distances.len() % 2 == 0 {
            return 0.0;
        }

        // Depth = distance to first surface + gaps for re-entry
        // This matches Java: depth = d[0] + (d[2]-d[1]) + (d[4]-d[3]) + ...
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
        format!(
            "Cuboid crystal of size [{:.0}, {:.0}, {:.0}] um [x, y, z] at a resolution of {:.2} microns per voxel edge.",
            self.cryst_size_um[0],
            self.cryst_size_um[1],
            self.cryst_size_um[2],
            1.0 / self.pix_per_um
        )
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

    fn expose(&mut self, beam: &mut dyn Beam, wedge: &Wedge) {
        eprintln!("{}", self.crystal_info());

        match self.subprogram.as_str() {
            "RD3D" | "" => {
                // We need to take container out temporarily to avoid borrow issues
                let mut container: Box<dyn Container> =
                    std::mem::replace(&mut self.container, Box::new(ContainerTransparent));
                super::expose_rd3d(self, beam, wedge, &mut *container);
                self.container = container;
            }
            "MONTECARLO" | "GOS" => {
                let gos = self.subprogram == "GOS";
                let vertices: Vec<[f64; 3]> = self.vertices.to_vec();
                let indices: Vec<[usize; 3]> = self.indices.to_vec();
                let n = self.cryst_size_voxels[0]
                    * self.cryst_size_voxels[1]
                    * self.cryst_size_voxels[2];
                let cryst_occ = vec![true; n];
                for run_idx in 0..self.runs {
                    let mut mc = MonteCarloSimulation::new(
                        vertices.clone(),
                        indices.clone(),
                        self.cryst_coord.clone(),
                        self.pix_per_um,
                        self.cryst_size_voxels,
                        cryst_occ.clone(),
                        run_idx + 1,
                        self.vertical_goniometer,
                        false, // xfel = false for MC
                        gos,
                        self.surrounding_thickness,
                        self.vertical_polarisation,
                    );
                    mc.populate_auger_linewidths();
                    mc.populate_fluorescence_linewidths();
                    mc.populate_angular_emission_probs();
                    mc.calculate_xfel(beam, wedge, &mut *self.coefcalc);
                }
            }
            "XFEL" => {
                let vertices: Vec<[f64; 3]> = self.vertices.to_vec();
                let indices: Vec<[usize; 3]> = self.indices.to_vec();
                let n = self.cryst_size_voxels[0]
                    * self.cryst_size_voxels[1]
                    * self.cryst_size_voxels[2];
                let cryst_occ = vec![true; n];
                for run_idx in 0..self.runs {
                    let mut xfel = XfelSimulation::new(
                        vertices.clone(),
                        indices.clone(),
                        self.cryst_coord.clone(),
                        self.pix_per_um,
                        self.cryst_size_voxels,
                        cryst_occ.clone(),
                        run_idx + 1,
                        self.vertical_goniometer,
                        true, // xfel = true
                        true, // gos = true
                        wedge.total_sec(),
                        self.vertical_polarisation,
                    );
                    xfel.populate_auger_linewidths();
                    xfel.populate_fluorescence_linewidths();
                    xfel.calculate_xfel(beam, wedge, &mut *self.coefcalc);
                }
            }
            "EMSP" | "EMED" => {
                let vertices: Vec<[f64; 3]> = self.vertices.to_vec();
                let indices: Vec<[usize; 3]> = self.indices.to_vec();
                let n = self.cryst_size_voxels[0]
                    * self.cryst_size_voxels[1]
                    * self.cryst_size_voxels[2];
                // MicroED uses Vec<(bool, bool)> for occupancy
                let cryst_occ = vec![(true, true); n];
                let mut micro_ed = MicroEdSimulation::new(
                    vertices,
                    indices,
                    self.cryst_coord.clone(),
                    self.pix_per_um,
                    self.cryst_size_voxels,
                    cryst_occ,
                    self.crystal_type.clone(),
                );
                micro_ed.calculate_em(beam, wedge, &mut *self.coefcalc);
            }
            other => {
                eprintln!("Subprogram '{}' not yet implemented, using RD3D", other);
                let mut container: Box<dyn Container> =
                    std::mem::replace(&mut self.container, Box::new(ContainerTransparent));
                super::expose_rd3d(self, beam, wedge, &mut *container);
                self.container = container;
            }
        }
    }
}

/// Calculate triangle normal and signed distance from origin.
fn calc_triangle_normal_and_dist(vertices: &[[f64; 3]], tri: &[usize; 3]) -> ([f64; 3], f64) {
    let v0 = vertices[tri[0]];
    let v1 = vertices[tri[1]];
    let v2 = vertices[tri[2]];

    // Edge vectors
    let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];

    // Cross product
    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];

    // Normalize
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    let normal = if len > 0.0 {
        [n[0] / len, n[1] / len, n[2] / len]
    } else {
        [0.0, 0.0, 0.0]
    };

    // Distance from origin: d = -(normal · v0)  (matches Java convention)
    let d = -(normal[0] * v0[0] + normal[1] * v0[1] + normal[2] * v0[2]);

    (normal, d)
}

/// Check if a point lies inside a triangle (3D, assumes point is on the triangle plane).
/// 2D polygon inclusion test (ray-casting in X-Y plane).
/// Matches Java's `Vector.polygonInclusionTest`.
fn point_in_triangle(p: &[f64; 3], a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> bool {
    let verts: [[f64; 3]; 3] = [*a, *b, *c];
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

// Overriding the is_crystal_at for cuboid: we can use the indexed normals approach
// from CrystalPolyhedron, but for a simple cuboid we use the same ray-casting approach
// as the Java code.

// For the 8-vertex cuboid, note that occupancy checking using the polyhedron approach
// with 12 triangles is equivalent to checking if a point is within the rotated bounding box.
