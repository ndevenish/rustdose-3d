// Many variables in this module are intermediate simulation-state trackers
// that will be connected to outputs once the full MC pipeline is wired up.
#![allow(dead_code, unused_variables, unused_assignments)]
/// Monte Carlo simulation engine for X-ray damage calculation.
///
/// Ports Java MC.java (~4808 lines).
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;

use crate::beam::Beam;
use crate::coefcalc::CoefCalc;
use crate::element::database::Element;
use crate::wedge::Wedge;

/// Speed of light (m/s).
pub const C: f64 = 299_792_458.0;
/// Electron rest mass (kg).
pub const M_E: f64 = 9.109_383_56e-31;
/// Planck constant (J·s).
pub const H: f64 = 6.626_070_04e-34;
/// keV to Joules conversion.
pub const KEV_TO_JOULES: f64 = 1.602_176_565e-16;

/// Atomic numbers supporting Auger cascade.
const AUGER_ELEMENTS: [u32; 19] = [
    6, 7, 8, 11, 12, 14, 15, 16, 17, 19, 20, 25, 26, 27, 28, 29, 30, 33, 34,
];

/// Return number of Auger shells for a given element (by index into AUGER_ELEMENTS).
fn get_num_auger_shells(idx: usize) -> u32 {
    match AUGER_ELEMENTS[idx] {
        6..=8 => 0,
        11 | 12 => 1,
        14 | 15 | 16 | 17 | 19 | 20 => 2,
        25..=30 => 2,
        33 | 34 => 3,
        _ => 0,
    }
}

/// Auger/FL transition data for one element+shell.
#[derive(Default, Clone)]
pub struct TransitionData {
    pub linewidths: Vec<f64>,
    pub probs: Vec<f64>,
    pub energies: Vec<f64>,
    pub cumulative_probs: Vec<f64>,
    pub exit_index: Vec<f64>,
    pub drop_index: Vec<f64>,
}

/// Per-element transition map: shell → TransitionData.
type ShellTransitionMap = HashMap<u32, TransitionData>;
/// Full transition map: Z → ShellTransitionMap.
type TransitionMap = HashMap<u32, ShellTransitionMap>;

/// Monte Carlo simulation state.
pub struct MonteCarloSimulation {
    // Crystal geometry
    pub vertices: Vec<[f64; 3]>,
    pub rotated_vertices: Vec<[f64; 3]>,
    pub expanded_rotated_vertices: Vec<[[f64; 3]; 3]>,
    pub indices: Vec<[usize; 3]>,
    /// Flat voxel coordinates [nx*ny*nz][3] in µm.
    pub cryst_coord: Vec<[f64; 3]>,
    pub cryst_pix_per_um: f64,
    pub cryst_size_voxels: [usize; 3],
    /// Flat occupancy array.
    pub cryst_occ: Vec<bool>,
    pub run_number: usize,

    // Cached normals for ray-casting
    normals: Option<Vec<[f64; 3]>>,
    origin_distances: Option<Vec<f64>>,
    rotated_normals: Option<Vec<[f64; 3]>>,
    rotated_origin_distances: Option<Vec<f64>>,

    // Surrounding (cryoprotectant) volume
    pub vertices_surrounding: Vec<[f64; 3]>,
    pub indices_surrounding: Vec<[usize; 3]>,
    normals_surrounding: Option<Vec<[f64; 3]>>,
    origin_distances_surrounding: Option<Vec<f64>>,
    rotated_normals_surrounding: Option<Vec<[f64; 3]>>,
    rotated_origin_distances_surrounding: Option<Vec<f64>>,

    // Crystal dimensions in nm
    pub x_dimension: f64,
    pub y_dimension: f64,
    pub z_dimension: f64,

    // Bounding box
    max_dims: [f64; 3],
    min_dims: [f64; 3],

    // Aggregate dose metrics (keV, converted to MGy in process_dose)
    pub dose: f64,
    pub photon_dose: f64,
    pub electron_dose: f64,
    pub gos_electron_dose: f64,
    pub photon_dose_v_resolved: f64,
    pub gos_electron_dose_v_resolved: f64,
    pub electron_dose_surrounding: f64,
    pub raddose_style_dose: f64,
    pub raddose_style_dose_compton: f64,
    pub escaped_energy: f64,

    // Angular emission distribution
    angular_emission_probs: Vec<f64>,
    number_angular_emission_bins: usize,

    // Elastic angle lookup tables: [atomic_number] → BTreeMap<energy_bits, angle_probs>
    low_energy_angles: Vec<Option<BTreeMap<u64, Vec<f64>>>>,
    high_energy_angles: Vec<Option<BTreeMap<u64, Vec<f64>>>>,

    // Auger transition data
    pub auger_data: TransitionMap,
    pub tot_k_auger_prob: HashMap<u32, f64>,
    // Fluorescence transition data
    pub fl_data: TransitionMap,
    pub fl_cumulative: HashMap<u32, HashMap<u32, Vec<f64>>>,

    // Voxel arrays (flat, sized cryst_size_voxels[0]*[1]*[2])
    pub voxel_energy_v_resolved: Vec<f64>,
    pub voxel_ionisations_v_resolved: Vec<f64>,
    pub voxel_elastic: Vec<f64>,
    pub dose_simple: Vec<f64>,
    pub last_time_vox: Vec<f64>,

    // Ionisation tracking
    pub total_ionisation_events: u64,
    pub total_ionisation_events_v_resolved: u64,
    pub low_energy_ionisations: u64,
    ionisations_old: f64,
    pub atomic_ionisations: HashMap<String, u64>,
    pub atomic_ionisations_exposed: HashMap<String, u64>,
    pub atomic_ionisations_per_atom: HashMap<String, f64>,
    pub atomic_ionisations_per_atom_exposed: HashMap<String, f64>,

    // Energy loss straggling
    energy_per_inel: BTreeMap<u64, f64>,

    // Averages
    avg_w: f64,
    avg_w_num: u32,
    avg_uk: f64,
    avg_uk_num: u32,

    // Configuration
    pub vertical_goni: bool,
    pub vertical_pol: bool,
    pub simple_mc: bool,
    pub do_xfel: bool,
    tot_elastic: f64,

    // Simulation parameters
    pub num_photons: u64,
    pub pulse_length: f64,     // fs
    pub pulse_bin_length: f64, // fs
    pub pulse_energy: f64,     // J
    pub mean_energy: f64,      // keV
    pub energy_fwhm: Option<f64>,
    pub photon_energy_array: Vec<f64>,
    pub surrounding_thickness: [f64; 3], // nm

    pub last_time: f64,
}

// ── Helper: flat 3-D index ─────────────────────────────────────────────────────
fn vox3(size: &[usize; 3], i: usize, j: usize, k: usize) -> usize {
    i * size[1] * size[2] + j * size[2] + k
}

// ── f64 BTreeMap key helpers ───────────────────────────────────────────────────
fn f64_key(x: f64) -> u64 {
    x.to_bits()
}
fn btree_floor(map: &BTreeMap<u64, Vec<f64>>, key: f64) -> Option<(&u64, &Vec<f64>)> {
    map.range(..=f64_key(key)).next_back()
}
fn btree_ceil(map: &BTreeMap<u64, Vec<f64>>, key: f64) -> Option<(&u64, &Vec<f64>)> {
    map.range(f64_key(key)..).next()
}

impl MonteCarloSimulation {
    /// Construct a new MC simulation.
    ///
    /// Parameters match Java's `MC(vertices, indices, crystCoord, crystalPixPerUM,
    /// crystSizeVoxels, crystOcc, runNum, verticalGoniometer, xfel, gos,
    /// surrThickness, verticalPolarisation)`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        vertices: Vec<[f64; 3]>,
        indices: Vec<[usize; 3]>,
        cryst_coord: Vec<[f64; 3]>,
        cryst_pix_per_um: f64,
        cryst_size_voxels: [usize; 3],
        cryst_occ: Vec<bool>,
        run_number: usize,
        vertical_goniometer: bool,
        xfel: bool,
        gos: bool,
        surr_thickness: [f64; 3],
        vertical_polarisation: bool,
    ) -> Self {
        let surrounding_thickness = [
            surr_thickness[0] * 1000.0,
            surr_thickness[1] * 1000.0,
            surr_thickness[2] * 1000.0,
        ];

        let x_min_max = Self::min_max_vertices_static(0, &vertices);
        let y_min_max = Self::min_max_vertices_static(1, &vertices);
        let z_min_max = Self::min_max_vertices_static(2, &vertices);
        let x_dimension = 1000.0 * (x_min_max[1] - x_min_max[0]);
        let y_dimension = 1000.0 * (y_min_max[1] - y_min_max[0]);
        let z_dimension = 1000.0 * (z_min_max[1] - z_min_max[0]);

        let n = cryst_size_voxels[0] * cryst_size_voxels[1] * cryst_size_voxels[2];
        let nz = cryst_size_voxels[2];

        let vertices_surrounding = vec![[0.0_f64; 3]; vertices.len()];
        let indices_surrounding = indices.clone();

        MonteCarloSimulation {
            vertices,
            rotated_vertices: Vec::new(),
            expanded_rotated_vertices: Vec::new(),
            indices,
            cryst_coord,
            cryst_pix_per_um,
            cryst_size_voxels,
            cryst_occ,
            run_number,
            normals: None,
            origin_distances: None,
            rotated_normals: None,
            rotated_origin_distances: None,
            vertices_surrounding,
            indices_surrounding,
            normals_surrounding: None,
            origin_distances_surrounding: None,
            rotated_normals_surrounding: None,
            rotated_origin_distances_surrounding: None,
            x_dimension,
            y_dimension,
            z_dimension,
            max_dims: [0.0; 3],
            min_dims: [0.0; 3],
            dose: 0.0,
            photon_dose: 0.0,
            electron_dose: 0.0,
            gos_electron_dose: 0.0,
            photon_dose_v_resolved: 0.0,
            gos_electron_dose_v_resolved: 0.0,
            electron_dose_surrounding: 0.0,
            raddose_style_dose: 0.0,
            raddose_style_dose_compton: 0.0,
            escaped_energy: 0.0,
            angular_emission_probs: vec![0.0; 50],
            number_angular_emission_bins: 50,
            low_energy_angles: vec![None; 95],
            high_energy_angles: vec![None; 95],
            auger_data: HashMap::new(),
            tot_k_auger_prob: HashMap::new(),
            fl_data: HashMap::new(),
            fl_cumulative: HashMap::new(),
            voxel_energy_v_resolved: vec![0.0; n],
            voxel_ionisations_v_resolved: vec![0.0; n],
            voxel_elastic: vec![0.0; n],
            dose_simple: vec![0.0; n],
            last_time_vox: vec![0.0; nz],
            total_ionisation_events: 0,
            total_ionisation_events_v_resolved: 0,
            low_energy_ionisations: 0,
            ionisations_old: 0.0,
            atomic_ionisations: HashMap::new(),
            atomic_ionisations_exposed: HashMap::new(),
            atomic_ionisations_per_atom: HashMap::new(),
            atomic_ionisations_per_atom_exposed: HashMap::new(),
            energy_per_inel: BTreeMap::new(),
            avg_w: 0.0,
            avg_w_num: 0,
            avg_uk: 0.0,
            avg_uk_num: 0,
            vertical_goni: vertical_goniometer,
            vertical_pol: vertical_polarisation,
            simple_mc: !gos,
            do_xfel: xfel,
            tot_elastic: 0.0,
            num_photons: 1_000_000,
            pulse_length: 30.0,
            pulse_bin_length: 0.1,
            pulse_energy: 1.4e-3,
            mean_energy: 12.4,
            energy_fwhm: None,
            photon_energy_array: Vec::new(),
            surrounding_thickness,
            last_time: 0.0,
        }
    }

    // ── Static geometry helpers ──────────────────────────────────────────────

    pub fn min_max_vertices_static(dimension: usize, vertices: &[[f64; 3]]) -> [f64; 2] {
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        for v in vertices {
            let x = v[dimension];
            if x < min {
                min = x;
            }
            if x > max {
                max = x;
            }
        }
        [min, max]
    }

    pub fn min_max_vertices(&self, dimension: usize) -> [f64; 2] {
        Self::min_max_vertices_static(dimension, &self.vertices)
    }

    // ── Normal calculation ────────────────────────────────────────────────────

    fn calculate_normals_for(
        vertices: &[[f64; 3]],
        indices: &[[usize; 3]],
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        let mut normals = Vec::with_capacity(indices.len());
        let mut dists = Vec::with_capacity(indices.len());
        for tri in indices {
            let v0 = vertices[tri[0]];
            let v1 = vertices[tri[1]];
            let v2 = vertices[tri[2]];
            let a = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            let b = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            let n = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let norm = if mag > 0.0 {
                [n[0] / mag, n[1] / mag, n[2] / mag]
            } else {
                [0.0, 0.0, 1.0]
            };
            let d = -(norm[0] * v0[0] + norm[1] * v0[1] + norm[2] * v0[2]);
            normals.push(norm);
            dists.push(d);
        }
        (normals, dists)
    }

    pub fn calculate_normals(&mut self, rotated: bool) {
        if rotated {
            if self.rotated_normals.is_none() {
                let (n, d) = Self::calculate_normals_for(&self.rotated_vertices, &self.indices);
                self.rotated_normals = Some(n);
                self.rotated_origin_distances = Some(d);
            }
        } else if self.normals.is_none() {
            let (n, d) = Self::calculate_normals_for(&self.vertices, &self.indices);
            self.normals = Some(n);
            self.origin_distances = Some(d);
        }
    }

    fn calculate_normals_surrounding(&mut self) {
        if self.normals_surrounding.is_none() {
            let (n, d) =
                Self::calculate_normals_for(&self.vertices_surrounding, &self.indices_surrounding);
            self.normals_surrounding = Some(n);
            self.origin_distances_surrounding = Some(d);
        }
    }

    // ── Point-in-polyhedron (ray-casting along +Z) ────────────────────────────

    fn point_in_polyhedron_static(
        x: f64,
        y: f64,
        z: f64,
        normals: &[[f64; 3]],
        dists: &[f64],
        indices: &[[usize; 3]],
        vertices: &[[f64; 3]],
    ) -> bool {
        let ray_dir = [0.0_f64, 0.0, 1.0];
        let mut count = 0i32;
        for (i, tri) in indices.iter().enumerate() {
            let n = normals[i];
            let d = dists[i];
            let denom = n[0] * ray_dir[0] + n[1] * ray_dir[1] + n[2] * ray_dir[2];
            if denom.abs() < 1e-12 {
                continue;
            }
            let t = -(n[0] * x + n[1] * y + n[2] * z + d) / denom;
            if t < 0.0 {
                continue;
            }
            let hit = [x + t * ray_dir[0], y + t * ray_dir[1], z + t * ray_dir[2]];
            let v0 = vertices[tri[0]];
            let v1 = vertices[tri[1]];
            let v2 = vertices[tri[2]];
            if Self::point_in_triangle_2d(&hit, &v0, &v1, &v2) {
                count += 1;
            }
        }
        count % 2 == 1
    }

    fn point_in_triangle_2d(p: &[f64; 3], a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> bool {
        let px = p[0];
        let py = p[1];
        let ax = a[0];
        let ay = a[1];
        let bx = b[0];
        let by = b[1];
        let cx = c[0];
        let cy = c[1];
        let d1 = (px - bx) * (ay - by) - (ax - bx) * (py - by);
        let d2 = (px - cx) * (by - cy) - (bx - cx) * (py - cy);
        let d3 = (px - ax) * (cy - ay) - (cx - ax) * (py - ay);
        let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
        let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
        !(has_neg && has_pos)
    }

    pub fn calculate_crystal_occupancy(&mut self, x: f64, y: f64, z: f64) -> bool {
        self.calculate_normals(false);
        let normals = self.normals.as_ref().unwrap();
        let dists = self.origin_distances.as_ref().unwrap();
        let vertices = &self.vertices;
        let indices = &self.indices;
        Self::point_in_polyhedron_static(x, y, z, normals, dists, indices, vertices)
    }

    pub fn calculate_crystal_occupancy_surrounding(&mut self, x: f64, y: f64, z: f64) -> bool {
        self.calculate_normals_surrounding();
        let normals = self.normals_surrounding.as_ref().unwrap();
        let dists = self.origin_distances_surrounding.as_ref().unwrap();
        let vertices = &self.vertices_surrounding;
        let indices = &self.indices_surrounding;
        Self::point_in_polyhedron_static(x, y, z, normals, dists, indices, vertices)
    }

    // ── Surrounding volume setup ──────────────────────────────────────────────

    fn setup_surrounding_volume(&mut self) {
        let verts: Vec<[f64; 3]> = self
            .vertices
            .iter()
            .map(|v| {
                let mut sv = *v;
                for (j, component) in sv.iter_mut().enumerate() {
                    if *component < 0.0 {
                        *component -= self.surrounding_thickness[j] / 1000.0;
                    } else {
                        *component += self.surrounding_thickness[j] / 1000.0;
                    }
                }
                sv
            })
            .collect();
        self.vertices_surrounding = verts;
        self.normals_surrounding = None;
        self.origin_distances_surrounding = None;
    }

    // ── Vertex rotation for wedge angle ──────────────────────────────────────

    pub fn set_up_rotated_vertices(&mut self, ang_rad: f64, _wedge: &Wedge) {
        let cos_a = ang_rad.cos();
        let sin_a = ang_rad.sin();
        self.rotated_vertices = self
            .vertices
            .iter()
            .map(|v| {
                [
                    v[0] * cos_a + v[2] * sin_a,
                    v[1],
                    -v[0] * sin_a + v[2] * cos_a,
                ]
            })
            .collect();
        self.expanded_rotated_vertices = self
            .indices
            .iter()
            .map(|tri| {
                [
                    self.rotated_vertices[tri[0]],
                    self.rotated_vertices[tri[1]],
                    self.rotated_vertices[tri[2]],
                ]
            })
            .collect();
        self.rotated_normals = None;
        self.rotated_origin_distances = None;
    }

    fn set_up_rotated_vertices_surrounding(&mut self, ang_rad: f64, _wedge: &Wedge) {
        let cos_a = ang_rad.cos();
        let sin_a = ang_rad.sin();
        // rotate surrounding vertices
        let rotated: Vec<[f64; 3]> = self
            .vertices_surrounding
            .iter()
            .map(|v| {
                [
                    v[0] * cos_a + v[2] * sin_a,
                    v[1],
                    -v[0] * sin_a + v[2] * cos_a,
                ]
            })
            .collect();
        // Update the surrounding normals when we set new vertices
        let (n, d) = Self::calculate_normals_for(&rotated, &self.indices_surrounding);
        self.rotated_normals_surrounding = Some(n);
        self.rotated_origin_distances_surrounding = Some(d);
    }

    // ── Coordinate helpers ────────────────────────────────────────────────────

    /// Check if a point (nm) is inside the crystal at this rotation angle.
    fn is_microcrystal_at(
        &mut self,
        x_nm: f64,
        y_nm: f64,
        z_nm: f64,
        angle: f64,
        _wedge: &Wedge,
    ) -> bool {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        // Convert nm → µm and unrotate
        let x = x_nm / 1000.0;
        let y = y_nm / 1000.0;
        let z = z_nm / 1000.0;
        let ux = x * cos_a - z * sin_a;
        let uy = y;
        let uz = x * sin_a + z * cos_a;
        self.calculate_crystal_occupancy(ux, uy, uz)
    }

    /// Convert nm position to voxel coordinates.
    fn convert_to_pixel_coordinates(
        &self,
        x_nm: f64,
        y_nm: f64,
        z_nm: f64,
        angle: f64,
        _wedge: &Wedge,
    ) -> [usize; 3] {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let x = x_nm / 1000.0;
        let y = y_nm / 1000.0;
        let z = z_nm / 1000.0;
        let ux = x * cos_a - z * sin_a;
        let uy = y;
        let uz = x * sin_a + z * cos_a;
        let pix = self.cryst_pix_per_um;
        let [nx, ny, nz] = self.cryst_size_voxels;
        let i = ((ux * pix) + (nx as f64 / 2.0)).clamp(0.0, nx as f64 - 1.0) as usize;
        let j = ((uy * pix) + (ny as f64 / 2.0)).clamp(0.0, ny as f64 - 1.0) as usize;
        let k = ((uz * pix) + (nz as f64 / 2.0)).clamp(0.0, nz as f64 - 1.0) as usize;
        [i, j, k]
    }

    /// Convert pixel coords back to Cartesian nm.
    fn convert_to_cartesian_coordinates(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let [nx, ny, nz] = self.cryst_size_voxels;
        let pix = self.cryst_pix_per_um;
        let x = (i as f64 - nx as f64 / 2.0) / pix * 1000.0;
        let y = (j as f64 - ny as f64 / 2.0) / pix * 1000.0;
        let z = (k as f64 - nz as f64 / 2.0) / pix * 1000.0;
        [x, y, z]
    }

    // ── Starting Z for photon injection ─────────────────────────────────────

    fn get_starting_z(&mut self, _angle: f64, _wedge: &Wedge, xy: &[f64; 3], front: bool) -> f64 {
        // Ray-trace from (x_um, y_um, 0) in +z direction through the rotated crystal,
        // matching Java's MC.getStartingZ ray-tracing approach.
        // xy is in µm; returns z in nm.
        self.calculate_normals(true);
        let normals = self.rotated_normals.as_ref().unwrap().clone();
        let dists = self.rotated_origin_distances.as_ref().unwrap().clone();
        let tris = self.expanded_rotated_vertices.clone();

        let x = xy[0]; // µm
        let y = xy[1]; // µm

        let mut ts: Vec<f64> = Vec::new();
        for (i, tri) in tris.iter().enumerate() {
            let n = normals[i];
            let d = dists[i];
            let denom = n[2]; // ray direction is (0,0,1)
            if denom.abs() < 1e-12 {
                continue;
            }
            let t = -(n[0] * x + n[1] * y + d) / denom;
            let hit = [x, y, t];
            if Self::point_in_triangle_2d(&hit, &tri[0], &tri[1], &tri[2]) {
                ts.push(t);
            }
        }

        if ts.is_empty() {
            return 0.0;
        }
        ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // Deduplicate near-equal t values (shared triangle edges)
        ts.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

        let idx = if front { 0 } else { ts.len().saturating_sub(1) };
        ts.get(idx).copied().unwrap_or(0.0) * 1000.0 // µm → nm
    }

    fn get_starting_z_surrounding(&mut self, _angle: f64, _wedge: &Wedge, xy: &[f64; 3]) -> f64 {
        let mut z = -(self.z_dimension / 2.0 + self.surrounding_thickness[2] + 50.0);
        for _ in 0..6000 {
            let xnm = xy[0] * 1000.0;
            let ynm = xy[1] * 1000.0;
            if self.calculate_crystal_occupancy_surrounding(xnm / 1000.0, ynm / 1000.0, z / 1000.0)
            {
                return z;
            }
            z += 1.0;
        }
        0.0
    }

    // ── Intersection distance ─────────────────────────────────────────────────

    #[allow(clippy::too_many_arguments)]
    fn get_intersection_distance(
        &mut self,
        px: f64,
        py: f64,
        pz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        from_crystal: bool,
        angle: f64,
        _wedge: &Wedge,
    ) -> f64 {
        // March along direction until we leave the crystal, return distance in µm.
        let step_nm = 1.0_f64;
        let mut dist = 0.0;
        let max_dist = 100_000.0; // 100 µm max
        loop {
            dist += step_nm;
            let x = px + dist * dx;
            let y = py + dist * dy;
            let z = pz + dist * dz;
            let inside = self.is_microcrystal_at(x, y, z, angle, _wedge);
            if from_crystal && !inside {
                break;
            }
            if !from_crystal && inside {
                break;
            }
            if dist > max_dist {
                break;
            }
        }
        dist / 1000.0 // µm
    }

    fn get_intersection_point(
        dist_um: f64,
        px: f64,
        py: f64,
        pz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> [f64; 3] {
        let dist_nm = dist_um * 1000.0;
        [px + dist_nm * dx, py + dist_nm * dy, pz + dist_nm * dz]
    }

    // ── Sample volume ─────────────────────────────────────────────────────────

    fn get_sample_volume(&self) -> f64 {
        // Return volume in nm^3
        self.x_dimension * self.y_dimension * self.z_dimension
    }

    // ── Beam area / exposure ──────────────────────────────────────────────────

    pub fn get_exposed_area(&self, beam: &dyn Beam) -> f64 {
        // Returns area in nm²
        let bx = beam.beam_size_x(); // µm
        let by = beam.beam_size_y(); // µm
        if beam.is_circular() {
            PI * (bx * 1000.0) * (by * 1000.0)
        } else {
            (bx * 1000.0) * (by * 1000.0)
        }
    }

    fn test_if_inside_exposed_area(&self, x_nm: f64, y_nm: f64, beam: &dyn Beam) -> bool {
        let bx = beam.beam_size_x();
        let by = beam.beam_size_y();
        if beam.is_circular() {
            let ex = (x_nm / 1000.0) / (bx / 2.0);
            let ey = (y_nm / 1000.0) / (by / 2.0);
            ex * ex + ey * ey <= 1.0
        } else {
            (x_nm / 1000.0).abs() <= bx / 2.0 && (y_nm / 1000.0).abs() <= by / 2.0
        }
    }

    /// Sample a random (x, y) photon start position from the beam profile.
    fn get_photon_beam_xy_pos(&self, beam: &dyn Beam, _wedge: &Wedge) -> [f64; 2] {
        let bx = beam.beam_size_x();
        let by = beam.beam_size_y();
        let rx: f64 = rand::random();
        let ry: f64 = rand::random();
        let x_nm = 1000.0 * bx * (rx - 0.5);
        let y_nm = 1000.0 * ((ry * by) - by / 2.0);
        [x_nm, y_nm]
    }

    // ── Auger + FL linewidths data loading ────────────────────────────────────

    pub fn populate_auger_linewidths(&mut self) {
        for (i, &z) in AUGER_ELEMENTS.iter().enumerate() {
            let num_shells = get_num_auger_shells(i);
            let mut element_map: ShellTransitionMap = HashMap::new();
            for shell in 0..=num_shells {
                let filename = format!("{}-{}.csv", z, shell);
                let data = load_auger_fl_csv("auger_linewidths", &filename);
                if let Some(td) = data {
                    if shell == 0 {
                        self.tot_k_auger_prob
                            .insert(z, td.probs.iter().sum::<f64>());
                    }
                    element_map.insert(shell, td);
                }
            }
            self.auger_data.insert(z, element_map);
        }
    }

    pub fn populate_fluorescence_linewidths(&mut self) {
        for (i, &z) in AUGER_ELEMENTS.iter().enumerate() {
            let num_shells = get_num_auger_shells(i);
            let mut element_map: ShellTransitionMap = HashMap::new();
            let mut cumulative_map: HashMap<u32, Vec<f64>> = HashMap::new();
            for shell in 0..=num_shells {
                let filename = format!("{}-{}.csv", z, shell);
                let data = load_fl_csv("fl_linewidths", &filename);
                if let Some(td) = data {
                    cumulative_map.insert(shell, td.cumulative_probs.clone());
                    element_map.insert(shell, td);
                }
            }
            self.fl_data.insert(z, element_map);
            self.fl_cumulative.insert(z, cumulative_map);
        }
    }

    // ── Angular emission probabilities ────────────────────────────────────────

    pub fn populate_angular_emission_probs(&mut self) {
        let n = self.number_angular_emission_bins;
        self.angular_emission_probs = vec![0.0; n];
        // Integrate under polarisation equation
        let mut total_area = 0.0;
        let mut last_height = 0.0;
        for i in 0..=100 {
            let angle = (PI / 100.0) * i as f64;
            let height = Self::solve_polarisation_equation_for_angle(angle, 1.0, 2.0);
            if i > 0 {
                total_area += (last_height + height) / 2.0 * (PI / 100.0);
            }
            last_height = height;
        }
        let mut cum_prob = 0.0;
        last_height = 0.0;
        for i in 0..=n {
            let angle = (PI / n as f64) * i as f64;
            let height = Self::solve_polarisation_equation_for_angle(angle, 1.0, 2.0);
            if i > 0 {
                let area = (last_height + height) / 2.0 * (PI / n as f64);
                cum_prob += area / total_area;
                self.angular_emission_probs[i - 1] = cum_prob;
            }
            last_height = height;
        }
    }

    fn solve_polarisation_equation_for_angle(phi: f64, photo_electric: f64, beta: f64) -> f64 {
        (photo_electric / (4.0 * PI)) * (1.0 + beta * 0.5 * (3.0 * phi.cos().powi(2) - 1.0))
    }

    fn get_cos_angle_to_x(&self) -> f64 {
        let rnd: f64 = rand::random();
        let n = self.number_angular_emission_bins;
        let mut last_prob = 0.0;
        let mut angle = 0.0;
        for i in 0..n {
            if rnd < self.angular_emission_probs[i] {
                let start = i as f64 * PI / n as f64;
                let end = (i + 1) as f64 * PI / n as f64;
                let prop = (rnd - last_prob) / (self.angular_emission_probs[i] - last_prob);
                angle = start + prop * (end - start);
                break;
            }
            last_prob = self.angular_emission_probs[i];
        }
        angle.cos()
    }

    fn pos_or_neg() -> f64 {
        let r: f64 = rand::random();
        if r < 0.5 {
            1.0
        } else {
            -1.0
        }
    }

    // ── Elastic scattering angle (ELSEPA tables) ──────────────────────────────

    fn get_primary_elastic_scattering_angle(
        &mut self,
        electron_energy: f64,
        atomic_number: usize,
    ) -> f64 {
        let high_energy = electron_energy > 20.0;
        if atomic_number >= 95 {
            return 0.0;
        }

        // Load table if needed
        if high_energy {
            if self.high_energy_angles[atomic_number].is_none() {
                self.high_energy_angles[atomic_number] = load_angle_file(true, atomic_number);
            }
        } else if self.low_energy_angles[atomic_number].is_none() {
            self.low_energy_angles[atomic_number] = load_angle_file(false, atomic_number);
        }

        let table = if high_energy {
            self.high_energy_angles[atomic_number].as_ref()
        } else {
            self.low_energy_angles[atomic_number].as_ref()
        };

        if let Some(map) = table {
            if map.is_empty() {
                return 0.0;
            }
            let energy_key = Self::nearest_energy_key(map, electron_energy);
            if let Some(probs) = map.get(&energy_key) {
                return Self::return_deflection_angle(high_energy, probs);
            }
        }
        0.0
    }

    fn nearest_energy_key(map: &BTreeMap<u64, Vec<f64>>, energy: f64) -> u64 {
        let before = map.range(..=f64_key(energy)).next_back().map(|(k, _)| *k);
        let after = map.range(f64_key(energy)..).next().map(|(k, _)| *k);
        match (before, after) {
            (Some(b), Some(a)) => {
                let be = f64::from_bits(b);
                let ae = f64::from_bits(a);
                if (energy - be).abs() <= (ae - energy).abs() {
                    b
                } else {
                    a
                }
            }
            (Some(b), None) => b,
            (None, Some(a)) => a,
            (None, None) => f64_key(energy),
        }
    }

    fn return_deflection_angle(high_energy: bool, probs: &[f64]) -> f64 {
        let total: f64 = probs.iter().sum();
        if total == 0.0 {
            return 0.0;
        }
        let mut cum = 0.0;
        let rnd: f64 = rand::random();
        let mut index = 0usize;
        for (k, &p) in probs.iter().enumerate() {
            cum += p / total;
            if cum >= rnd {
                index = k;
                break;
            }
        }
        let angle_degrees = if high_energy {
            if index == 0 {
                0.0
            } else if index < 146 {
                let factor = (index as f64 - 1.0) / 36.0;
                let start = 10f64.powf(factor.floor()) * 0.0001;
                start + ((index - 1) % 36) as f64 * (start / 4.0)
            } else if index < 236 {
                1.0 + (index - 146) as f64 / 10.0
            } else if index < 297 {
                10.0 + (index - 236) as f64 / 4.0 * (10.0 / 4.0)
            } else {
                25.0 + (index - 296) as f64
            }
        } else {
            index as f64
        };
        angle_degrees * PI / 180.0
    }

    // ── FSE cross section ─────────────────────────────────────────────────────

    fn get_fse_x_section(electron_energy: f64) -> f64 {
        let e_charge = 4.803_204_25e-10_f64;
        let m_cgs = 9.109_383_56e-28_f64;
        let c_cgs = 2.997_924_58e10_f64;
        let c2 = (c_cgs / 100.0).powi(2);
        let vo = electron_energy * KEV_TO_JOULES;
        let beta2 = 1.0 - (m_cgs / 1000.0 * c2 / (vo + m_cgs / 1000.0 * c2)).powi(2);
        let v2 = beta2 * c2 * 10_000.0;
        let constant = (2.0 * PI * e_charge.powi(4)) / (m_cgs * v2 * (vo * 1000.0 * 10_000.0));
        let energy_cutoff = (14.0 / 1000.0) / electron_energy;
        let tau = electron_energy / 511.0;
        let cs = (((2.0 * tau + 1.0) / (tau + 1.0).powi(2)) * ((1.0_f64 / 0.5 - 1.0).ln())
            + (tau / (tau + 1.0)).powi(2)
            - 1.0 / 0.5
            - 1.0 / (0.5 - 1.0))
            - (((2.0 * tau + 1.0) / (tau + 1.0).powi(2)) * ((1.0 / energy_cutoff - 1.0).ln())
                + (tau / (tau + 1.0)).powi(2)
                - 1.0 / energy_cutoff
                - 1.0 / (energy_cutoff - 1.0));
        cs * constant
    }

    fn get_time_to_distance(electron_energy: f64, s_nm: f64) -> f64 {
        let c2 = C * C;
        let vo = electron_energy * KEV_TO_JOULES;
        let beta2 = 1.0 - (M_E * c2 / (vo + M_E * c2)).powi(2);
        let v = (beta2 * c2).sqrt() * 1.0e9 / 1.0e15; // nm/fs
        s_nm / v
    }

    // ── Physics: GOS model ────────────────────────────────────────────────────

    pub fn wk_to_wak(e: f64, wk: f64, uk: f64) -> f64 {
        if e * 1000.0 > 3.0 * wk - 2.0 * uk {
            wk
        } else {
            (e * 1000.0 + 2.0 * uk) / 3.0
        }
    }

    pub fn get_qak(e: f64, wk: f64, uk: f64) -> f64 {
        if e * 1000.0 > 3.0 * wk - 2.0 * uk {
            uk
        } else {
            uk * (e * 1000.0 / (3.0 * wk - 2.0 * uk))
        }
    }

    pub fn get_energy_loss_distant(wdis: f64, uk: f64) -> f64 {
        let rnd: f64 = rand::random();
        wdis - (rnd * (wdis - uk).powi(2)).sqrt()
    }

    fn get_gos_primary_scatter_long(e_kev: f64, q: f64, _wak_ev: f64) -> f64 {
        // Longitudinal GOS deflection angle
        let mc2 = M_E * C * C;
        let e_j = e_kev * KEV_TO_JOULES;
        let _t = e_j / mc2;
        if q == 0.0 {
            return 0.0;
        }

        (q / (e_j + mc2)).asin().abs()
    }

    fn get_gos_primary_scatter_close(e_kev: f64, w_kev: f64) -> f64 {
        let mc2 = M_E * C * C;
        let e_j = e_kev * KEV_TO_JOULES;
        let _t = e_j / mc2;
        let w_j = w_kev * KEV_TO_JOULES;
        let theta = (w_j / (2.0 * e_j + mc2)).asin().abs();
        theta.min(PI)
    }

    fn secondary_scatter_distant(e_kev: f64, wak_ev: f64, _q: f64) -> f64 {
        let mc2 = M_E * C * C;
        let e_j = e_kev * KEV_TO_JOULES;
        let w_j = wak_ev / 1000.0 * KEV_TO_JOULES;
        let theta = (w_j / (2.0 * e_j + mc2)).asin().abs();
        theta.min(PI)
    }

    fn secondary_scatter_close(e_kev: f64, w_kev: f64) -> f64 {
        let mc2 = M_E * C * C;
        let e_j = e_kev * KEV_TO_JOULES;
        let w_j = w_kev * KEV_TO_JOULES;
        let theta = (w_j / (e_j + mc2)).asin().abs();
        theta.min(PI)
    }

    pub fn sample_k(_e: f64, qk: f64) -> f64 {
        let rnd: f64 = rand::random();
        // Simple approximation: k proportional to qk
        qk * (1.0 + rnd)
    }

    fn get_new_direction_vector(
        xn: f64,
        yn: f64,
        zn: f64,
        scatter_theta: f64,
        scatter_phi: f64,
    ) -> [f64; 3] {
        if yn.abs() < 1e-10 {
            let xn2 = xn * scatter_theta.cos() + zn * scatter_phi.sin() * scatter_theta.sin();
            let yn2 = yn * scatter_theta.cos() + scatter_phi.cos() * scatter_theta.sin();
            let zn2 = -xn * scatter_theta.sin() * scatter_phi.sin() + zn * scatter_theta.cos();
            let mag = (xn2 * xn2 + yn2 * yn2 + zn2 * zn2).sqrt().max(1e-15);
            return [xn2 / mag, yn2 / mag, zn2 / mag];
        }
        let an = -(xn / yn);
        let am = 1.0 / (1.0 + an * an).sqrt();
        let v1 = an * scatter_theta.sin();
        let v2 = an * am * scatter_theta.sin();
        let v3 = scatter_phi.cos();
        let v4 = scatter_phi.sin();
        let ca = xn * scatter_theta.cos() + v1 * v3 + yn * v2 * v4;
        let cb = yn * scatter_theta.cos() + v4 * (zn * v1 - xn * v2);
        let cc = zn * scatter_theta.cos() + v2 * v3 - yn * v1 * v4;
        let mag = (ca * ca + cb * cb + cc * cc).sqrt().max(1e-15);
        [ca / mag, cb / mag, cc / mag]
    }

    fn get_scattering_theta(
        &mut self,
        electron_energy: f64,
        elastic_probs: &HashMap<String, f64>,
    ) -> f64 {
        let z = self.get_elastic_element_z(elastic_probs);
        self.get_primary_elastic_scattering_angle(electron_energy, z as usize)
    }

    fn get_scattering_phi() -> f64 {
        let r: f64 = rand::random();
        2.0 * PI * r
    }

    fn get_elastic_element_z(&self, probs: &HashMap<String, f64>) -> u32 {
        let rnd: f64 = rand::random();
        for (name, &p) in probs {
            if p > rnd {
                // Try to parse name as atomic number
                if let Ok(z) = name.parse::<u32>() {
                    return z;
                }
                return 6; // default to carbon
            }
        }
        6
    }

    // ── Shell binding energy ──────────────────────────────────────────────────

    fn get_shell_binding_energy(element: &Element, shell: usize) -> f64 {
        match shell {
            0 => element.k_edge().unwrap_or(0.0),
            1 => element.l1_edge().unwrap_or(0.0),
            2 => element.l2_edge().unwrap_or(0.0),
            3 => element.l3_edge().unwrap_or(0.0),
            4 => element.m1_edge(),
            5 => element.m2_edge().unwrap_or(0.0),
            6 => element.m3_edge().unwrap_or(0.0),
            7 => element.m4_edge().unwrap_or(0.0),
            8 => element.m5_edge().unwrap_or(0.0),
            _ => 0.0,
        }
    }

    fn get_shell_fluorescence_yield(element: &Element, shell: usize) -> f64 {
        match shell {
            0 => element.k_fluorescence_yield().unwrap_or(0.0),
            1 => element.l1_fluorescence_yield().unwrap_or(0.0),
            2 => element.l2_fluorescence_yield().unwrap_or(0.0),
            3 => element.l3_fluorescence_yield().unwrap_or(0.0),
            4 => element.m1_fluorescence_yield().unwrap_or(0.0),
            5 => element.m2_fluorescence_yield().unwrap_or(0.0),
            6 => element.m3_fluorescence_yield().unwrap_or(0.0),
            7 => element.m4_fluorescence_yield().unwrap_or(0.0),
            8 => element.m5_fluorescence_yield().unwrap_or(0.0),
            _ => 0.0,
        }
    }

    // ── Auger transition index ────────────────────────────────────────────────

    fn get_transition_index_auger(&self, z: u32, shell: u32) -> usize {
        if let Some(em) = self.auger_data.get(&z) {
            if let Some(td) = em.get(&shell) {
                let rnd: f64 = rand::random();
                for (i, &p) in td.cumulative_probs.iter().enumerate() {
                    if rnd < p {
                        return i;
                    }
                }
                return td.cumulative_probs.len().saturating_sub(1);
            }
        }
        0
    }

    fn get_transition_index_fl(&self, z: u32, shell: u32) -> usize {
        if let Some(em) = self.fl_data.get(&z) {
            if let Some(td) = em.get(&shell) {
                let rnd: f64 = rand::random();
                for (i, &p) in td.cumulative_probs.iter().enumerate() {
                    if rnd < p {
                        return i;
                    }
                }
                return td.cumulative_probs.len().saturating_sub(1);
            }
        }
        0
    }

    // ── Relative shell probabilities ──────────────────────────────────────────

    fn get_relative_shell_probs(
        element_abs_probs: &HashMap<String, f64>,
        beam_energy: f64,
        elements: &HashMap<String, Element>,
    ) -> HashMap<String, Vec<f64>> {
        let mut result = HashMap::new();
        for name in element_abs_probs.keys() {
            if let Some(e) = elements.get(name) {
                let mut shell_probs = vec![0.0_f64; 9];
                let mut running = 0.0;
                let k = e.k_ionisation_prob();
                if beam_energy > e.k_edge().unwrap_or(f64::INFINITY) {
                    running += k;
                    shell_probs[0] = running;
                }
                if beam_energy > e.l1_edge().unwrap_or(f64::INFINITY) && e.atomic_number() >= 12 {
                    let p = e.l1_ionisation_prob() * (1.0 - k);
                    running += p;
                    shell_probs[1] = running;
                }
                if beam_energy > e.l2_edge().unwrap_or(f64::INFINITY) && e.atomic_number() >= 12 {
                    let p = e.l2_ionisation_prob() * (1.0 - shell_probs[0] - shell_probs[1]);
                    running += p;
                    shell_probs[2] = running;
                }
                if beam_energy > e.l3_edge().unwrap_or(f64::INFINITY) && e.atomic_number() >= 12 {
                    let p = e.l3_ionisation_prob()
                        * (1.0 - shell_probs[0] - shell_probs[1] - shell_probs[2]);
                    running += p;
                    shell_probs[3] = running;
                }
                result.insert(name.clone(), shell_probs);
            }
        }
        result
    }

    // ── Ionisation helpers ────────────────────────────────────────────────────

    fn get_ionised_element(element_probs: &HashMap<String, f64>) -> Option<String> {
        let rnd: f64 = rand::random();
        // Sort by cumulative probability (ascending) so the sample is correct
        // regardless of HashMap iteration order.
        let mut sorted: Vec<(&String, f64)> = element_probs.iter().map(|(k, &v)| (k, v)).collect();
        sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (name, p) in sorted {
            if p > rnd {
                return Some(name.clone());
            }
        }
        None
    }

    fn get_ionised_shell(shell_probs: &[f64]) -> usize {
        let rnd: f64 = rand::random();
        for (i, &p) in shell_probs.iter().enumerate() {
            if p > rnd {
                return i;
            }
        }
        0
    }

    fn find_if_element_ionised(shell_probs: &[f64], rnd: f64) -> i32 {
        for (k, &p) in shell_probs.iter().enumerate() {
            if p > rnd {
                return k as i32;
            }
        }
        -1
    }

    fn add_ionisation_vox(
        &mut self,
        pixel_coord: [usize; 3],
        element_name: &str,
        x_nm: f64,
        y_nm: f64,
        beam: &dyn Beam,
        count: u64,
    ) {
        let idx = vox3(
            &self.cryst_size_voxels,
            pixel_coord[0],
            pixel_coord[1],
            pixel_coord[2],
        );
        self.total_ionisation_events += count;
        self.total_ionisation_events_v_resolved += count;
        if idx < self.voxel_ionisations_v_resolved.len() {
            self.voxel_ionisations_v_resolved[idx] += count as f64;
        }
        *self
            .atomic_ionisations
            .entry(element_name.to_string())
            .or_insert(0) += count;
        if self.test_if_inside_exposed_area(x_nm, y_nm, beam) {
            *self
                .atomic_ionisations_exposed
                .entry(element_name.to_string())
                .or_insert(0) += count;
        }
    }

    // ── GOS collision type ────────────────────────────────────────────────────

    fn get_gos_inelastic_type(shell_probs: &[f64; 4], _shell_index: usize) -> usize {
        if shell_probs[3] == 0.0 {
            return 0;
        }
        let rnd: f64 = rand::random();
        let mut running = 0.0;
        for i in 0..3 {
            running += shell_probs[i] / shell_probs[3];
            if running > rnd {
                return i;
            }
        }
        0
    }

    // ── Core photon loop ──────────────────────────────────────────────────────

    /// Main Monte Carlo loop.
    pub fn start_monte_carlo_xfel(
        &mut self,
        beam: &dyn Beam,
        wedge: &Wedge,
        coef_calc: &mut dyn CoefCalc,
    ) {
        self.populate_auger_linewidths();
        self.populate_fluorescence_linewidths();
        self.populate_angular_emission_probs();
        coef_calc.populate_cross_section_coefficients();

        let num_photons = self.num_photons;
        let photon_divisions = num_photons as f64 / (self.pulse_length / self.pulse_bin_length);
        let mut miss_count = 0u64;
        let mut last_progress = 0.0_f64;

        for i in 0..num_photons {
            let energy_of_photon = if (i as usize) < self.photon_energy_array.len() {
                self.photon_energy_array[i as usize]
            } else {
                self.mean_energy
            };

            let progress = i as f64 / num_photons as f64;
            if progress - last_progress >= 0.05 {
                last_progress = progress;
                eprint!("{}% ", (progress * 100.0) as u32);
            }

            // Determine rotation angle for this photon
            let angle = if (wedge.start_ang - wedge.end_ang).abs() < wedge.ang_res
                || wedge.end_ang == 0.0
            {
                0.0
            } else {
                let sign = if wedge.end_ang < wedge.start_ang {
                    -1.0
                } else {
                    1.0
                };
                let a = wedge.start_ang + progress * sign * (wedge.end_ang - wedge.start_ang);
                let turns = (a / (2.0 * PI)).floor();
                a - turns * 2.0 * PI
            };

            self.set_up_rotated_vertices(angle, wedge);

            // Update coefficients for this photon energy
            coef_calc.update_coefficients(energy_of_photon);
            let abs_coef = coef_calc.absorption_coefficient();
            let compton_coef = coef_calc.inelastic_coefficient();
            let elastic_coef = coef_calc.elastic_coefficient();

            let photon_mfpl = (1.0 / (abs_coef + compton_coef)) * 1000.0; // nm
            let prob_compton = 1.0 - photon_mfpl / ((1.0 / abs_coef) * 1000.0);
            let total_mfpl = (1.0 / (abs_coef + compton_coef + elastic_coef)) * 1000.0;
            let elastic_prob = elastic_coef / (abs_coef + compton_coef);

            let element_abs_probs = coef_calc.photo_electric_probs_element(energy_of_photon);
            let element_compton_probs = coef_calc.compton_probs_element(energy_of_photon);
            let ionisation_probs = coef_calc.relative_shell_probs(energy_of_photon, false);

            // Time stamp for this photon
            let time_stamp = (i as f64 / photon_divisions) * self.pulse_bin_length;

            let xy_pos = self.get_photon_beam_xy_pos(beam, wedge);
            let previous_x = xy_pos[0];
            let previous_y = xy_pos[1];

            let coord = [previous_x / 1000.0, previous_y / 1000.0, 0.0];
            let previous_z = self.get_starting_z(angle, wedge, &coord, true);

            if previous_z == 0.0 {
                miss_count += 1;
                continue;
            }

            // Sample photon free path in crystal
            let rnd_s: f64 = rand::random();
            let s = -total_mfpl * rnd_s.ln();
            let xn = previous_x + s * 0.0;
            let yn = previous_y + s * 0.0;
            let zn = previous_z + s;

            // Check if it hit the crystal
            if self.is_microcrystal_at(xn, yn, zn, angle, wedge) {
                let time_to_point = (1.0 / C) * (s / 1.0e9);
                let t_stamp = time_stamp + time_to_point * 1.0e15;
                let dose_time = (t_stamp / self.pulse_bin_length).max(0.0) as usize;

                let rnd_el: f64 = rand::random();
                if rnd_el < elastic_prob {
                    let pixel = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
                    let idx = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);
                    if idx < self.voxel_elastic.len() {
                        self.voxel_elastic[idx] += 1.0;
                    }
                    self.tot_elastic += 1.0;
                } else {
                    let rnd_c: f64 = rand::random();
                    if rnd_c < prob_compton {
                        self.ionisations_old += 1.0;
                        self.produce_compton(
                            beam,
                            coef_calc,
                            t_stamp,
                            xn,
                            yn,
                            zn,
                            false,
                            energy_of_photon,
                            &element_compton_probs,
                            angle,
                            wedge,
                        );
                    } else {
                        self.ionisations_old += 1.0;
                        let ion_probs_clone = ionisation_probs.clone();
                        let el_probs_clone = element_abs_probs.clone();
                        self.produce_photo_electron(
                            beam,
                            coef_calc,
                            &el_probs_clone,
                            &ion_probs_clone,
                            t_stamp,
                            dose_time,
                            xn,
                            yn,
                            zn,
                            false,
                            energy_of_photon,
                            angle,
                            wedge,
                        );
                    }
                }
            } else {
                miss_count += 1;
            }
        }
        println!();

        self.last_time = (1.0 / C) * (self.z_dimension / 1.0e9) * 1.0e15 + self.pulse_length;
        for k in 0..self.last_time_vox.len() {
            let nz = self.last_time_vox.len();
            self.last_time_vox[k] =
                (1.0 / C) * ((self.z_dimension / 1.0e9) / nz as f64) * (k + 1) as f64 * 1.0e15
                    + self.pulse_length;
        }
    }

    // ── Compton production ────────────────────────────────────────────────────

    #[allow(clippy::too_many_arguments)]
    fn produce_compton(
        &mut self,
        beam: &dyn Beam,
        coef_calc: &mut dyn CoefCalc,
        time_stamp: f64,
        xn: f64,
        yn: f64,
        zn: f64,
        surrounding: bool,
        photon_energy: f64,
        element_compton_probs: &HashMap<String, f64>,
        angle: f64,
        wedge: &Wedge,
    ) {
        let pixel_coord = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
        let photon_theta: f64 = PI * rand::random::<f64>();
        let mc_squared = M_E * C * C;
        let incident_energy = photon_energy * KEV_TO_JOULES;
        let e_comp = ((incident_energy.powi(2) * (1.0 - photon_theta.cos()))
            / (mc_squared * (1.0 + (incident_energy / mc_squared) * (1.0 - photon_theta.cos()))))
            / KEV_TO_JOULES;

        if !surrounding {
            self.raddose_style_dose_compton += e_comp;
        }

        let electron_phi =
            (1.0 / (photon_theta / 2.0).tan() / (1.0 + incident_energy / mc_squared)).atan();
        let z_norm = electron_phi.cos();
        let x_norm = Self::pos_or_neg() * rand::random::<f64>() * (1.0 - z_norm.powi(2)).sqrt();
        let y_norm = Self::pos_or_neg() * (1.0 - x_norm.powi(2) - z_norm.powi(2)).max(0.0).sqrt();
        let theta = z_norm.acos();
        let phi = if theta.sin().abs() > 1e-10 {
            (x_norm / theta.sin()).acos()
        } else {
            0.0
        };

        if !surrounding {
            if let Some(elem_name) = Self::get_ionised_element(element_compton_probs) {
                self.add_ionisation_vox(pixel_coord, &elem_name, xn, yn, beam, 1);
            }
        }

        self.track_photoelectron(
            coef_calc,
            time_stamp,
            e_comp,
            xn,
            yn,
            zn,
            x_norm,
            y_norm,
            z_norm,
            theta,
            phi,
            surrounding,
            true,
            beam,
            angle,
            wedge,
        );
    }

    // ── Photoelectron production ──────────────────────────────────────────────

    #[allow(clippy::too_many_arguments)]
    fn produce_photo_electron(
        &mut self,
        beam: &dyn Beam,
        coef_calc: &mut dyn CoefCalc,
        element_abs_probs: &HashMap<String, f64>,
        ionisation_probs: &HashMap<String, Vec<f64>>,
        time_stamp: f64,
        _dose_time: usize,
        xn: f64,
        yn: f64,
        zn: f64,
        surrounding: bool,
        photon_energy: f64,
        angle: f64,
        wedge: &Wedge,
    ) {
        let pixel_coord = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
        let ionised_elem_name = match Self::get_ionised_element(element_abs_probs) {
            Some(n) => n,
            None => return,
        };

        if !surrounding {
            self.total_ionisation_events += 1;
            self.add_ionisation_vox(pixel_coord, &ionised_elem_name, xn, yn, beam, 1);
        }

        // Get shell ionisation probabilities for this element
        let shell_probs = match ionisation_probs.get(&ionised_elem_name) {
            Some(p) => p.clone(),
            None => return,
        };
        let shell_index = Self::get_ionised_shell(&shell_probs);

        // We need Element for binding energy - use coefcalc
        let shell_binding_energy = coef_calc.shell_binding_energy(&ionised_elem_name, shell_index);
        let photoelectron_energy = (photon_energy - shell_binding_energy).max(0.0);

        if !surrounding {
            self.photon_dose += shell_binding_energy;
            self.raddose_style_dose += photon_energy;
        }

        // Sample photoelectron direction
        let this_angle = 2.0 * PI - angle;
        let polarised: f64 = rand::random();
        let (x_norm, y_norm, z_norm, theta, phi) = if shell_index == 0 && polarised > 0.25 {
            let (xni, zni) = if (self.vertical_goni && !self.vertical_pol)
                || (!self.vertical_goni && self.vertical_pol)
            {
                let xc = self.get_cos_angle_to_x();
                let y = Self::pos_or_neg() * rand::random::<f64>() * (1.0 - xc.powi(2)).sqrt();
                let zi = Self::pos_or_neg() * (1.0 - xc.powi(2) - y.powi(2)).max(0.0).sqrt();
                (xc, zi)
            } else {
                let y = self.get_cos_angle_to_x();
                let xc = Self::pos_or_neg() * rand::random::<f64>() * (1.0 - y.powi(2)).sqrt();
                let zi = Self::pos_or_neg() * (1.0 - xc.powi(2) - y.powi(2)).max(0.0).sqrt();
                (xc, zi)
            };
            let y = Self::pos_or_neg() * rand::random::<f64>() * (1.0 - xni.powi(2)).sqrt();
            let xn2 = xni * this_angle.cos() + zni * this_angle.sin();
            let zn2 = -xni * this_angle.sin() + zni * this_angle.cos();
            let th = zn2.acos();
            let ph = if th.sin().abs() > 1e-10 {
                (xn2 / th.sin()).acos()
            } else {
                0.0
            };
            (xn2, y, zn2, th, ph)
        } else {
            let th: f64 = rand::random::<f64>() * 2.0 * PI;
            let ph: f64 = rand::random::<f64>() * 2.0 * PI;
            (th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos(), th, ph)
        };

        self.track_photoelectron(
            coef_calc,
            time_stamp,
            photoelectron_energy,
            xn,
            yn,
            zn,
            x_norm,
            y_norm,
            z_norm,
            theta,
            phi,
            surrounding,
            true,
            beam,
            angle,
            wedge,
        );

        if !surrounding {
            self.ionisations_old += photoelectron_energy * 1000.0 / 21.0;
        }

        // Auger cascade
        if !surrounding {
            self.produce_auger_electron(
                coef_calc,
                time_stamp,
                shell_index,
                &ionised_elem_name,
                xn,
                yn,
                zn,
                surrounding,
                beam,
                angle,
                wedge,
            );
        }
    }

    // ── Auger electron production ─────────────────────────────────────────────

    #[allow(clippy::too_many_arguments)]
    fn produce_auger_electron(
        &mut self,
        coef_calc: &mut dyn CoefCalc,
        time_stamp: f64,
        shell_index: usize,
        element_name: &str,
        xn: f64,
        yn: f64,
        zn: f64,
        surrounding: bool,
        beam: &dyn Beam,
        angle: f64,
        wedge: &Wedge,
    ) {
        let shell_binding_energy = coef_calc.shell_binding_energy(element_name, shell_index);
        let pixel_coord = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
        let z = coef_calc.atomic_number_of(element_name) as u32;
        let shell = shell_index as u32;

        // Check if this element/shell supports Auger
        let can_auger = AUGER_ELEMENTS.contains(&z) && {
            match z {
                6..=8 => shell == 0,
                11 | 12 => shell <= 1,
                14 | 15 | 16 | 17 | 19 | 20 => shell <= 2,
                25..=30 => shell <= 2,
                33 | 34 => shell <= 3,
                _ => false,
            }
        };

        if !can_auger {
            // Deposit energy directly
            let idx = vox3(
                &self.cryst_size_voxels,
                pixel_coord[0],
                pixel_coord[1],
                pixel_coord[2],
            );
            self.gos_electron_dose_v_resolved += shell_binding_energy;
            if idx < self.voxel_energy_v_resolved.len() {
                self.voxel_energy_v_resolved[idx] += shell_binding_energy;
            }
            self.dose += shell_binding_energy;
            if idx < self.dose_simple.len() {
                self.dose_simple[idx] += shell_binding_energy;
            }
            return;
        }

        let fl_yield = coef_calc.shell_fluorescence_yield(element_name, shell_index);
        let rnd: f64 = rand::random();
        let idx = vox3(
            &self.cryst_size_voxels,
            pixel_coord[0],
            pixel_coord[1],
            pixel_coord[2],
        );

        if rnd > fl_yield {
            // Auger emission
            if !self.simple_mc {
                if let Some(em) = self.auger_data.get(&z) {
                    if let Some(td) = em.get(&shell) {
                        let ti = self.get_transition_index_auger(z, shell);
                        if ti < td.energies.len() {
                            let auger_energy = td.energies[ti];
                            let lw = td.linewidths[ti];
                            let exit_idx = td.exit_index[ti] as u32;
                            let drop_idx = td.drop_index[ti] as u32;

                            self.gos_electron_dose_v_resolved +=
                                shell_binding_energy - auger_energy;
                            if idx < self.voxel_energy_v_resolved.len() {
                                self.voxel_energy_v_resolved[idx] +=
                                    shell_binding_energy - auger_energy;
                            }

                            let lifetime =
                                1.0e15 * (H / (2.0 * PI)) / ((lw / 1000.0) * KEV_TO_JOULES);
                            let new_ts = time_stamp + lifetime;

                            let th: f64 = rand::random::<f64>() * 2.0 * PI;
                            let ph: f64 = rand::random::<f64>() * 2.0 * PI;
                            let xn2 = th.sin() * ph.cos();
                            let yn2 = th.sin() * ph.sin();
                            let zn2 = th.cos();

                            if !surrounding {
                                self.add_ionisation_vox(pixel_coord, element_name, xn, yn, beam, 1);
                            }

                            self.track_photoelectron(
                                coef_calc,
                                new_ts,
                                auger_energy,
                                xn,
                                yn,
                                zn,
                                xn2,
                                yn2,
                                zn2,
                                th,
                                ph,
                                surrounding,
                                false,
                                beam,
                                angle,
                                wedge,
                            );

                            if !surrounding {
                                self.ionisations_old += 1.0;
                            }

                            // Cascade
                            let exit_name = element_name.to_string();
                            self.produce_auger_electron(
                                coef_calc,
                                new_ts,
                                exit_idx as usize,
                                &exit_name,
                                xn,
                                yn,
                                zn,
                                surrounding,
                                beam,
                                angle,
                                wedge,
                            );
                            self.produce_auger_electron(
                                coef_calc,
                                new_ts,
                                drop_idx as usize,
                                &exit_name,
                                xn,
                                yn,
                                zn,
                                surrounding,
                                beam,
                                angle,
                                wedge,
                            );
                        }
                    }
                }
            } else {
                // simple MC: deposit
                self.dose += shell_binding_energy;
                if idx < self.dose_simple.len() {
                    self.dose_simple[idx] += shell_binding_energy;
                }
            }
        } else {
            // Fluorescence
            if let Some(em) = self.fl_data.get(&z) {
                if let Some(td) = em.get(&shell) {
                    let ti = self.get_transition_index_fl(z, shell);
                    let fl_energy = if ti < td.energies.len() {
                        td.energies[ti]
                    } else {
                        0.0
                    };
                    if !self.simple_mc {
                        let drop_idx = if ti < td.drop_index.len() {
                            td.drop_index[ti] as u32
                        } else {
                            0
                        };
                        let lw = if ti < td.linewidths.len() {
                            td.linewidths[ti]
                        } else {
                            0.0
                        };
                        self.gos_electron_dose_v_resolved += shell_binding_energy;
                        if idx < self.voxel_energy_v_resolved.len() {
                            self.voxel_energy_v_resolved[idx] += shell_binding_energy;
                        }
                        let lifetime = if lw > 0.0 {
                            1.0e15 * (H / (2.0 * PI)) / ((lw / 1000.0) * KEV_TO_JOULES)
                        } else {
                            0.0
                        };
                        let new_ts = time_stamp + lifetime;
                        let en = element_name.to_string();
                        self.produce_auger_electron(
                            coef_calc,
                            new_ts,
                            drop_idx as usize,
                            &en,
                            xn,
                            yn,
                            zn,
                            surrounding,
                            beam,
                            angle,
                            wedge,
                        );
                    } else {
                        self.dose += shell_binding_energy - fl_energy;
                        if idx < self.dose_simple.len() {
                            self.dose_simple[idx] += shell_binding_energy - fl_energy;
                        }
                    }
                }
            }
        }
    }

    // ── Electron tracking ─────────────────────────────────────────────────────

    #[allow(clippy::too_many_arguments)]
    fn track_photoelectron(
        &mut self,
        coef_calc: &mut dyn CoefCalc,
        starting_ts: f64,
        starting_energy: f64,
        prev_x: f64,
        prev_y: f64,
        prev_z: f64,
        mut x_norm: f64,
        mut y_norm: f64,
        mut z_norm: f64,
        mut theta: f64,
        mut phi: f64,
        mut surrounding: bool,
        primary_electron: bool,
        beam: &dyn Beam,
        angle: f64,
        wedge: &Wedge,
    ) {
        if starting_energy < 0.0 {
            return;
        }

        let energy_loss_to_update = 0.0_f64;
        let mut loss_since_update = 0.0_f64;

        let mut electron_energy = starting_energy;
        let mut time_stamp = starting_ts;
        let mut prev_x = prev_x;
        let mut prev_y = prev_y;
        let mut prev_z = prev_z;
        let mut prev_theta = theta;
        let mut prev_phi = phi;

        // Exit immediately if energy too low
        if electron_energy < 0.05 {
            let pixel = self.convert_to_pixel_coordinates(prev_x, prev_y, prev_z, angle, wedge);
            let idx = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);
            if self.is_microcrystal_at(prev_x, prev_y, prev_z, angle, wedge) {
                if primary_electron {
                    self.dose += electron_energy;
                    self.electron_dose += electron_energy;
                    self.gos_electron_dose += electron_energy;
                    if idx < self.dose_simple.len() {
                        self.dose_simple[idx] += electron_energy;
                    }
                }
                if !self.simple_mc {
                    self.gos_electron_dose_v_resolved += electron_energy;
                    if idx < self.voxel_energy_v_resolved.len() {
                        self.voxel_energy_v_resolved[idx] += electron_energy;
                    }
                }
            }
            return;
        }

        let mut stopping_power = coef_calc.stopping_power(electron_energy, surrounding);
        let mut lambda_el = coef_calc.electron_elastic_mfpl(electron_energy, surrounding);
        let mut elastic_probs = coef_calc.elastic_probs(surrounding);
        let fse_xsection = Self::get_fse_x_section(electron_energy);
        let _fse_lambda = coef_calc.fse_lambda(fse_xsection, surrounding);

        coef_calc.populate_cross_section_coefficients();
        let mut inner_shell_lambda =
            coef_calc.bethe_ionisation_x_section(electron_energy, surrounding);
        let ionisation_probs = coef_calc.all_shell_probs(surrounding);

        let mut gos_inel_lambda = 0.0;
        let mut gos_inner_lambda = 0.0;
        let mut gos_outer_lambda = 0.0;
        let mut gos_ion_probs: HashMap<String, Vec<f64>> = HashMap::new();
        let mut gos_outer_ion_probs: HashMap<String, f64> = HashMap::new();
        let mut p_inner = 0.0_f64;

        if !self.simple_mc {
            gos_inel_lambda = coef_calc.gos_inel(surrounding, electron_energy);
            gos_inner_lambda = coef_calc.gos_inner_lambda(surrounding);
            gos_outer_lambda = coef_calc.gos_outer_lambda(surrounding);
            gos_outer_ion_probs =
                coef_calc.gos_outer_shell_probs_simple(surrounding, gos_outer_lambda);
            gos_ion_probs = coef_calc.gos_shell_probs(surrounding, gos_inel_lambda);
            if inner_shell_lambda > 0.0 {
                gos_inel_lambda = 1.0 / (1.0 / gos_outer_lambda + 1.0 / inner_shell_lambda);
            } else {
                gos_inel_lambda = 1.0 / gos_outer_lambda;
            }
            if inner_shell_lambda > 0.0 {
                p_inner = gos_inel_lambda / inner_shell_lambda;
            }
        }

        let mut lambda_t = if gos_inel_lambda > 0.0 {
            1.0 / (1.0 / lambda_el + 1.0 / gos_inel_lambda)
        } else {
            lambda_el
        };
        let mut p_inel = 1.0 - lambda_t / lambda_el;

        let rnd: f64 = rand::random();
        let mut s = -lambda_t * rnd.ln();
        let mut xn = prev_x + s * x_norm;
        let mut yn = prev_y + s * y_norm;
        let mut zn = prev_z + s * z_norm;

        let mut exited = false;
        let mut entered = false;

        while !exited {
            if self.is_microcrystal_at(xn, yn, zn, angle, wedge) {
                // In crystal: deposit energy
                if surrounding {
                    entered = true;
                    surrounding = false;
                    // Update coefficients for crystal
                    stopping_power = coef_calc.stopping_power(electron_energy, false);
                    lambda_el = coef_calc.electron_elastic_mfpl(electron_energy, false);
                    elastic_probs = coef_calc.elastic_probs(false);
                    if !self.simple_mc {
                        inner_shell_lambda =
                            coef_calc.bethe_ionisation_x_section(electron_energy, false);
                        gos_inel_lambda = coef_calc.gos_inel(false, electron_energy);
                        gos_outer_lambda = coef_calc.gos_outer_lambda(false);
                        gos_outer_ion_probs =
                            coef_calc.gos_outer_shell_probs_simple(false, gos_outer_lambda);
                        if inner_shell_lambda > 0.0 {
                            gos_inel_lambda =
                                1.0 / (1.0 / gos_outer_lambda + 1.0 / inner_shell_lambda);
                            p_inner = gos_inel_lambda / inner_shell_lambda;
                        } else {
                            gos_inel_lambda = 1.0 / gos_outer_lambda;
                            p_inner = 0.0;
                        }
                        gos_ion_probs = coef_calc.gos_shell_probs(false, gos_inel_lambda);
                    }
                    lambda_t = if gos_inel_lambda > 0.0 {
                        1.0 / (1.0 / lambda_el + 1.0 / gos_inel_lambda)
                    } else {
                        lambda_el
                    };
                    p_inel = 1.0 - lambda_t / lambda_el;
                    let r: f64 = rand::random();
                    s = -lambda_t * r.ln();
                    xn = prev_x + s * x_norm;
                    yn = prev_y + s * y_norm;
                    zn = prev_z + s * z_norm;
                    continue;
                }

                let energy_lost = s * stopping_power;
                if self.simple_mc {
                    loss_since_update += energy_lost;
                }
                let time_to_dist = Self::get_time_to_distance(electron_energy, s);
                time_stamp += time_to_dist;
                let energy_to_add = energy_lost;

                if primary_electron {
                    let pixel = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
                    let idx = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);
                    self.dose += energy_to_add;
                    if idx < self.dose_simple.len() {
                        self.dose_simple[idx] += energy_to_add;
                    }
                    if entered {
                        self.electron_dose_surrounding += energy_to_add;
                    } else {
                        self.electron_dose += energy_to_add;
                    }
                }

                prev_theta = theta;
                prev_phi = phi;
                prev_x = xn;
                prev_y = yn;
                prev_z = zn;

                // Elastic or inelastic
                let rnd_inel: f64 = rand::random();
                if rnd_inel < p_inel && !self.simple_mc {
                    // GOS inelastic collision
                    let pixel = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
                    let idx = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);

                    // Determine inner vs outer shell
                    let rnd_inner: f64 = rand::random();
                    let elem_rnd: f64 = rand::random();
                    let mut collided_elem: Option<String> = None;
                    let mut collided_shell: i32 = -1;
                    let mut plasmon = false;

                    if rnd_inner < p_inner {
                        // Inner shell
                        for (name, probs) in &ionisation_probs {
                            let s = Self::find_if_element_ionised(probs, elem_rnd);
                            if s >= 0 {
                                collided_elem = Some(name.clone());
                                collided_shell = s;
                                break;
                            }
                        }
                    } else {
                        // Outer shell
                        for (name, &p) in &gos_outer_ion_probs {
                            if p > elem_rnd {
                                collided_elem = Some(name.clone());
                                collided_shell = 999;
                                break;
                            }
                        }
                    }
                    if collided_shell == -1 {
                        plasmon = true;
                    }

                    let shell_binding = if plasmon || collided_shell == 999 {
                        0.0
                    } else {
                        let cs = collided_shell as usize;
                        coef_calc.shell_binding_energy(collided_elem.as_deref().unwrap_or(""), cs)
                    };

                    let wk = if plasmon {
                        coef_calc.wcb_all(surrounding)
                    } else {
                        coef_calc.wk_molecule(
                            coef_calc.return_adjustment(),
                            collided_elem.as_deref().unwrap_or(""),
                            collided_shell.max(0) as usize,
                            surrounding,
                        )
                    };

                    let uk = if plasmon { 0.0 } else { shell_binding * 1000.0 };
                    let wak = Self::wk_to_wak(electron_energy, wk, uk);
                    let wdis = 3.0 * wak - 2.0 * uk;
                    let qak = Self::get_qak(electron_energy, wk, uk);

                    // Get collision type (long/transverse/close)
                    let gos_var = coef_calc.gos_variable(surrounding);
                    let col_type = if !plasmon {
                        if let Some(elem_var) = gos_var.get(collided_elem.as_deref().unwrap_or(""))
                        {
                            let arr: [f64; 4] = [
                                elem_var
                                    .get(collided_shell.max(0) as usize)
                                    .copied()
                                    .unwrap_or(0.0),
                                elem_var
                                    .get(collided_shell.max(0) as usize + 1)
                                    .copied()
                                    .unwrap_or(0.0),
                                elem_var
                                    .get(collided_shell.max(0) as usize + 2)
                                    .copied()
                                    .unwrap_or(0.0),
                                elem_var.iter().sum(),
                            ];
                            Self::get_gos_inelastic_type(&arr, 0)
                        } else {
                            2
                        }
                    } else {
                        let pv = coef_calc.plasmon_variable(surrounding);
                        let arr = [pv, pv, pv, pv * 3.0];
                        Self::get_gos_inelastic_type(&arr, 0)
                    };

                    let w = if col_type == 0 || col_type == 1 {
                        if plasmon {
                            wk / 1000.0
                        } else {
                            Self::get_energy_loss_distant(wdis, uk) / 1000.0
                        }
                    } else if plasmon {
                        wk / 1000.0
                    } else {
                        let k = Self::sample_k(electron_energy, uk);
                        k * (electron_energy + uk / 1000.0)
                    };

                    let scatter_theta = match col_type {
                        0 => {
                            let q = coef_calc.recoil_energy_distant(electron_energy, wak, qak);
                            Self::get_gos_primary_scatter_long(electron_energy, q, wak)
                        }
                        1 => 0.0,
                        _ => Self::get_gos_primary_scatter_close(electron_energy, w),
                    };

                    let se_energy = w - uk / 1000.0;
                    let min_track = 0.05;
                    if se_energy > 0.0 {
                        if !surrounding {
                            if let Some(ref en) = collided_elem {
                                self.add_ionisation_vox(pixel, en, xn, yn, beam, 1);
                            }
                        }
                        self.avg_uk += uk;
                        self.avg_uk_num += 1;

                        if se_energy > min_track {
                            let se_def_theta = if col_type == 0 || col_type == 1 {
                                let q = coef_calc.recoil_energy_distant(electron_energy, wak, qak);
                                Self::secondary_scatter_distant(electron_energy, wak, q)
                            } else {
                                Self::secondary_scatter_close(electron_energy, w)
                            };
                            let se_def_phi: f64 = rand::random::<f64>() * 2.0 * PI;
                            let new_vec = Self::get_new_direction_vector(
                                x_norm,
                                y_norm,
                                z_norm,
                                se_def_theta,
                                se_def_phi,
                            );
                            self.track_photoelectron(
                                coef_calc,
                                time_stamp,
                                se_energy,
                                xn,
                                yn,
                                zn,
                                new_vec[0],
                                new_vec[1],
                                new_vec[2],
                                new_vec[2].acos(),
                                0.0,
                                surrounding,
                                false,
                                beam,
                                angle,
                                wedge,
                            );
                            if primary_electron {
                                self.gos_electron_dose += w;
                                self.avg_w += w;
                                self.avg_w_num += 1;
                            }
                        } else {
                            if primary_electron {
                                self.gos_electron_dose += w;
                            }
                            self.gos_electron_dose_v_resolved += se_energy;
                            if idx < self.voxel_energy_v_resolved.len() {
                                self.voxel_energy_v_resolved[idx] += se_energy;
                            }
                            let avg_e = coef_calc.plasma_energy(surrounding) / 1000.0;
                            if avg_e > 0.0 && !surrounding {
                                let n_ion = (se_energy / avg_e) as u64;
                                if n_ion > 0 {
                                    self.low_energy_ionisations += n_ion;
                                }
                            }
                        }
                        loss_since_update += w;
                        if let Some(ref en) = collided_elem {
                            self.produce_auger_electron(
                                coef_calc,
                                time_stamp,
                                collided_shell.max(0) as usize,
                                en,
                                xn,
                                yn,
                                zn,
                                surrounding,
                                beam,
                                angle,
                                wedge,
                            );
                        }
                    }

                    let def_phi: f64 = rand::random::<f64>() * 2.0 * PI;
                    let new_vec = Self::get_new_direction_vector(
                        x_norm,
                        y_norm,
                        z_norm,
                        scatter_theta,
                        def_phi,
                    );
                    x_norm = new_vec[0];
                    y_norm = new_vec[1];
                    z_norm = new_vec[2];
                } else {
                    // Elastic
                    let def_theta =
                        self.get_scattering_theta(electron_energy, &elastic_probs.clone());
                    let def_phi = Self::get_scattering_phi();
                    let new_vec =
                        Self::get_new_direction_vector(x_norm, y_norm, z_norm, def_theta, def_phi);
                    x_norm = new_vec[0];
                    y_norm = new_vec[1];
                    z_norm = new_vec[2];
                }

                theta = z_norm.acos();
                phi = if theta.sin().abs() > 1e-10 {
                    (x_norm / theta.sin()).acos()
                } else {
                    0.0
                };

                // Update energy
                if electron_energy > 0.5 {
                    if loss_since_update > energy_loss_to_update {
                        electron_energy -= loss_since_update;
                        if !self.simple_mc {
                            gos_inel_lambda = coef_calc.gos_inel(false, electron_energy);
                            inner_shell_lambda =
                                coef_calc.bethe_ionisation_x_section(electron_energy, false);
                            gos_outer_lambda = coef_calc.gos_outer_lambda(surrounding);
                            if inner_shell_lambda > 0.0 {
                                gos_inel_lambda =
                                    1.0 / (1.0 / gos_outer_lambda + 1.0 / inner_shell_lambda);
                            } else {
                                gos_inel_lambda = 1.0 / gos_outer_lambda;
                            }
                        }
                        loss_since_update = 0.0;
                    }
                } else {
                    electron_energy -= loss_since_update;
                    loss_since_update = 0.0;
                }

                if electron_energy <= 0.0 {
                    exited = true;
                    break;
                }

                stopping_power = coef_calc.stopping_power(electron_energy, false);
                lambda_el = coef_calc.electron_elastic_mfpl(electron_energy, false);
                elastic_probs = coef_calc.elastic_probs(false);
                lambda_t = if gos_inel_lambda > 0.0 {
                    1.0 / (1.0 / lambda_el + 1.0 / gos_inel_lambda)
                } else {
                    lambda_el
                };
                p_inel = 1.0 - lambda_t / lambda_el;
                let r: f64 = rand::random();
                s = -lambda_t * r.ln();
                xn = prev_x + s * x_norm;
                yn = prev_y + s * y_norm;
                zn = prev_z + s * z_norm;
            } else {
                // Left the crystal
                if !surrounding {
                    if !coef_calc.is_cryo() {
                        exited = true;
                        if self.simple_mc {
                            let esc_dist = 1000.0
                                * self.get_intersection_distance(
                                    prev_x, prev_y, prev_z, x_norm, y_norm, z_norm, true, angle,
                                    wedge,
                                );
                            let energy_lost = esc_dist * stopping_power;
                            electron_energy -= energy_lost;
                            self.dose += energy_lost;
                        }
                    } else {
                        surrounding = true;
                        exited = true; // simplified: don't track in surrounding
                    }
                } else {
                    exited = true;
                }
            }

            if electron_energy < 0.05 {
                exited = true;
                let _px = if !exited { xn } else { prev_x };
                let _py = if !exited { yn } else { prev_y };
                let _pz = if !exited { zn } else { prev_z };
                if self.is_microcrystal_at(prev_x, prev_y, prev_z, angle, wedge) {
                    let pixel =
                        self.convert_to_pixel_coordinates(prev_x, prev_y, prev_z, angle, wedge);
                    let idx = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);
                    if primary_electron {
                        self.dose += electron_energy;
                        if idx < self.dose_simple.len() {
                            self.dose_simple[idx] += electron_energy;
                        }
                        self.gos_electron_dose += electron_energy;
                        if entered {
                            self.electron_dose_surrounding += electron_energy;
                        } else {
                            self.electron_dose += electron_energy;
                        }
                    }
                    if !self.simple_mc {
                        self.gos_electron_dose_v_resolved += electron_energy;
                        if idx < self.voxel_energy_v_resolved.len() {
                            self.voxel_energy_v_resolved[idx] += electron_energy;
                        }
                        let avg_e = coef_calc.plasma_energy(surrounding) / 1000.0;
                        if avg_e > 0.0 {
                            let n_ion = (electron_energy / avg_e) as u64;
                            if n_ion > 0 {
                                self.low_energy_ionisations += n_ion;
                            }
                        }
                    }
                }
            }
        }
    }

    // ── Main entry point ──────────────────────────────────────────────────────

    /// Main entry point: sets up simulation parameters and runs Monte Carlo.
    pub fn calculate_xfel(&mut self, beam: &dyn Beam, wedge: &Wedge, coef_calc: &mut dyn CoefCalc) {
        if self.do_xfel {
            self.pulse_energy = beam.pulse_energy() / 1000.0;
        } else {
            let flux = beam.photons_per_sec();
            let n_photons = flux * wedge.total_sec();
            self.pulse_energy = n_photons * (beam.photon_energy() * KEV_TO_JOULES);
        }

        self.pulse_length = wedge.total_sec().round();
        self.pulse_bin_length = if self.pulse_length <= 1.0 {
            0.01
        } else if self.pulse_length <= 50.0 {
            0.1
        } else {
            1.0
        };

        let num = coef_calc.number_simulated_electrons();
        self.num_photons = if num == 0 { 1_000_000 } else { num };

        self.mean_energy = beam.photon_energy();
        self.energy_fwhm = beam.energy_fwhm();

        // Build photon energy array
        self.photon_energy_array = if let Some(fwhm) = self.energy_fwhm {
            sample_normal_energies(self.mean_energy, fwhm, self.num_photons)
        } else {
            vec![self.mean_energy; self.num_photons as usize]
        };

        self.last_time = (1.0 / C) * (self.z_dimension / 1.0e9) * 1.0e15 + self.pulse_length;

        let start = std::time::Instant::now();
        self.start_monte_carlo_xfel(beam, wedge, coef_calc);
        self.process_dose(beam, coef_calc);
        println!(
            "Monte Carlo simulation complete, runtime in seconds was: {:.8e}",
            start.elapsed().as_secs_f64()
        );
    }

    // ── Dose processing ───────────────────────────────────────────────────────

    pub fn process_dose(&mut self, beam: &dyn Beam, coef_calc: &dyn CoefCalc) {
        let _vox_length = self.cryst_pix_per_um.recip() * 1000.0; // nm
        let vox_vol_cm3 = self.cryst_pix_per_um.recip().powi(3) * 1.0e-12;
        let vox_mass_kg = (coef_calc.density() * vox_vol_cm3) / 1000.0;
        let sample_vol_cm3 = self.get_sample_volume() * 1.0e-21;
        let sample_mass_kg = (coef_calc.density() * sample_vol_cm3) / 1000.0;
        let mean_e_joules = self.mean_energy * KEV_TO_JOULES;
        let n_photons = self.pulse_energy / mean_e_joules;
        let scale = n_photons / self.num_photons as f64;

        self.tot_elastic *= scale;
        self.dose = (self.dose * scale * KEV_TO_JOULES / sample_mass_kg) / 1.0e6;
        self.electron_dose = (self.electron_dose * scale * KEV_TO_JOULES / sample_mass_kg) / 1.0e6;
        self.gos_electron_dose =
            (self.gos_electron_dose * scale * KEV_TO_JOULES / sample_mass_kg) / 1.0e6;
        self.photon_dose = (self.photon_dose * scale * KEV_TO_JOULES / sample_mass_kg) / 1.0e6;
        self.raddose_style_dose =
            (self.raddose_style_dose * scale * KEV_TO_JOULES / sample_mass_kg) / 1.0e6;
        self.raddose_style_dose_compton =
            (self.raddose_style_dose_compton * scale * KEV_TO_JOULES / sample_mass_kg) / 1.0e6;

        let size = self.cryst_size_voxels;
        for a in 0..size[0] {
            for b in 0..size[1] {
                for c in 0..size[2] {
                    let idx = vox3(&size, a, b, c);
                    if idx < self.voxel_energy_v_resolved.len() {
                        self.voxel_energy_v_resolved[idx] *= scale * KEV_TO_JOULES;
                    }
                    if idx < self.voxel_ionisations_v_resolved.len() {
                        self.voxel_ionisations_v_resolved[idx] *= scale;
                    }
                    if idx < self.voxel_elastic.len() {
                        self.voxel_elastic[idx] *= scale;
                    }
                }
            }
        }

        // Compute summary statistics
        let exposed_area = self.get_exposed_area(beam);
        let exposed_vol = exposed_area * self.z_dimension * 1.0e-21;
        let exposed_vol = exposed_vol.min(sample_vol_cm3);
        let _vol_fraction = exposed_vol / sample_vol_cm3;
        let exposed_mass = (coef_calc.density() * exposed_vol) / 1000.0;

        let mut vox_dose_v_resolved = 0.0_f64;
        let mut vox_dose_exposed = 0.0_f64;
        let mut dwd = 0.0_f64;
        let total_elastic = self.tot_elastic;

        for a in 0..size[0] {
            for b in 0..size[1] {
                for c in 0..size[2] {
                    let idx = vox3(&size, a, b, c);
                    if idx >= self.voxel_energy_v_resolved.len() {
                        continue;
                    }
                    let cart = self.convert_to_cartesian_coordinates(a, b, c);
                    if c < self.last_time_vox.len() {
                        vox_dose_v_resolved += self.voxel_energy_v_resolved[idx];
                        if self.test_if_inside_exposed_area(cart[0], cart[1], beam) {
                            vox_dose_exposed += self.voxel_energy_v_resolved[idx];
                        }
                        if idx < self.voxel_elastic.len()
                            && self.voxel_elastic[idx] > 0.0
                            && total_elastic > 0.0
                        {
                            dwd += self.voxel_energy_v_resolved[idx]
                                * (self.voxel_elastic[idx] / total_elastic);
                        }
                    }
                }
            }
        }

        let vdvr_mgy = (vox_dose_v_resolved / sample_mass_kg) / 1.0e6;
        let vde_mgy = if exposed_mass > 0.0 {
            (vox_dose_exposed / exposed_mass) / 1.0e6
        } else {
            0.0
        };
        let _dwd_mgy = if vox_mass_kg > 0.0 {
            (dwd / vox_mass_kg) / 1.0e6
        } else {
            0.0
        };

        let _frac_elastic =
            1.0 - (-coef_calc.elastic_coefficient() * (self.z_dimension / 1000.0)).exp();
        let _frac_compton =
            1.0 - (-coef_calc.inelastic_coefficient() * (self.z_dimension / 1000.0)).exp();

        let tot_raddose = self.raddose_style_dose + self.raddose_style_dose_compton;
        let exposed_mass_sample = (coef_calc.density() * exposed_vol) / 1000.0;
        let rd_exposed = if exposed_mass_sample > 0.0 {
            tot_raddose * sample_mass_kg / exposed_mass_sample
        } else {
            0.0
        };

        // Print summary
        println!();
        println!(
            "RADDOSE-3D style average dose whole crystal: {:.3}",
            tot_raddose
        );
        println!(
            "RADDOSE-3D style average dose exposed region: {:.3}",
            rd_exposed
        );
        if self.do_xfel {
            println!(
                "RADDOSE-XFEL average dose whole crystal (ADWC): {:.3}",
                vdvr_mgy
            );
            println!(
                "RADDOSE-XFEL average dose exposed region (ADER): {:.3}",
                vde_mgy
            );
        } else if !self.simple_mc {
            println!(
                "RADDOSE-GOS average dose whole crystal (ADWC): {:.3}",
                vdvr_mgy
            );
            println!(
                "RADDOSE-GOS average dose exposed region (ADER): {:.3}",
                vde_mgy
            );
        } else {
            println!(
                "Monte Carlo average dose whole crystal (ADWC): {:.3}",
                self.dose
            );
        }

        // Write output CSV
        if let Err(e) =
            self.write_output_csv("outputXFEL.CSV", tot_raddose, rd_exposed, vdvr_mgy, vde_mgy)
        {
            println!("Could not write outputXFEL.CSV: {}", e);
        }
    }

    fn write_output_csv(
        &self,
        filename: &str,
        rd_adwc: f64,
        rd_ader: f64,
        xfel_adwc: f64,
        xfel_ader: f64,
    ) -> std::io::Result<()> {
        use std::io::Write;
        let append = self.run_number > 1;
        let mut f = if append {
            std::fs::OpenOptions::new().append(true).open(filename)?
        } else {
            std::fs::File::create(filename)?
        };
        if !append {
            writeln!(f, "Run Number,RD3D-ADWC,RD3D-ADER,XFEL-ADWC,XFEL-ADER")?;
        }
        writeln!(
            f,
            " {},{:.6},{:.6},{:.6},{:.6}",
            self.run_number, rd_adwc, rd_ader, xfel_adwc, xfel_ader
        )?;
        Ok(())
    }
}

// ── Pre-parsed binary data (generated by build.rs) ────────────────────────────

include!(concat!(env!("OUT_DIR"), "/elsepa_arrays.rs"));
include!(concat!(env!("OUT_DIR"), "/transition_data.rs"));

// ── Binary deserializers ───────────────────────────────────────────────────────

/// Deserialize a pre-parsed ELSEPA binary blob.
///
/// Format: [n_rows: u32 LE][n_cols: u32 LE]
///         [n_rows × u64 LE]  energy keys (f64::to_bits)
///         [n_rows × n_cols × f32 LE]  probabilities, row-major
pub(super) fn parse_elsepa_bin(data: &[u8]) -> BTreeMap<u64, Vec<f64>> {
    if data.len() < 8 {
        return BTreeMap::new();
    }
    let n_rows = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    let n_cols = u32::from_le_bytes(data[4..8].try_into().unwrap()) as usize;
    let keys_base = 8;
    let probs_base = keys_base + n_rows * 8;
    if data.len() < probs_base + n_rows * n_cols * 4 {
        return BTreeMap::new();
    }
    let mut map = BTreeMap::new();
    for row in 0..n_rows {
        let k_off = keys_base + row * 8;
        let key = u64::from_le_bytes(data[k_off..k_off + 8].try_into().unwrap());
        let p_base = probs_base + row * n_cols * 4;
        let probs: Vec<f64> = (0..n_cols)
            .map(|col| {
                let off = p_base + col * 4;
                f32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as f64
            })
            .collect();
        map.insert(key, probs);
    }
    map
}

/// Deserialize a pre-parsed transition binary blob.
///
/// Format: [n: u32 LE]
///         [n × f32 LE] × 6  linewidths, probs, energies,
///                            cumulative_probs, exit_index, drop_index
pub(super) fn parse_transition_bin(data: &[u8]) -> Option<TransitionData> {
    if data.len() < 4 {
        return None;
    }
    let n = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    if n == 0 || data.len() < 4 + n * 6 * 4 {
        return None;
    }
    let read_vec = |base: usize| -> Vec<f64> {
        (0..n)
            .map(|i| {
                let off = base + i * 4;
                f32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as f64
            })
            .collect()
    };
    Some(TransitionData {
        linewidths: read_vec(4),
        probs: read_vec(4 + n * 4),
        energies: read_vec(4 + n * 8),
        cumulative_probs: read_vec(4 + n * 12),
        exit_index: read_vec(4 + n * 16),
        drop_index: read_vec(4 + n * 20),
    })
}

// ── Data loading (thin wrappers over binary deserializers) ─────────────────────

fn load_auger_fl_csv(_dir: &str, filename: &str) -> Option<TransitionData> {
    let stem = filename.trim_end_matches(".csv");
    parse_transition_bin(get_transition_bin("auger", stem)?)
}

fn load_fl_csv(_dir: &str, filename: &str) -> Option<TransitionData> {
    let stem = filename.trim_end_matches(".csv");
    parse_transition_bin(get_transition_bin("fl", stem)?)
}

fn load_angle_file(high_energy: bool, atomic_num: usize) -> Option<BTreeMap<u64, Vec<f64>>> {
    let data: &[u8] = if high_energy {
        ELSEPA_ABOVE.get(atomic_num)?
    } else {
        ELSEPA_BELOW.get(atomic_num)?
    };
    if data.is_empty() {
        return None;
    }
    Some(parse_elsepa_bin(data))
}

/// Sample photon energies from a normal distribution.
fn sample_normal_energies(mean: f64, fwhm: f64, n: u64) -> Vec<f64> {
    crate::energy_distribution::sample_normal_energies(mean, fwhm, n as usize)
}
