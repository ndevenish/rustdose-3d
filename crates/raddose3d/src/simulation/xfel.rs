// Many variables in this module are intermediate simulation-state trackers
// that will be connected to outputs once the full XFEL pipeline is wired up.
#![allow(dead_code, unused_variables, unused_assignments)]
/// XFEL (X-ray Free Electron Laser) simulation engine.
///
/// Ports Java XFEL.java (~4387 lines).
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;

use crate::beam::Beam;
use crate::coefcalc::CoefCalc;
use crate::wedge::Wedge;

// Re-use binary data loaders from mc module.
use super::mc::{get_transition_bin, parse_elsepa_bin, parse_transition_bin};

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

fn get_num_auger_shells(idx: usize) -> u32 {
    match AUGER_ELEMENTS[idx] {
        6..=8 => 0,
        11 => 1,
        12 => 2,
        14 | 15 | 16 | 17 | 19 | 20 => 3,
        25..=30 => 6,
        33 | 34 => 8,
        _ => 0,
    }
}

/// Auger/FL transition data for one element+shell.
#[derive(Default, Clone)]
pub struct XfelTransitionData {
    pub linewidths: Vec<f64>,
    pub probs: Vec<f64>,
    pub energies: Vec<f64>,
    pub cumulative_probs: Vec<f64>,
    pub exit_index: Vec<f64>,
    pub drop_index: Vec<f64>,
}

type XfelShellMap = HashMap<u32, XfelTransitionData>;
type XfelTransitionMap = HashMap<u32, XfelShellMap>;

// ── Flat 3-D index ─────────────────────────────────────────────────────────
fn vox3(size: &[usize; 3], i: usize, j: usize, k: usize) -> usize {
    i * size[1] * size[2] + j * size[2] + k
}

fn f64_key(x: f64) -> u64 {
    x.to_bits()
}

/// XFEL simulation state — time-resolved variant of MC.
pub struct XfelSimulation {
    // Crystal geometry
    pub vertices: Vec<[f64; 3]>,
    pub indices: Vec<[usize; 3]>,
    pub cryst_coord: Vec<[f64; 3]>,
    pub cryst_pix_per_um: f64,
    pub cryst_size_voxels: [usize; 3],
    pub cryst_occ: Vec<bool>,
    pub run_number: usize,

    // Cached normals
    normals: Option<Vec<[f64; 3]>>,
    origin_distances: Option<Vec<f64>>,
    rotated_normals: Option<Vec<[f64; 3]>>,
    rotated_origin_distances: Option<Vec<f64>>,

    // Dimensions in nm
    pub x_dimension: f64,
    pub y_dimension: f64,
    pub z_dimension: f64,

    // Time-resolved dose arrays (indexed by time bin)
    pub dose: Vec<f64>,
    pub photon_dose: Vec<f64>,
    pub electron_dose: Vec<f64>,
    pub gos_electron_dose: Vec<f64>,
    pub photon_dose_vr: Vec<f64>,
    pub gos_electron_dose_vr: Vec<f64>,
    pub electron_dose_surrounding: Vec<f64>,

    // Scalar metrics
    pub raddose_style_dose: f64,
    pub raddose_style_dose_compton: f64,
    pub escaped_energy: f64,

    // Ionisation tracking (time-resolved)
    pub total_ionisation_events: Vec<i64>,
    pub total_ionisation_events_vr: Vec<i64>,
    pub low_energy_ionisations: Vec<i64>,
    pub total_ionisation_events_per_atom: Vec<f64>,
    pub low_energy_ionisations_per_atom: Vec<f64>,
    pub total_ionisation_events_per_non_h_atom: Vec<f64>,
    pub total_ionisation_events_per_atom_vr: Vec<f64>,
    pub total_ionisation_events_per_non_h_atom_vr: Vec<f64>,

    // Per-element ionisation tracking
    pub atomic_ionisations: HashMap<String, i64>,
    pub atomic_ionisations_per_atom: HashMap<String, f64>,
    pub atomic_ionisations_exposed: HashMap<String, i64>,
    pub atomic_ionisations_per_atom_exposed: HashMap<String, f64>,

    // 4D voxel arrays: outer = flat voxel index, inner = time bin
    pub voxel_energy_vr: Vec<Vec<f64>>,
    pub voxel_ionisations_vr: Vec<Vec<f64>>,
    pub dose_simple: Vec<Vec<f64>>,
    pub voxel_elastic: Vec<f64>,
    tot_elastic: f64,

    // Timing
    pub pulse_length: f64,     // fs
    pub pulse_bin_length: f64, // fs
    pub pulse_energy: f64,     // J
    pub last_time: f64,
    pub last_time_vox: Vec<f64>,
    pub num_photons: u64,

    // Angular emission distribution
    angular_emission_probs: Vec<f64>,
    number_angular_emission_bins: usize,

    // Elastic angle lookup tables
    low_energy_angles: Vec<Option<BTreeMap<u64, Vec<f64>>>>,
    high_energy_angles: Vec<Option<BTreeMap<u64, Vec<f64>>>>,

    // Auger transition data
    pub auger_data: XfelTransitionMap,
    pub tot_k_auger_prob: HashMap<u32, f64>,
    // Fluorescence transition data
    pub fl_data: XfelTransitionMap,
    pub fl_cumulative: HashMap<u32, HashMap<u32, Vec<f64>>>,

    // Misc ionisation counters
    ionisations_old: f64,
    ionisations_per_photoelectron: i32,
    total_shell_binding_energy: f64,
    avg_w: f64,
    avg_w_num: u32,
    avg_uk: f64,
    avg_uk_num: u32,

    // Configuration
    pub simple_mc: bool,
    pub do_xfel: bool,
    pub vertical_goni: bool,
    pub vertical_pol: bool,

    // Energy / beam
    pub mean_energy: f64,
    pub energy_fwhm: Option<f64>,
    pub photon_energy_array: Vec<f64>,

    // Time bin count (pre-computed)
    num_time_bins: usize,
}

impl XfelSimulation {
    /// Construct a new XFEL simulation.
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
        wedge_tot_sec: f64,
        vertical_polarisation: bool,
    ) -> Self {
        let x_mm = Self::min_max_static(0, &vertices);
        let y_mm = Self::min_max_static(1, &vertices);
        let z_mm = Self::min_max_static(2, &vertices);
        let x_dimension = 1000.0 * (x_mm[1] - x_mm[0]);
        let y_dimension = 1000.0 * (y_mm[1] - y_mm[0]);
        let z_dimension = 1000.0 * (z_mm[1] - z_mm[0]);

        // Determine pulse length / bin length from wedge total seconds
        let pulse_length = wedge_tot_sec.round().max(1.0);
        let pulse_bin_length = if pulse_length <= 1.0 {
            0.01
        } else if pulse_length <= 50.0 {
            0.1
        } else {
            1.0
        };

        let last_time = (1.0 / C) * (z_dimension / 1.0e9) * 1.0e15 + pulse_length;
        let num_time_bins = (last_time / pulse_bin_length).ceil() as usize + 2;

        let nz = cryst_size_voxels[2];
        let n = cryst_size_voxels[0] * cryst_size_voxels[1] * cryst_size_voxels[2];

        // 4D voxel arrays: each voxel gets a vector of time bins
        let voxel_energy_vr = vec![vec![0.0_f64; num_time_bins]; n];
        let voxel_ionisations_vr = vec![vec![0.0_f64; num_time_bins]; n];
        let dose_simple = vec![vec![0.0_f64; num_time_bins]; n];

        XfelSimulation {
            vertices,
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
            x_dimension,
            y_dimension,
            z_dimension,
            dose: vec![0.0; num_time_bins],
            photon_dose: vec![0.0; num_time_bins],
            electron_dose: vec![0.0; num_time_bins],
            gos_electron_dose: vec![0.0; num_time_bins],
            photon_dose_vr: vec![0.0; num_time_bins],
            gos_electron_dose_vr: vec![0.0; num_time_bins],
            electron_dose_surrounding: vec![0.0; num_time_bins],
            raddose_style_dose: 0.0,
            raddose_style_dose_compton: 0.0,
            escaped_energy: 0.0,
            total_ionisation_events: vec![0; num_time_bins],
            total_ionisation_events_vr: vec![0; num_time_bins],
            low_energy_ionisations: vec![0; num_time_bins],
            total_ionisation_events_per_atom: vec![0.0; num_time_bins],
            low_energy_ionisations_per_atom: vec![0.0; num_time_bins],
            total_ionisation_events_per_non_h_atom: vec![0.0; num_time_bins],
            total_ionisation_events_per_atom_vr: vec![0.0; num_time_bins],
            total_ionisation_events_per_non_h_atom_vr: vec![0.0; num_time_bins],
            atomic_ionisations: HashMap::new(),
            atomic_ionisations_per_atom: HashMap::new(),
            atomic_ionisations_exposed: HashMap::new(),
            atomic_ionisations_per_atom_exposed: HashMap::new(),
            voxel_energy_vr,
            voxel_ionisations_vr,
            dose_simple,
            voxel_elastic: vec![0.0; n],
            tot_elastic: 0.0,
            pulse_length,
            pulse_bin_length,
            pulse_energy: 1.4e-3,
            last_time,
            last_time_vox: vec![0.0; nz],
            num_photons: 1_000_000,
            angular_emission_probs: vec![0.0; 50],
            number_angular_emission_bins: 50,
            low_energy_angles: vec![None; 95],
            high_energy_angles: vec![None; 95],
            auger_data: HashMap::new(),
            tot_k_auger_prob: HashMap::new(),
            fl_data: HashMap::new(),
            fl_cumulative: HashMap::new(),
            ionisations_old: 0.0,
            ionisations_per_photoelectron: 0,
            total_shell_binding_energy: 0.0,
            avg_w: 0.0,
            avg_w_num: 0,
            avg_uk: 0.0,
            avg_uk_num: 0,
            simple_mc: !gos,
            do_xfel: xfel,
            vertical_goni: vertical_goniometer,
            vertical_pol: vertical_polarisation,
            mean_energy: 12.4,
            energy_fwhm: None,
            photon_energy_array: Vec::new(),
            num_time_bins,
        }
    }

    // ── Static geometry helpers ──────────────────────────────────────────────

    fn min_max_static(dim: usize, verts: &[[f64; 3]]) -> [f64; 2] {
        let mut mn = f64::INFINITY;
        let mut mx = f64::NEG_INFINITY;
        for v in verts {
            if v[dim] < mn {
                mn = v[dim];
            }
            if v[dim] > mx {
                mx = v[dim];
            }
        }
        [mn, mx]
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

    fn ensure_normals(&mut self) {
        if self.normals.is_none() {
            let (n, d) = Self::calculate_normals_for(&self.vertices, &self.indices);
            self.normals = Some(n);
            self.origin_distances = Some(d);
        }
    }

    // ── Point-in-polyhedron (ray-casting along +Z) ────────────────────────────

    fn point_in_triangle_2d(p: &[f64; 3], a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> bool {
        let (px, py) = (p[0], p[1]);
        let (ax, ay) = (a[0], a[1]);
        let (bx, by) = (b[0], b[1]);
        let (cx, cy) = (c[0], c[1]);
        let d1 = (px - bx) * (ay - by) - (ax - bx) * (py - by);
        let d2 = (px - cx) * (by - cy) - (bx - cx) * (py - cy);
        let d3 = (px - ax) * (cy - ay) - (cx - ax) * (py - ay);
        let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
        let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
        !(has_neg && has_pos)
    }

    fn point_in_polyhedron(
        x: f64,
        y: f64,
        z: f64,
        normals: &[[f64; 3]],
        dists: &[f64],
        indices: &[[usize; 3]],
        vertices: &[[f64; 3]],
    ) -> bool {
        let ray = [0.0_f64, 0.0, 1.0];
        let mut count = 0i32;
        for (i, tri) in indices.iter().enumerate() {
            let n = normals[i];
            let d = dists[i];
            let denom = n[0] * ray[0] + n[1] * ray[1] + n[2] * ray[2];
            if denom.abs() < 1e-12 {
                continue;
            }
            let t = -(n[0] * x + n[1] * y + n[2] * z + d) / denom;
            if t < 0.0 {
                continue;
            }
            let hit = [x + t * ray[0], y + t * ray[1], z + t * ray[2]];
            if Self::point_in_triangle_2d(
                &hit,
                &vertices[tri[0]],
                &vertices[tri[1]],
                &vertices[tri[2]],
            ) {
                count += 1;
            }
        }
        count % 2 == 1
    }

    fn calculate_crystal_occupancy(&mut self, x: f64, y: f64, z: f64) -> bool {
        self.ensure_normals();
        let normals = self.normals.as_ref().unwrap().clone();
        let dists = self.origin_distances.as_ref().unwrap().clone();
        let vertices = self.vertices.clone();
        let indices = self.indices.clone();
        Self::point_in_polyhedron(x, y, z, &normals, &dists, &indices, &vertices)
    }

    // ── Coordinate helpers ────────────────────────────────────────────────────

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
        let x = x_nm / 1000.0;
        let y = y_nm / 1000.0;
        let z = z_nm / 1000.0;
        let ux = x * cos_a - z * sin_a;
        let uy = y;
        let uz = x * sin_a + z * cos_a;
        self.calculate_crystal_occupancy(ux, uy, uz)
    }

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

    fn convert_to_cartesian_coordinates(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let [nx, ny, nz] = self.cryst_size_voxels;
        let pix = self.cryst_pix_per_um;
        let x = (i as f64 - nx as f64 / 2.0) / pix * 1000.0;
        let y = (j as f64 - ny as f64 / 2.0) / pix * 1000.0;
        let z = (k as f64 - nz as f64 / 2.0) / pix * 1000.0;
        [x, y, z]
    }

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
        wedge: &Wedge,
    ) -> f64 {
        let step = 1.0_f64;
        let mut dist = 0.0;
        let max_dist = 100_000.0;
        loop {
            dist += step;
            let x = px + dist * dx;
            let y = py + dist * dy;
            let z = pz + dist * dz;
            let inside = self.is_microcrystal_at(x, y, z, angle, wedge);
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
        dist / 1000.0
    }

    // ── Beam area helpers ─────────────────────────────────────────────────────

    fn get_exposed_area(&self, beam: &dyn Beam) -> f64 {
        let bx = beam.beam_size_x();
        let by = beam.beam_size_y();
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

    fn get_photon_beam_xy_pos(&self, beam: &dyn Beam, _wedge: &Wedge) -> [f64; 2] {
        let bx = beam.beam_size_x();
        let by = beam.beam_size_y();
        let rx: f64 = rand::random();
        let ry: f64 = rand::random();
        let x_nm = 1000.0 * bx * (rx - 0.5);
        let y_nm = 1000.0 * ((ry * by) - by / 2.0);
        [x_nm, y_nm]
    }

    // ── Angular emission ──────────────────────────────────────────────────────

    fn populate_angular_emission_probs(&mut self) {
        let n = self.number_angular_emission_bins;
        self.angular_emission_probs = vec![0.0; n];
        let mut total_area = 0.0;
        let mut last_h = 0.0;
        for i in 0..=100 {
            let a = (PI / 100.0) * i as f64;
            let h = Self::solve_polarisation(a, 1.0, 2.0);
            if i > 0 {
                total_area += (last_h + h) / 2.0 * (PI / 100.0);
            }
            last_h = h;
        }
        let mut cum = 0.0;
        last_h = 0.0;
        for i in 0..=n {
            let a = (PI / n as f64) * i as f64;
            let h = Self::solve_polarisation(a, 1.0, 2.0);
            if i > 0 {
                let area = (last_h + h) / 2.0 * (PI / n as f64);
                cum += area / total_area;
                self.angular_emission_probs[i - 1] = cum;
            }
            last_h = h;
        }
    }

    fn solve_polarisation(phi: f64, pe: f64, beta: f64) -> f64 {
        (pe / (4.0 * PI)) * (1.0 + beta * 0.5 * (3.0 * phi.cos().powi(2) - 1.0))
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
        if rand::random::<f64>() < 0.5 {
            1.0
        } else {
            -1.0
        }
    }

    // ── Auger/FL data loading ─────────────────────────────────────────────────

    pub fn populate_auger_linewidths(&mut self) {
        for (i, &z) in AUGER_ELEMENTS.iter().enumerate() {
            let num_shells = get_num_auger_shells(i);
            let mut element_map: XfelShellMap = HashMap::new();
            for shell in 0..=num_shells {
                let filename = format!("{}-{}.csv", z, shell);
                if let Some(td) = load_xfel_auger_csv("auger_linewidths", &filename) {
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
            let mut element_map: XfelShellMap = HashMap::new();
            let mut cumulative_map: HashMap<u32, Vec<f64>> = HashMap::new();
            for shell in 0..=num_shells {
                let filename = format!("{}-{}.csv", z, shell);
                if let Some(td) = load_xfel_fl_csv("fl_linewidths", &filename) {
                    cumulative_map.insert(shell, td.cumulative_probs.clone());
                    element_map.insert(shell, td);
                }
            }
            self.fl_data.insert(z, element_map);
            self.fl_cumulative.insert(z, cumulative_map);
        }
    }

    // ── Elastic angle lookup (ELSEPA) ─────────────────────────────────────────

    fn get_primary_elastic_scattering_angle(
        &mut self,
        electron_energy: f64,
        atomic_number: usize,
    ) -> f64 {
        if atomic_number >= 95 {
            return 0.0;
        }
        let high_energy = electron_energy > 20.0;
        if high_energy {
            if self.high_energy_angles[atomic_number].is_none() {
                self.high_energy_angles[atomic_number] = load_xfel_angle_file(true, atomic_number);
            }
        } else if self.low_energy_angles[atomic_number].is_none() {
            self.low_energy_angles[atomic_number] = load_xfel_angle_file(false, atomic_number);
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
            let key = Self::nearest_energy_key(map, electron_energy);
            if let Some(probs) = map.get(&key) {
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
        let deg = if high_energy {
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
        deg * PI / 180.0
    }

    // ── Physics helpers ───────────────────────────────────────────────────────

    fn get_fse_x_section(electron_energy: f64) -> f64 {
        let e_charge = 4.803_204_25e-10_f64;
        let m_cgs = 9.109_383_56e-28_f64;
        let c_cgs = 2.997_924_58e10_f64;
        let c2 = (c_cgs / 100.0).powi(2);
        let vo = electron_energy * KEV_TO_JOULES;
        let beta2 = 1.0 - (m_cgs / 1000.0 * c2 / (vo + m_cgs / 1000.0 * c2)).powi(2);
        let v2 = beta2 * c2 * 10_000.0;
        let constant = (2.0 * PI * e_charge.powi(4)) / (m_cgs * v2 * (vo * 1000.0 * 10_000.0));
        let ec = (14.0 / 1000.0) / electron_energy;
        let tau = electron_energy / 511.0;
        let cs = (((2.0 * tau + 1.0) / (tau + 1.0).powi(2)) * ((1.0_f64 / 0.5 - 1.0).ln())
            + (tau / (tau + 1.0)).powi(2)
            - 1.0 / 0.5
            - 1.0 / (0.5 - 1.0))
            - (((2.0 * tau + 1.0) / (tau + 1.0).powi(2)) * ((1.0 / ec - 1.0).ln())
                + (tau / (tau + 1.0)).powi(2)
                - 1.0 / ec
                - 1.0 / (ec - 1.0));
        cs * constant
    }

    fn get_time_to_distance(electron_energy: f64, s_nm: f64) -> f64 {
        let c2 = C * C;
        let vo = electron_energy * KEV_TO_JOULES;
        let beta2 = 1.0 - (M_E * c2 / (vo + M_E * c2)).powi(2);
        let v = (beta2 * c2).sqrt() * 1.0e9 / 1.0e15; // nm/fs
        s_nm / v
    }

    // ── GOS model helpers ─────────────────────────────────────────────────────

    fn wk_to_wak(e: f64, wk: f64, uk: f64) -> f64 {
        if e * 1000.0 > 3.0 * wk - 2.0 * uk {
            wk
        } else {
            (e * 1000.0 + 2.0 * uk) / 3.0
        }
    }

    fn get_qak(e: f64, wk: f64, uk: f64) -> f64 {
        if e * 1000.0 > 3.0 * wk - 2.0 * uk {
            uk
        } else {
            uk * (e * 1000.0 / (3.0 * wk - 2.0 * uk))
        }
    }

    fn get_energy_loss_distant(wdis: f64, uk: f64) -> f64 {
        let rnd: f64 = rand::random();
        wdis - (rnd * (wdis - uk).powi(2)).sqrt()
    }

    fn get_gos_primary_scatter_long(e_kev: f64, q: f64, _wak: f64) -> f64 {
        let mc2 = M_E * C * C;
        let e_j = e_kev * KEV_TO_JOULES;
        if q == 0.0 {
            return 0.0;
        }
        (q / (e_j + mc2)).asin().abs()
    }

    fn get_gos_primary_scatter_close(e_kev: f64, w_kev: f64) -> f64 {
        let mc2 = M_E * C * C;
        let e_j = e_kev * KEV_TO_JOULES;
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

    fn sample_k(_e: f64, qk: f64) -> f64 {
        let rnd: f64 = rand::random();
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
        let z = self.get_elastic_element_z(elastic_probs) as usize;
        self.get_primary_elastic_scattering_angle(electron_energy, z)
    }

    fn get_scattering_phi() -> f64 {
        rand::random::<f64>() * 2.0 * PI
    }

    fn get_elastic_element_z(&self, probs: &HashMap<String, f64>) -> u32 {
        let rnd: f64 = rand::random();
        for (name, &p) in probs {
            if p > rnd {
                if let Ok(z) = name.parse::<u32>() {
                    return z;
                }
                return 6;
            }
        }
        6
    }

    // ── Transition index helpers ──────────────────────────────────────────────

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

    // ── Ionisation helpers ────────────────────────────────────────────────────

    fn get_ionised_element(probs: &HashMap<String, f64>) -> Option<String> {
        let rnd: f64 = rand::random();
        for (name, &p) in probs {
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

    fn gos_inelastic_type(shell_probs: &[f64; 4]) -> usize {
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

    /// Clamp dose_time index to valid range.
    fn safe_dose_time(&self, t: f64) -> usize {
        let idx = (t / self.pulse_bin_length).max(0.0) as usize;
        idx.min(self.num_time_bins.saturating_sub(1))
    }

    /// Add ionisation event (time-resolved).
    #[allow(clippy::too_many_arguments)]
    fn add_ionisation(
        &mut self,
        time_stamp: f64,
        pixel_coord: [usize; 3],
        element_name: &str,
        ionisation_time: usize,
        x_nm: f64,
        y_nm: f64,
        beam: &dyn Beam,
        count: i64,
    ) {
        let it = ionisation_time.min(self.num_time_bins.saturating_sub(1));
        self.total_ionisation_events[it] += count;
        let ltv = if pixel_coord[2] < self.last_time_vox.len() {
            self.last_time_vox[pixel_coord[2]]
        } else {
            self.last_time
        };
        if time_stamp < ltv - self.pulse_bin_length {
            self.total_ionisation_events_vr[it] += count;
            let idx = vox3(
                &self.cryst_size_voxels,
                pixel_coord[0],
                pixel_coord[1],
                pixel_coord[2],
            );
            if idx < self.voxel_ionisations_vr.len() && it < self.voxel_ionisations_vr[idx].len() {
                self.voxel_ionisations_vr[idx][it] += count as f64;
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
    }

    // ── Main entry point ──────────────────────────────────────────────────────

    /// Main XFEL calculation entry point, mirroring Java CalculateXFEL.
    pub fn calculate_xfel(&mut self, beam: &dyn Beam, wedge: &Wedge, coef_calc: &mut dyn CoefCalc) {
        if self.do_xfel {
            self.pulse_energy = beam.pulse_energy() / 1000.0;
        } else {
            let flux = beam.photons_per_sec();
            let n = flux * wedge.total_sec();
            self.pulse_energy = n * (beam.photon_energy() * KEV_TO_JOULES);
        }

        let num = coef_calc.number_simulated_electrons();
        self.num_photons = if num == 0 { 1_000_000 } else { num };

        self.mean_energy = beam.photon_energy();
        self.energy_fwhm = beam.energy_fwhm();

        self.photon_energy_array = if let Some(fwhm) = self.energy_fwhm {
            sample_normal_energies(self.mean_energy, fwhm, self.num_photons)
        } else {
            vec![self.mean_energy; self.num_photons as usize]
        };

        // Compute lastTimeVox
        let nz = self.last_time_vox.len();
        for i in 0..nz {
            self.last_time_vox[i] =
                (1.0 / C) * ((self.z_dimension / 1.0e9) / nz as f64) * (i + 1) as f64 * 1.0e15
                    + self.pulse_length;
        }

        let start = std::time::Instant::now();
        self.start_monte_carlo_xfel(beam, wedge, coef_calc);
        self.process_dose(beam, coef_calc);
        println!(
            "RADDOSE-XFEL simulation complete, runtime in seconds was: {:.8e}",
            start.elapsed().as_secs_f64()
        );
    }

    // ── Photon loop ───────────────────────────────────────────────────────────

    /// Monte Carlo photon loop for XFEL, mirroring Java startMonteCarloXFEL.
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

            // Determine angle
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

            coef_calc.update_coefficients(energy_of_photon);
            let abs_coef = coef_calc.absorption_coefficient();
            let compton_coef = coef_calc.inelastic_coefficient();
            let elastic_coef = coef_calc.elastic_coefficient();

            let photon_mfpl = (1.0 / (abs_coef + compton_coef)) * 1000.0;
            let prob_compton = 1.0 - photon_mfpl / ((1.0 / abs_coef) * 1000.0);
            let total_mfpl = (1.0 / (abs_coef + compton_coef + elastic_coef)) * 1000.0;
            let elastic_prob = elastic_coef / (abs_coef + compton_coef);

            let element_abs_probs = coef_calc.photo_electric_probs_element(energy_of_photon);
            let element_compton_probs = coef_calc.compton_probs_element(energy_of_photon);
            let ionisation_probs = coef_calc.relative_shell_probs(energy_of_photon, false);

            // Time stamp for this photon (which bin it belongs to)
            let time_stamp = (i as f64 / photon_divisions) * self.pulse_bin_length;

            let xy_pos = self.get_photon_beam_xy_pos(beam, wedge);
            let previous_x = xy_pos[0];
            let previous_y = xy_pos[1];
            let previous_z = -self.z_dimension / 2.0;

            // Surrounding check (simplified — no cryo for XFEL stub)
            let surrounding = !self.is_microcrystal_at(previous_x, previous_y, 0.0, angle, wedge);
            let exited = false;

            let rnd_s: f64 = rand::random();
            let s = if surrounding {
                -photon_mfpl * rnd_s.ln()
            } else {
                -total_mfpl * rnd_s.ln()
            };

            let xn = previous_x + s * 0.0;
            let yn = previous_y + s * 0.0;
            let zn = previous_z + s * 1.0; // direction along +Z

            if !exited && self.is_microcrystal_at(xn, yn, zn, angle, wedge) {
                let time_to_point = (1.0 / C) * (s / 1.0e9);
                let t_stamp = time_stamp + time_to_point * 1.0e15;
                let dose_time = self.safe_dose_time(t_stamp);

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
                        let ec_probs = element_compton_probs.clone();
                        self.produce_compton(
                            beam,
                            coef_calc,
                            t_stamp,
                            dose_time,
                            xn,
                            yn,
                            zn,
                            surrounding,
                            energy_of_photon,
                            &ec_probs,
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
                            surrounding,
                            energy_of_photon,
                            angle,
                            wedge,
                        );
                    }
                }
            }
        }
        println!();

        // Update last_time after loop
        self.last_time = (1.0 / C) * (self.z_dimension / 1.0e9) * 1.0e15 + self.pulse_length;
        let nz = self.last_time_vox.len();
        for k in 0..nz {
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
        dose_time: usize,
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
        let mc2 = M_E * C * C;
        let incident_energy = photon_energy * KEV_TO_JOULES;
        let e_comp = ((incident_energy.powi(2) * (1.0 - photon_theta.cos()))
            / (mc2 * (1.0 + (incident_energy / mc2) * (1.0 - photon_theta.cos()))))
            / KEV_TO_JOULES;

        if !surrounding {
            self.raddose_style_dose_compton += e_comp;
        }

        let electron_phi =
            (1.0 / (photon_theta / 2.0).tan() / (1.0 + incident_energy / mc2)).atan();
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
                let ion_time = self.safe_dose_time(time_stamp);
                self.add_ionisation(
                    time_stamp,
                    pixel_coord,
                    &elem_name,
                    ion_time,
                    xn,
                    yn,
                    beam,
                    1,
                );
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

        if !surrounding {
            self.ionisations_old += (e_comp * 1000.0 / 21.0).floor();
        }
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
        dose_time: usize,
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

        let ion_time = self.safe_dose_time(time_stamp);
        if !surrounding {
            self.add_ionisation(
                time_stamp,
                pixel_coord,
                &ionised_elem_name,
                ion_time,
                xn,
                yn,
                beam,
                1,
            );
        }

        let shell_probs = match ionisation_probs.get(&ionised_elem_name) {
            Some(p) => p.clone(),
            None => return,
        };
        let shell_index = Self::get_ionised_shell(&shell_probs);
        let shell_binding_energy = coef_calc.shell_binding_energy(&ionised_elem_name, shell_index);
        let photoelectron_energy = (photon_energy - shell_binding_energy).max(0.0);

        let dt = dose_time.min(self.num_time_bins.saturating_sub(1));

        if !surrounding {
            let idx = vox3(
                &self.cryst_size_voxels,
                pixel_coord[0],
                pixel_coord[1],
                pixel_coord[2],
            );
            self.dose[dt] += shell_binding_energy;
            self.photon_dose[dt] += shell_binding_energy;
            if idx < self.dose_simple.len() && dt < self.dose_simple[idx].len() {
                self.dose_simple[idx][dt] += shell_binding_energy;
            }
            self.raddose_style_dose += photon_energy;
        }

        // Sample photoelectron direction
        let this_angle = 2.0 * PI - angle;
        let polarised: f64 = rand::random();
        let (x_norm, y_norm, z_norm, theta, phi) = if shell_index == 0 && polarised > 0.0 {
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
            self.ionisations_old += (photoelectron_energy * 1000.0 / 21.0).floor();
        }

        if !self.simple_mc && !surrounding {
            self.produce_auger_electron(
                coef_calc,
                time_stamp,
                shell_index,
                &ionised_elem_name.clone(),
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
        if time_stamp >= self.last_time {
            return;
        }
        let shell_binding_energy = coef_calc.shell_binding_energy(element_name, shell_index);
        let pixel_coord = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
        let idx = vox3(
            &self.cryst_size_voxels,
            pixel_coord[0],
            pixel_coord[1],
            pixel_coord[2],
        );
        let z = coef_calc.atomic_number_of(element_name) as u32;
        let shell = shell_index as u32;
        let dose_time = self.safe_dose_time(time_stamp);
        let dt = dose_time.min(self.num_time_bins.saturating_sub(1));

        // Check element/shell supports Auger
        let z_ok = matches!(
            z,
            6 | 7 | 8 | 11 | 12 | 14 | 16 | 17 | 19 | 20 | 25 | 26 | 27 | 28 | 29 | 30 | 33 | 34
        );
        let shell_ok = match z {
            6..=8 => shell == 0,
            11 => shell <= 1,
            12 => shell <= 2,
            14..=20 => shell <= 3,
            25..=30 => shell <= 6,
            33 | 34 => shell <= 8,
            _ => false,
        };

        if !z_ok || !shell_ok {
            // Deposit energy directly
            self.gos_electron_dose_vr_add(dt, shell_binding_energy);
            if idx < self.voxel_energy_vr.len() && dt < self.voxel_energy_vr[idx].len() {
                self.voxel_energy_vr[idx][dt] += shell_binding_energy;
            }
            return;
        }

        let fl_yield = coef_calc.shell_fluorescence_yield(element_name, shell_index);
        let rnd: f64 = rand::random();

        if rnd > fl_yield {
            // Auger emission
            if !self.simple_mc {
                if let Some(em) = self.auger_data.get(&z) {
                    if let Some(td) = em.get(&shell) {
                        let ti = {
                            let rnd2: f64 = rand::random();
                            let mut res = td.cumulative_probs.len().saturating_sub(1);
                            for (k, &p) in td.cumulative_probs.iter().enumerate() {
                                if rnd2 < p {
                                    res = k;
                                    break;
                                }
                            }
                            res
                        };
                        if ti < td.energies.len() {
                            let auger_energy = td.energies[ti];
                            let lw = td.linewidths[ti];
                            let exit_idx = td.exit_index[ti] as u32;
                            let drop_idx = td.drop_index[ti] as u32;

                            self.gos_electron_dose_vr_add(dt, shell_binding_energy - auger_energy);
                            if idx < self.voxel_energy_vr.len()
                                && dt < self.voxel_energy_vr[idx].len()
                            {
                                self.voxel_energy_vr[idx][dt] +=
                                    shell_binding_energy - auger_energy;
                            }

                            let lifetime =
                                1.0e15 * (H / (2.0 * PI)) / ((lw / 1000.0) * KEV_TO_JOULES);
                            let new_ts = time_stamp + lifetime;
                            let new_ion_time = self.safe_dose_time(new_ts);

                            let th: f64 = rand::random::<f64>() * 2.0 * PI;
                            let ph: f64 = rand::random::<f64>() * 2.0 * PI;
                            let xn2 = th.sin() * ph.cos();
                            let yn2 = th.sin() * ph.sin();
                            let zn2 = th.cos();

                            if !surrounding {
                                self.add_ionisation(
                                    new_ts,
                                    pixel_coord,
                                    element_name,
                                    new_ion_time,
                                    xn,
                                    yn,
                                    beam,
                                    1,
                                );
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

                            let en = element_name.to_string();
                            self.produce_auger_electron(
                                coef_calc,
                                new_ts,
                                exit_idx as usize,
                                &en,
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
                                &en,
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
                // Simple MC: deposit directly
                self.dose[dt] += shell_binding_energy;
                if idx < self.dose_simple.len() && dt < self.dose_simple[idx].len() {
                    self.dose_simple[idx][dt] += shell_binding_energy;
                }
            }
        } else {
            // Fluorescence
            if let Some(em) = self.fl_data.get(&z) {
                if let Some(td) = em.get(&shell) {
                    let ti = {
                        let rnd2: f64 = rand::random();
                        let mut res = td.cumulative_probs.len().saturating_sub(1);
                        for (k, &p) in td.cumulative_probs.iter().enumerate() {
                            if rnd2 < p {
                                res = k;
                                break;
                            }
                        }
                        res
                    };
                    let fl_energy = if ti < td.energies.len() {
                        td.energies[ti]
                    } else {
                        0.0
                    };
                    let lw = if ti < td.linewidths.len() {
                        td.linewidths[ti]
                    } else {
                        0.0
                    };
                    let drop_idx = if ti < td.drop_index.len() {
                        td.drop_index[ti] as u32
                    } else {
                        0
                    };

                    self.gos_electron_dose_vr_add(dt, shell_binding_energy - fl_energy);
                    if idx < self.voxel_energy_vr.len() && dt < self.voxel_energy_vr[idx].len() {
                        self.voxel_energy_vr[idx][dt] += shell_binding_energy - fl_energy;
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
                }
            }
        }
    }

    // Helper to add to gos_electron_dose_vr at a specific time bin
    fn gos_electron_dose_vr_add(&mut self, dt: usize, val: f64) {
        if dt < self.gos_electron_dose_vr.len() {
            self.gos_electron_dose_vr[dt] += val;
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

        // Exit immediately if energy too low
        if electron_energy < 0.05 {
            let pixel = self.convert_to_pixel_coordinates(prev_x, prev_y, prev_z, angle, wedge);
            let idx = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);
            let dt = self.safe_dose_time(time_stamp);
            if self.is_microcrystal_at(prev_x, prev_y, prev_z, angle, wedge)
                && time_stamp < self.last_time
            {
                if primary_electron {
                    if dt < self.dose.len() {
                        self.dose[dt] += electron_energy;
                    }
                    if dt < self.electron_dose.len() {
                        self.electron_dose[dt] += electron_energy;
                    }
                    if dt < self.gos_electron_dose.len() {
                        self.gos_electron_dose[dt] += electron_energy;
                    }
                    self.gos_electron_dose_vr_add(dt, electron_energy);
                    if idx < self.dose_simple.len() && dt < self.dose_simple[idx].len() {
                        self.dose_simple[idx][dt] += electron_energy;
                    }
                    if idx < self.voxel_energy_vr.len() && dt < self.voxel_energy_vr[idx].len() {
                        self.voxel_energy_vr[idx][dt] += electron_energy;
                    }
                }
                // Low energy ionisations
                let avg_e = coef_calc.plasma_energy(surrounding) / 1000.0;
                if avg_e > 0.0 && !self.simple_mc {
                    let n_ion = (electron_energy / avg_e) as i64;
                    if n_ion > 0 && dt < self.low_energy_ionisations.len() {
                        self.low_energy_ionisations[dt] += n_ion;
                    }
                }
            }
            return;
        }

        let mut stopping_power = coef_calc.stopping_power(electron_energy, surrounding);
        let mut lambda_el = coef_calc.electron_elastic_mfpl(electron_energy, surrounding);
        let mut elastic_probs = coef_calc.elastic_probs(surrounding);
        let _fse_xsection = Self::get_fse_x_section(electron_energy);

        coef_calc.populate_cross_section_coefficients();
        let ionisation_probs = coef_calc.all_shell_probs(surrounding);
        let mut inner_shell_lambda =
            coef_calc.bethe_ionisation_x_section(electron_energy, surrounding);

        let mut gos_inel_lambda = 0.0_f64;
        let mut gos_outer_lambda = 0.0_f64;
        let mut gos_ion_probs: HashMap<String, Vec<f64>> = HashMap::new();
        let mut gos_outer_ion_probs: HashMap<String, f64> = HashMap::new();
        let mut p_inner = 0.0_f64;

        if !self.simple_mc {
            gos_inel_lambda = coef_calc.gos_inel(surrounding, electron_energy);
            let _gos_inner_lambda = coef_calc.gos_inner_lambda(surrounding);
            gos_outer_lambda = coef_calc.gos_outer_lambda(surrounding);
            gos_outer_ion_probs =
                coef_calc.gos_outer_shell_probs_simple(surrounding, gos_outer_lambda);
            gos_ion_probs = coef_calc.gos_shell_probs(surrounding, gos_inel_lambda);
            if inner_shell_lambda > 0.0 {
                gos_inel_lambda = 1.0 / (1.0 / gos_outer_lambda + 1.0 / inner_shell_lambda);
                p_inner = gos_inel_lambda / inner_shell_lambda;
            } else {
                gos_inel_lambda = 1.0 / gos_outer_lambda;
                p_inner = 0.0;
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
                // Entered crystal from surrounding
                if surrounding {
                    entered = true;
                    surrounding = false;
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
                let dose_time_mid = self.safe_dose_time(time_stamp + time_to_dist / 2.0);
                time_stamp += time_to_dist;

                if time_stamp > self.last_time {
                    exited = true;
                    break;
                }

                let dt = dose_time_mid;
                let energy_to_add = energy_lost;

                if primary_electron {
                    let pixel = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
                    let idx2 = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);
                    if dt < self.dose.len() {
                        self.dose[dt] += energy_to_add;
                    }
                    if idx2 < self.dose_simple.len() && dt < self.dose_simple[idx2].len() {
                        self.dose_simple[idx2][dt] += energy_to_add;
                    }
                    if entered {
                        if dt < self.electron_dose_surrounding.len() {
                            self.electron_dose_surrounding[dt] += energy_to_add;
                        }
                    } else if dt < self.electron_dose.len() {
                        self.electron_dose[dt] += energy_to_add;
                    }
                }

                prev_x = xn;
                prev_y = yn;
                prev_z = zn;

                let rnd_inel: f64 = rand::random();
                if rnd_inel < p_inel && !self.simple_mc {
                    // GOS inelastic collision
                    let pixel = self.convert_to_pixel_coordinates(xn, yn, zn, angle, wedge);
                    let idx2 = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);

                    let rnd_inner: f64 = rand::random();
                    let elem_rnd: f64 = rand::random();
                    let mut collided_elem: Option<String> = None;
                    let mut collided_shell: i32 = -1;
                    let mut plasmon = false;

                    if rnd_inner < p_inner {
                        for (name, probs) in &ionisation_probs {
                            let s2 = Self::find_if_element_ionised(probs, elem_rnd);
                            if s2 >= 0 {
                                collided_elem = Some(name.clone());
                                collided_shell = s2;
                                break;
                            }
                        }
                    } else {
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
                        coef_calc.shell_binding_energy(
                            collided_elem.as_deref().unwrap_or(""),
                            collided_shell as usize,
                        )
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

                    let gos_var = coef_calc.gos_variable(surrounding);
                    let col_type = if !plasmon {
                        if let Some(ev) = gos_var.get(collided_elem.as_deref().unwrap_or("")) {
                            let arr: [f64; 4] = [
                                ev.get(collided_shell.max(0) as usize)
                                    .copied()
                                    .unwrap_or(0.0),
                                ev.get(collided_shell.max(0) as usize + 1)
                                    .copied()
                                    .unwrap_or(0.0),
                                ev.get(collided_shell.max(0) as usize + 2)
                                    .copied()
                                    .unwrap_or(0.0),
                                ev.iter().sum(),
                            ];
                            Self::gos_inelastic_type(&arr)
                        } else {
                            2
                        }
                    } else {
                        let pv = coef_calc.plasmon_variable(surrounding);
                        let arr = [pv, pv, pv, pv * 3.0];
                        Self::gos_inelastic_type(&arr)
                    };

                    let w = if col_type == 0 || col_type == 1 {
                        if plasmon {
                            wk / 1000.0
                        } else {
                            Self::get_energy_loss_distant(wdis, uk) / 1000.0
                        }
                    } else {
                        if plasmon {
                            wk / 1000.0
                        } else {
                            let k = Self::sample_k(electron_energy, uk);
                            k * (electron_energy + uk / 1000.0)
                        }
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
                        let dose_time_gos = self.safe_dose_time(time_stamp);
                        if !surrounding {
                            if let Some(ref en) = collided_elem {
                                self.add_ionisation(
                                    time_stamp,
                                    pixel,
                                    en,
                                    dose_time_gos,
                                    xn,
                                    yn,
                                    beam,
                                    1,
                                );
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
                                if dt < self.gos_electron_dose.len() {
                                    self.gos_electron_dose[dt] += w;
                                }
                                self.avg_w += w;
                                self.avg_w_num += 1;
                            }
                        } else {
                            if primary_electron && dt < self.gos_electron_dose.len() {
                                self.gos_electron_dose[dt] += w;
                            }
                            self.gos_electron_dose_vr_add(dt, se_energy);
                            if idx2 < self.voxel_energy_vr.len()
                                && dt < self.voxel_energy_vr[idx2].len()
                            {
                                self.voxel_energy_vr[idx2][dt] += se_energy;
                            }
                            let avg_e = coef_calc.plasma_energy(surrounding) / 1000.0;
                            if avg_e > 0.0 && !surrounding {
                                let n_ion = (se_energy / avg_e) as i64;
                                if n_ion > 0 && dt < self.low_energy_ionisations.len() {
                                    self.low_energy_ionisations[dt] += n_ion;
                                }
                            }
                        }
                        loss_since_update += w;
                        if let Some(ref en) = collided_elem {
                            let en2 = en.clone();
                            self.produce_auger_electron(
                                coef_calc,
                                time_stamp,
                                collided_shell.max(0) as usize,
                                &en2,
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
                    surrounding = true;
                    exited = true; // Don't track in surrounding (simplified)
                } else {
                    // In surrounding, check if worth continuing
                    let max_dist = electron_energy
                        / coef_calc
                            .stopping_power(electron_energy, surrounding)
                            .max(1e-30);
                    if xn.abs() > max_dist + self.x_dimension / 2.0
                        || yn.abs() > max_dist + self.y_dimension / 2.0
                        || zn.abs() > max_dist + self.z_dimension / 2.0
                    {
                        exited = true;
                    } else {
                        let energy_lost = s * stopping_power;
                        let time_to_dist = Self::get_time_to_distance(electron_energy, s);
                        time_stamp += time_to_dist;
                        if time_stamp > self.last_time {
                            exited = true;
                            break;
                        }
                        prev_x = xn;
                        prev_y = yn;
                        prev_z = zn;
                        let def_theta =
                            self.get_scattering_theta(electron_energy, &elastic_probs.clone());
                        let def_phi = Self::get_scattering_phi();
                        let new_vec = Self::get_new_direction_vector(
                            x_norm, y_norm, z_norm, def_theta, def_phi,
                        );
                        x_norm = new_vec[0];
                        y_norm = new_vec[1];
                        z_norm = new_vec[2];
                        electron_energy -= energy_lost;
                        stopping_power = coef_calc.stopping_power(electron_energy, surrounding);
                        lambda_el = coef_calc.electron_elastic_mfpl(electron_energy, surrounding);
                        lambda_t = lambda_el;
                        let r: f64 = rand::random();
                        s = -lambda_t * r.ln();
                        elastic_probs = coef_calc.elastic_probs(surrounding);
                        xn = prev_x + s * x_norm;
                        yn = prev_y + s * y_norm;
                        zn = prev_z + s * z_norm;
                    }
                }
            }

            if electron_energy < 0.05 {
                exited = true;
                let dt2 = self.safe_dose_time(time_stamp);
                if self.is_microcrystal_at(prev_x, prev_y, prev_z, angle, wedge) {
                    let pixel =
                        self.convert_to_pixel_coordinates(prev_x, prev_y, prev_z, angle, wedge);
                    let idx2 = vox3(&self.cryst_size_voxels, pixel[0], pixel[1], pixel[2]);
                    if primary_electron {
                        if dt2 < self.dose.len() {
                            self.dose[dt2] += electron_energy;
                        }
                        if dt2 < self.dose_simple.len().min(idx2 + 1) {
                            // bounds-safe
                        }
                        if idx2 < self.dose_simple.len() && dt2 < self.dose_simple[idx2].len() {
                            self.dose_simple[idx2][dt2] += electron_energy;
                        }
                        if dt2 < self.gos_electron_dose.len() {
                            self.gos_electron_dose[dt2] += electron_energy;
                        }
                        if entered {
                            if dt2 < self.electron_dose_surrounding.len() {
                                self.electron_dose_surrounding[dt2] += electron_energy;
                            }
                        } else if dt2 < self.electron_dose.len() {
                            self.electron_dose[dt2] += electron_energy;
                        }
                    }
                    if !self.simple_mc {
                        self.gos_electron_dose_vr_add(dt2, electron_energy);
                        if idx2 < self.voxel_energy_vr.len()
                            && dt2 < self.voxel_energy_vr[idx2].len()
                        {
                            self.voxel_energy_vr[idx2][dt2] += electron_energy;
                        }
                        let avg_e = coef_calc.plasma_energy(surrounding) / 1000.0;
                        if avg_e > 0.0 {
                            let n_ion = (electron_energy / avg_e) as i64;
                            if n_ion > 0 && dt2 < self.low_energy_ionisations.len() {
                                self.low_energy_ionisations[dt2] += n_ion;
                            }
                        }
                    }
                }
            }
        }
    }

    // ── Dose processing ───────────────────────────────────────────────────────

    pub fn process_dose(&mut self, beam: &dyn Beam, coef_calc: &dyn CoefCalc) {
        let vox_vol_cm3 = self.cryst_pix_per_um.recip().powi(3) * 1.0e-12;
        let _vox_mass_kg = (coef_calc.density() * vox_vol_cm3) / 1000.0;
        let sample_vol_cm3 = self.x_dimension * self.y_dimension * self.z_dimension * 1.0e-21;
        let sample_mass_kg = (coef_calc.density() * sample_vol_cm3) / 1000.0;
        let mean_e_joules = self.mean_energy * KEV_TO_JOULES;
        let n_photons = self.pulse_energy / mean_e_joules;
        let scale = n_photons / self.num_photons as f64;

        self.tot_elastic *= scale;

        let ntb = self.num_time_bins;

        // Scale all time-resolved dose arrays
        for i in 0..ntb {
            if i < self.dose.len() {
                self.dose[i] = ((self.dose[i] * scale * KEV_TO_JOULES) / sample_mass_kg) / 1.0e6;
            }
            if i < self.electron_dose.len() {
                self.electron_dose[i] =
                    ((self.electron_dose[i] * scale * KEV_TO_JOULES) / sample_mass_kg) / 1.0e6;
            }
            if i < self.gos_electron_dose.len() {
                self.gos_electron_dose[i] =
                    ((self.gos_electron_dose[i] * scale * KEV_TO_JOULES) / sample_mass_kg) / 1.0e6;
            }
            if i < self.photon_dose.len() {
                self.photon_dose[i] =
                    ((self.photon_dose[i] * scale * KEV_TO_JOULES) / sample_mass_kg) / 1.0e6;
            }
            if i < self.gos_electron_dose_vr.len() {
                self.gos_electron_dose_vr[i] =
                    ((self.gos_electron_dose_vr[i] * scale * KEV_TO_JOULES) / sample_mass_kg)
                        / 1.0e6;
            }
            if i < self.electron_dose_surrounding.len() {
                self.electron_dose_surrounding[i] =
                    ((self.electron_dose_surrounding[i] * scale * KEV_TO_JOULES) / sample_mass_kg)
                        / 1.0e6;
            }
            if i < self.total_ionisation_events.len() {
                self.total_ionisation_events[i] =
                    (self.total_ionisation_events[i] as f64 * scale).round() as i64;
            }
            if i < self.total_ionisation_events_vr.len() {
                self.total_ionisation_events_vr[i] =
                    (self.total_ionisation_events_vr[i] as f64 * scale).round() as i64;
            }
            if i < self.low_energy_ionisations.len() {
                self.low_energy_ionisations[i] =
                    (self.low_energy_ionisations[i] as f64 * scale).round() as i64;
            }
            // Make cumulative
            if i > 0 {
                let prev_ti = if i - 1 < self.total_ionisation_events.len() {
                    self.total_ionisation_events[i - 1]
                } else {
                    0
                };
                let prev_le = if i - 1 < self.low_energy_ionisations.len() {
                    self.low_energy_ionisations[i - 1]
                } else {
                    0
                };
                if i < self.total_ionisation_events.len() {
                    self.total_ionisation_events[i] += prev_ti;
                }
                if i < self.low_energy_ionisations.len() {
                    self.low_energy_ionisations[i] += prev_le;
                }
            }
        }

        // Scale voxel arrays
        let size = self.cryst_size_voxels;
        for a in 0..size[0] {
            for b in 0..size[1] {
                for c in 0..size[2] {
                    let idx = vox3(&size, a, b, c);
                    if idx < self.voxel_elastic.len() {
                        self.voxel_elastic[idx] *= scale;
                    }
                    if idx < self.voxel_energy_vr.len() {
                        for t in 0..self.voxel_energy_vr[idx].len() {
                            self.voxel_energy_vr[idx][t] *= scale * KEV_TO_JOULES;
                        }
                    }
                    if idx < self.voxel_ionisations_vr.len() {
                        for t in 0..self.voxel_ionisations_vr[idx].len() {
                            self.voxel_ionisations_vr[idx][t] *= scale;
                        }
                    }
                }
            }
        }

        // Scale scalars
        self.raddose_style_dose =
            ((self.raddose_style_dose * scale * KEV_TO_JOULES) / sample_mass_kg) / 1.0e6;
        self.raddose_style_dose_compton =
            ((self.raddose_style_dose_compton * scale * KEV_TO_JOULES) / sample_mass_kg) / 1.0e6;
        self.escaped_energy =
            ((self.escaped_energy * scale * KEV_TO_JOULES) / sample_mass_kg) / 1.0e6;

        // Compute summary metrics
        let exposed_area = self.get_exposed_area(beam);
        let exposed_vol = (exposed_area * self.z_dimension * 1.0e-21).min(sample_vol_cm3);
        let exposed_mass = (coef_calc.density() * exposed_vol) / 1000.0;

        let mut vox_dose_vr = 0.0_f64;
        let mut vox_dose_exposed = 0.0_f64;
        let mut dwd = 0.0_f64;
        let total_elastic = self.tot_elastic;

        for a in 0..size[0] {
            for b in 0..size[1] {
                for c in 0..size[2] {
                    let idx = vox3(&size, a, b, c);
                    if idx >= self.voxel_energy_vr.len() {
                        continue;
                    }
                    let cart = self.convert_to_cartesian_coordinates(a, b, c);
                    let ltv = if c < self.last_time_vox.len() {
                        self.last_time_vox[c]
                    } else {
                        self.last_time
                    };
                    // Sum all time bins within cutoff for this voxel
                    let mut energy_sum_vr = 0.0;
                    for (t, &en) in self.voxel_energy_vr[idx].iter().enumerate() {
                        if (t as f64) * self.pulse_bin_length < ltv - self.pulse_bin_length {
                            energy_sum_vr += en;
                        }
                    }
                    vox_dose_vr += energy_sum_vr;
                    if self.test_if_inside_exposed_area(cart[0], cart[1], beam) {
                        vox_dose_exposed += energy_sum_vr;
                    }
                    let elastic_here = if idx < self.voxel_elastic.len() {
                        self.voxel_elastic[idx]
                    } else {
                        0.0
                    };
                    if elastic_here > 0.0 && total_elastic > 0.0 {
                        dwd += energy_sum_vr * (elastic_here / total_elastic);
                    }
                }
            }
        }

        let vdvr_mgy = if sample_mass_kg > 0.0 {
            (vox_dose_vr / sample_mass_kg) / 1.0e6
        } else {
            0.0
        };
        let vde_mgy = if exposed_mass > 0.0 {
            (vox_dose_exposed / exposed_mass) / 1.0e6
        } else {
            0.0
        };

        let tot_raddose = self.raddose_style_dose + self.raddose_style_dose_compton;
        let rd_exposed = if exposed_mass > 0.0 && sample_mass_kg > 0.0 {
            tot_raddose * sample_mass_kg / exposed_mass
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
        }
        if vde_mgy > 400.0 {
            println!("Warning, damage may begin to be seen at these doses");
        }

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

// ── Data loading (delegates to mc binary loaders) ─────────────────────────────

fn load_xfel_auger_csv(_dir: &str, filename: &str) -> Option<XfelTransitionData> {
    let stem = filename.trim_end_matches(".csv");
    let td = parse_transition_bin(get_transition_bin("auger", stem)?)?;
    Some(XfelTransitionData {
        linewidths: td.linewidths,
        probs: td.probs,
        energies: td.energies,
        cumulative_probs: td.cumulative_probs,
        exit_index: td.exit_index,
        drop_index: td.drop_index,
    })
}

fn load_xfel_fl_csv(_dir: &str, filename: &str) -> Option<XfelTransitionData> {
    let stem = filename.trim_end_matches(".csv");
    let td = parse_transition_bin(get_transition_bin("fl", stem)?)?;
    Some(XfelTransitionData {
        linewidths: td.linewidths,
        probs: td.probs,
        energies: td.energies,
        cumulative_probs: td.cumulative_probs,
        exit_index: td.exit_index,
        drop_index: td.drop_index,
    })
}

fn load_xfel_angle_file(high_energy: bool, atomic_num: usize) -> Option<BTreeMap<u64, Vec<f64>>> {
    let data: &[u8] = if high_energy {
        super::mc::ELSEPA_ABOVE.get(atomic_num)?
    } else {
        super::mc::ELSEPA_BELOW.get(atomic_num)?
    };
    if data.is_empty() {
        return None;
    }
    Some(parse_elsepa_bin(data))
}

fn sample_normal_energies(mean: f64, fwhm: f64, n: u64) -> Vec<f64> {
    crate::energy_distribution::sample_normal_energies(mean, fwhm, n as usize)
}
