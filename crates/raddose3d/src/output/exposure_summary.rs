use std::collections::BTreeMap;

use crate::wedge::Wedge;

/// Generates summary statistics for a single wedge exposure.
/// This is the Rust port of Java's ExposureSummary class.
#[derive(Debug)]
pub struct ExposureSummary {
    /// Red-black tree of dose → count for quantile calculation.
    voxel_doses: BTreeMap<OrderedF64, u32>,

    // Per-voxel exposure variables
    total_absorbed_energy: f64,
    diff_num: f64,
    diff_denom: f64,
    wedge_elastic: f64,

    // Per-image variables
    running_sum_diff_dose: f64,
    images: usize,
    image_exposed_voxels: usize,

    // Summary variables
    total_dose: f64,
    total_energy: f64,
    exposed_voxels: usize,
    occupied_voxels: usize,

    // Results
    avg_diffracted_dose: f64,
    avg_dose_whole_crystal: f64,
    avg_dose_exposed_region: f64,
    used_volume_fraction: f64,
    dose_inefficiency: f64,
    /// Dose inefficiency computed from deposited energy (dose×mass per voxel).
    dose_inefficiency_pe: f64,

    // RDE tracking
    running_sum_rde: f64,
    fluence_weighted_running_sum_rde: f64,
    fluence_sum: f64,
    min_rde: f64,

    // Per-image arrays
    image_dwd: Vec<f64>,
    angle_dwd: Vec<f64>,
    image_vol: Vec<f64>,
    /// Per-image RDE at 5 resolution shells (max-res, 1Å, 2Å, 3Å, 4Å).
    image_rde: Vec<[f64; 5]>,
    /// Per-image [angle_rad, fluence-weighted avg RDE].
    fluence_weighted_rde_array: Vec<[f64; 2]>,
    /// Per-image [angle_rad, min RDE].
    min_rde_array: Vec<[f64; 2]>,

    // Last DWD tracking
    last_dwd_num: f64,
    last_dwd_denom: f64,
    last_dwd_tot: f64,
    last_dwd: f64,

    // Max resolution-related
    de: [f64; 5],
    q: [f64; 5],

    // Dose histogram for final-state reporting (10 bins, 0→30 MGy + overflow).
    /// counts[0] = underflow, counts[1..9] = 0.1…30 MGy, counts[10] = overflow (10+2 = 12 total)
    dose_hist_counts: [u32; 11],
    dose_hist_total: u32,

    // Per-image voxel dose/fluence snapshots (populated during exposure_observation).
    // crystal_size is set at exposure_start; None means no tracking.
    crystal_size: Option<[usize; 3]>,
    /// image_doses[i][flat_vox_idx] = cumulative dose at end of image i.
    pub vox_dose_image: Vec<Vec<f64>>,
    /// image_fluences[i][flat_vox_idx] = cumulative fluence at end of image i.
    pub vox_fluence_image: Vec<Vec<f64>>,
    /// Scratch buffer for current-image doses, flushed into vox_dose_image on image_complete.
    current_image_doses: Vec<f64>,
    /// Scratch buffer for current-image fluences.
    current_image_fluences: Vec<f64>,
}

/// Wrapper for f64 that implements Ord for use in BTreeMap.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderedF64(f64);

impl Eq for OrderedF64 {}

impl PartialOrd for OrderedF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl ExposureSummary {
    const ALPHA: f64 = 1.7;
    const K_CONST: f64 = 81.3;
}

impl Default for ExposureSummary {
    fn default() -> Self {
        Self::new()
    }
}

impl ExposureSummary {
    pub fn new() -> Self {
        let mut voxel_doses = BTreeMap::new();
        voxel_doses.insert(OrderedF64(0.0), 1);

        ExposureSummary {
            voxel_doses,
            total_absorbed_energy: 0.0,
            diff_num: 0.0,
            diff_denom: 0.0,
            wedge_elastic: 0.0,
            running_sum_diff_dose: 0.0,
            images: 0,
            image_exposed_voxels: 0,
            total_dose: 0.0,
            total_energy: 0.0,
            exposed_voxels: 0,
            occupied_voxels: 0,
            avg_diffracted_dose: 0.0,
            avg_dose_whole_crystal: 0.0,
            avg_dose_exposed_region: 0.0,
            used_volume_fraction: 0.0,
            dose_inefficiency: 0.0,
            dose_inefficiency_pe: 0.0,
            running_sum_rde: 0.0,
            fluence_weighted_running_sum_rde: 0.0,
            fluence_sum: 0.0,
            min_rde: 1.0,
            image_dwd: Vec::new(),
            angle_dwd: Vec::new(),
            image_vol: Vec::new(),
            image_rde: Vec::new(),
            fluence_weighted_rde_array: Vec::new(),
            min_rde_array: Vec::new(),
            last_dwd_num: 0.0,
            last_dwd_denom: 0.0,
            last_dwd_tot: 0.0,
            last_dwd: 0.0,
            de: [0.0; 5],
            q: [
                0.0,
                2.0 * std::f64::consts::PI,
                std::f64::consts::PI,
                2.0 * std::f64::consts::PI / 3.0,
                0.5 * std::f64::consts::PI,
            ],
            dose_hist_counts: [0; 11],
            dose_hist_total: 0,
            crystal_size: None,
            vox_dose_image: Vec::new(),
            vox_fluence_image: Vec::new(),
            current_image_doses: Vec::new(),
            current_image_fluences: Vec::new(),
        }
    }

    pub fn exposure_start(&mut self, image_count: usize, wedge: &Wedge, crystal_size: [usize; 3]) {
        self.total_absorbed_energy = 0.0;
        self.diff_num = 0.0;
        self.diff_denom = 0.0;
        self.wedge_elastic = 0.0;
        self.image_exposed_voxels = 0;

        self.image_vol = vec![0.0; image_count];

        self.last_dwd_num = 0.0;
        self.last_dwd_denom = 0.0;
        self.last_dwd_tot = 0.0;

        self.running_sum_rde = 0.0;
        self.fluence_weighted_running_sum_rde = 0.0;
        self.fluence_sum = 0.0;
        self.min_rde = 1.0;

        self.running_sum_diff_dose = 0.0;
        self.images = 0;
        self.image_dwd = vec![0.0; image_count];
        self.angle_dwd = vec![0.0; image_count];
        self.image_rde = vec![[0.0; 5]; image_count];
        self.fluence_weighted_rde_array = vec![[0.0; 2]; image_count];
        self.min_rde_array = vec![[0.0; 2]; image_count];

        // Resolution-dependent De values
        let max_res = wedge.max_resolution;
        self.q[0] = 2.0 * std::f64::consts::PI / max_res;
        for i in 0..5 {
            self.de[i] = Self::K_CONST / self.q[i].powf(Self::ALPHA);
        }

        self.total_dose = 0.0;
        self.total_energy = 0.0;
        self.exposed_voxels = 0;
        self.occupied_voxels = 0;

        self.voxel_doses.clear();
        self.voxel_doses.insert(OrderedF64(0.0), 1);

        // Reset dose histogram
        self.dose_hist_counts = [0; 11];
        self.dose_hist_total = 0;

        // Per-image voxel dose/fluence snapshots
        let nvox = crystal_size[0] * crystal_size[1] * crystal_size[2];
        self.crystal_size = Some(crystal_size);
        self.vox_dose_image = Vec::with_capacity(image_count);
        self.vox_fluence_image = Vec::with_capacity(image_count);
        self.current_image_doses = vec![0.0; nvox];
        self.current_image_fluences = vec![0.0; nvox];
    }

    #[allow(clippy::too_many_arguments)]
    pub fn exposure_observation(
        &mut self,
        wedge_image: usize,
        i: usize,
        j: usize,
        k: usize,
        added_dose: f64,
        total_vox_dose: f64,
        fluence: f64,
        dose_decay: f64,
        absorbed_energy: f64,
        elastic: f64,
        angle_count: f64,
    ) {
        self.diff_num += (total_vox_dose + added_dose / 2.0) * fluence * dose_decay;
        self.diff_denom += fluence * dose_decay;

        if wedge_image == (angle_count as usize).saturating_sub(1) {
            self.last_dwd_num += total_vox_dose * fluence * dose_decay;
            self.last_dwd_denom += fluence * dose_decay;
        }

        if fluence > 0.0 {
            self.image_exposed_voxels += 1;
            self.running_sum_rde += dose_decay;
            self.fluence_sum += fluence;
            self.fluence_weighted_running_sum_rde += fluence * dose_decay;
            if dose_decay < self.min_rde {
                self.min_rde = dose_decay;
            }
        }

        self.total_absorbed_energy += absorbed_energy;
        self.wedge_elastic += elastic;

        // Track per-voxel dose/fluence for this image
        if let Some(size) = self.crystal_size {
            let idx = i * size[1] * size[2] + j * size[2] + k;
            if idx < self.current_image_doses.len() {
                self.current_image_doses[idx] = total_vox_dose + added_dose;
                self.current_image_fluences[idx] = fluence;
            }
        }
    }

    pub fn image_complete(&mut self, image: usize, angle: f64, _last_angle: f64, vox_vol: f64) {
        if self.diff_denom != 0.0 {
            self.running_sum_diff_dose += self.diff_num / self.diff_denom;
            if image < self.image_dwd.len() {
                self.image_dwd[image] = self.diff_num / self.diff_denom;
                // Compute RDE at 5 resolution shells: exp(-DWD / De[i])
                for idx in 0..5 {
                    self.image_rde[image][idx] = (-self.image_dwd[image] / self.de[idx]).exp();
                }
            }
        }

        if image < self.angle_dwd.len() {
            self.angle_dwd[image] = angle;
        }
        if image < self.image_vol.len() {
            self.image_vol[image] = self.image_exposed_voxels as f64 * vox_vol;
        }

        if self.last_dwd_denom != 0.0 {
            self.last_dwd_tot += self.last_dwd_num / self.last_dwd_denom;
        }

        // Fluence-weighted and min RDE per image
        if self.fluence_sum > 0.0 {
            let fw_rde = self.fluence_weighted_running_sum_rde / self.fluence_sum;
            if image < self.fluence_weighted_rde_array.len() {
                self.fluence_weighted_rde_array[image] = [angle, fw_rde];
            }
            if image < self.min_rde_array.len() {
                self.min_rde_array[image] = [angle, self.min_rde];
            }
        }

        // Save per-voxel snapshot for this image
        if self.crystal_size.is_some() {
            let snap_dose = self.current_image_doses.clone();
            let snap_fluence = self.current_image_fluences.clone();
            self.vox_dose_image.push(snap_dose);
            self.vox_fluence_image.push(snap_fluence);
            // Reset scratch buffers (clear, keeping allocation)
            for v in &mut self.current_image_doses {
                *v = 0.0;
            }
            for v in &mut self.current_image_fluences {
                *v = 0.0;
            }
        }

        // Reset per-image state
        self.running_sum_rde = 0.0;
        self.fluence_weighted_running_sum_rde = 0.0;
        self.image_exposed_voxels = 0;
        self.fluence_sum = 0.0;
        self.min_rde = 1.0;
        self.diff_num = 0.0;
        self.diff_denom = 0.0;
        self.last_dwd_num = 0.0;
        self.last_dwd_denom = 0.0;
        self.images += 1;
    }

    pub fn summary_observation(
        &mut self,
        _i: usize,
        _j: usize,
        _k: usize,
        voxel_dose: f64,
        voxel_mass_kg: f64,
    ) {
        self.occupied_voxels += 1;

        if voxel_dose > 0.0 {
            let entry = self.voxel_doses.entry(OrderedF64(voxel_dose)).or_insert(0);
            *entry += 1;

            self.total_dose += voxel_dose;
            self.total_energy += (voxel_dose * 1e6) * voxel_mass_kg;
            self.exposed_voxels += 1;

            // Accumulate dose histogram (9 bins from 0.1 to 30.0 MGy)
            const HIST_MIN: f64 = 0.1;
            const HIST_MAX: f64 = 30.0;
            const HIST_BINS: usize = 9;
            const HIST_STEP: f64 = (HIST_MAX - HIST_MIN) / HIST_BINS as f64;
            let bucket = if voxel_dose < HIST_MIN {
                0usize
            } else if voxel_dose >= HIST_MAX {
                HIST_BINS + 1
            } else {
                1 + ((voxel_dose - HIST_MIN) / HIST_STEP) as usize
            };
            if bucket < self.dose_hist_counts.len() {
                self.dose_hist_counts[bucket] += 1;
            }
            self.dose_hist_total += 1;
        }
    }

    pub fn exposure_complete(&mut self) {
        if self.images > 0 {
            self.avg_diffracted_dose = self.running_sum_diff_dose / self.images as f64;
        }

        if self.exposed_voxels > 0 {
            self.avg_dose_exposed_region = self.total_dose / self.exposed_voxels as f64;
        }

        if self.occupied_voxels > 0 {
            self.used_volume_fraction =
                100.0 * self.exposed_voxels as f64 / self.occupied_voxels as f64;
            self.avg_dose_whole_crystal = self.total_dose / self.occupied_voxels as f64;
        }

        if self.total_absorbed_energy > 0.0 {
            self.dose_inefficiency =
                (self.max_dose() * 1e6) / (1000.0 * self.total_absorbed_energy);
        }

        if self.total_energy > 0.0 {
            self.dose_inefficiency_pe = (self.max_dose() * 1e6) / (1000.0 * self.total_energy);
        }

        if self.images > 0 {
            self.last_dwd = self.image_dwd[self.images - 1];
        }
    }

    // --- Getters ---

    pub fn avg_diffracted_dose(&self) -> f64 {
        self.avg_diffracted_dose
    }

    pub fn avg_dose_whole_crystal(&self) -> f64 {
        self.avg_dose_whole_crystal
    }

    pub fn avg_dose_exposed_region(&self) -> f64 {
        self.avg_dose_exposed_region
    }

    pub fn max_dose(&self) -> f64 {
        self.voxel_doses
            .keys()
            .next_back()
            .map(|k| k.0)
            .unwrap_or(0.0)
    }

    pub fn total_dose(&self) -> f64 {
        self.total_dose
    }

    pub fn total_energy(&self) -> f64 {
        self.total_energy
    }

    pub fn used_volume_fraction(&self) -> f64 {
        self.used_volume_fraction
    }

    pub fn abs_energy_total(&self) -> f64 {
        self.total_absorbed_energy
    }

    pub fn wedge_elastic(&self) -> f64 {
        self.wedge_elastic
    }

    pub fn dose_inefficiency(&self) -> f64 {
        self.dose_inefficiency
    }

    pub fn dose_inefficiency_pe(&self) -> f64 {
        self.dose_inefficiency_pe
    }

    pub fn image_dwds(&self) -> &[f64] {
        &self.image_dwd
    }

    pub fn last_dwd(&self) -> f64 {
        self.last_dwd
    }

    /// Get the absolute dose threshold for a given quantile.
    pub fn abs_dose_threshold(&self, dose_quantile: f64) -> f64 {
        let dose_cutoff = (1.0 - dose_quantile) * self.total_dose;
        let mut dose_seen = 0.0;

        for (dose_val, &count) in &self.voxel_doses {
            dose_seen += dose_val.0 * count as f64;
            if dose_seen >= dose_cutoff {
                return dose_val.0;
            }
        }
        f64::INFINITY
    }

    /// Average dose within a dose quantile volume.
    pub fn avg_dose_threshold(&self, dose_quantile: f64) -> f64 {
        let threshold = self.abs_dose_threshold(dose_quantile);
        let mut voxels_above = 0u32;
        for (dose_val, &count) in &self.voxel_doses {
            if dose_val.0 >= threshold {
                voxels_above += count;
            }
        }
        if voxels_above > 0 {
            dose_quantile * self.total_dose / voxels_above as f64
        } else {
            0.0
        }
    }

    /// Dose contrast: max dose / average dose in quantile volume.
    pub fn dose_contrast(&self, dose_quantile: f64) -> f64 {
        let avg = self.avg_dose_threshold(dose_quantile);
        if avg > 0.0 {
            self.max_dose() / avg
        } else {
            0.0
        }
    }

    pub fn occupied_voxels(&self) -> usize {
        self.occupied_voxels
    }

    pub fn exposed_voxels(&self) -> usize {
        self.exposed_voxels
    }

    /// Per-image DWD values (one per image).
    pub fn image_dwd_array(&self) -> &[[f64; 5]] {
        &self.image_rde
    }

    /// Per-image angle in radians.
    pub fn angle_dwd_array(&self) -> &[f64] {
        &self.angle_dwd
    }

    /// Per-image exposed volume (µm³).
    pub fn image_vol_array(&self) -> &[f64] {
        &self.image_vol
    }

    /// Per-image RDE at 5 resolution shells: [max_res, 1Å, 2Å, 3Å, 4Å].
    pub fn image_rde_array(&self) -> &[[f64; 5]] {
        &self.image_rde
    }

    /// Per-image [angle_rad, fluence-weighted avg RDE].
    pub fn fluence_weighted_rde_array(&self) -> &[[f64; 2]] {
        &self.fluence_weighted_rde_array
    }

    /// Per-image [angle_rad, min RDE].
    pub fn min_rde_array(&self) -> &[[f64; 2]] {
        &self.min_rde_array
    }

    /// Crystal size set at last exposure_start, or None.
    pub fn crystal_size(&self) -> Option<[usize; 3]> {
        self.crystal_size
    }

    /// Normalised dose histogram fractions (underflow, 9 bins, overflow = 11 entries).
    /// Bins span 0.1 → 30.0 MGy in equal steps.
    pub fn dose_hist_normalised(&self) -> [f64; 11] {
        let mut out = [0.0f64; 11];
        if self.dose_hist_total > 0 {
            let total = self.dose_hist_total as f64;
            for (i, &c) in self.dose_hist_counts.iter().enumerate() {
                out[i] = c as f64 / total;
            }
        }
        out
    }

    /// Dose histogram bin boundary for index i (1-based; i=0 is -∞, i=9 is 30.0).
    pub fn dose_hist_break(i: usize) -> f64 {
        const HIST_MIN: f64 = 0.1;
        const HIST_STEP: f64 = (30.0 - HIST_MIN) / 9.0;
        HIST_MIN + i as f64 * HIST_STEP
    }
}
