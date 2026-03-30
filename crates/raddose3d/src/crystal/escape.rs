//! Photoelectron (PE) and fluorescent (FL) escape from crystals.
//!
//! Ported from Java CrystalPolyhedron.java (addDoseAfterPE, addDoseAfterFL,
//! setPEparamsForCurrentBeam, setFLparamsForCurrentBeam, and helpers).

#[allow(unused_imports)]
use rand::thread_rng;
use rand::Rng;
use std::f64::consts::PI;

use crate::coefcalc::CoefCalc;
use crate::constants::KEV_TO_JOULES;

// ── Constants ─────────────────────────────────────────────────────────────────

/// Angular resolution limit for PE tracks (matches Java PE_ANGLE_RES_LIMIT).
const PE_ANGLE_RES_LIMIT: usize = 6;
/// Angular limit (radians) for PE tracks (matches Java PE_ANGLE_LIMIT).
const PE_ANGLE_LIMIT: f64 = PI;
/// Number of random tracks sampled per voxel in addDoseAfterPE.
const PE_ANGLE_RESOLUTION: usize = 1;
/// Public alias for use in mod.rs.
pub const PE_ANGLE_RESOLUTION_PUB: usize = PE_ANGLE_RESOLUTION;
/// Angular resolution limit for FL tracks (matches Java FL_ANGLE_RES_LIMIT).
const FL_ANGLE_RES_LIMIT: usize = 6;
/// Number of FL angular samples per voxel in addDoseAfterFL.
const FL_ANGLE_RESOLUTION: usize = 1;
/// Public alias for use in mod.rs.
pub const FL_ANGLE_RESOLUTION_PUB: usize = FL_ANGLE_RESOLUTION;
/// Number of FL distance bins (matches Java default of 8).
const FL_DIST_BINS: usize = 8;
/// Degree of beam polarisation (matches Java).
const DEGREE_OF_POLARISATION: f64 = 0.75;

// ── Energy coefficient data ───────────────────────────────────────────────────

const ENERGY_COEFS_LOW: &str = include_str!("../constants_data/EnergyCoefsLow.txt");
const ENERGY_COEFS_MED: &str = include_str!("../constants_data/EnergyCoefsMed.txt");
const ENERGY_COEFS_HIGH: &str = include_str!("../constants_data/EnergyCoefsHigh.txt");

// ── PE Escape ─────────────────────────────────────────────────────────────────

/// Pre-computed photoelectron escape parameters.
#[derive(Debug)]
pub struct PeEscape {
    /// Number of distance bins.
    pub pe_dist_bins: usize,
    /// Distance values in µm for each bin.
    pub pe_distances: Vec<f64>,
    /// Fraction of PE dose deposited at each distance bin.
    pub propn_dose_at_dist: Vec<f64>,
    /// Pre-computed relative voxel offsets: [dist_bin][track] -> [x,y,z].
    pub relative_vox_xyz: Vec<Vec<[f64; 3]>>,
    /// Track indices biased by angular distribution (beam polarization).
    pub track_bias: Vec<usize>,
    /// Binding energy to subtract from PE energy (keV).
    pub energy_to_subtract: f64,
    /// Auger energy per photoelectric event (Joules).
    pub auger_energy: f64,
    /// Total FL energy release per photon (Joules). Used to partition dose.
    pub fl_energy_release: f64,
    /// Fluorescence proportion per event: [element][K, L1, L2, L3].
    pub fl_proportion_event: Vec<[f64; 4]>,
    /// Total PE escaped dose tracking.
    pub total_escaped_pe: f64,
    /// Total FL escaped dose tracking.
    pub total_escaped_fl: f64,
}

/// Pre-computed fluorescent escape parameters.
#[derive(Debug)]
pub struct FlEscape {
    /// Number of distance bins.
    pub fl_dist_bins: usize,
    /// Distance distribution: [element][shell][dist_bin] — absorption probability.
    pub fl_dist_distribution: Vec<Vec<Vec<f64>>>,
    /// Distances travelled: [element][shell][dist_bin] in µm.
    pub fl_distances: Vec<Vec<Vec<f64>>>,
    /// Relative voxel offsets: [element][shell][dist_bin][angle] -> [x,y,z].
    pub fl_relative_vox_xyz: Vec<Vec<Vec<Vec<[f64; 3]>>>>,
    /// Number of angular tracks for FL.
    pub num_fl_tracks: usize,
}

/// Pre-computed cryo (surrounding) PE parameters.
#[derive(Debug)]
pub struct CryoEscape {
    /// Cryo PE distance bins.
    pub cryo_pe_distances: Vec<f64>,
    /// Cryo dose distribution along PE tracks.
    pub propn_dose_at_dist_cryo: Vec<f64>,
    /// Relative voxel offsets (cryo → crystal coords): [dist][track] -> [x,y,z].
    pub relative_vox_xyz_cryo: Vec<Vec<[f64; 3]>>,
    /// Track bias for cryo PE.
    pub cryo_track_bias: Vec<usize>,
    /// Cryo binding energy to subtract (keV).
    pub cryo_energy_to_subtract: f64,
    /// Cryo auger energy per event (Joules).
    pub cryo_auger_energy: f64,
    /// Cryo/crystal density blending factor.
    pub cryo_and_crystal_density: f64,
    /// Extra voxels of surrounding around crystal per side.
    pub cryo_extra_voxels: usize,
    /// Cryo voxel grid dimensions [nx, ny, nz].
    pub cryo_size_voxels: [usize; 3],
    /// Cryo pixels per µm (may differ from crystal).
    pub cryo_pix_per_um: f64,
    /// Cryo FL energy release per photon (Joules).
    pub cryo_fl_energy_release: f64,
    /// Accumulated dose from surrounding.
    pub total_dose_from_surrounding: f64,
}

// ── Setup functions ───────────────────────────────────────────────────────────

/// Calculate PE escape binding energy to subtract from PE energy.
/// Matches Java Crystal.calculatePEEnergySubtraction().
fn calculate_pe_energy_subtraction(fe_factors: &[Vec<f64>]) -> f64 {
    // Java layout: [0]=mu_ratio, [1]=edge, [2]=ion_prob, [3]=fl_yield, [4]=muabs, ...
    // Rust layout: [0]=mu_ratio, [1]=ion_prob, [2]=fl_yield, [3]=edge, [4]=muabs, ...
    //
    // Java calculates: sum(mu_ratio * edge * ion_prob) for each shell.
    // In Rust indices: edge=row[3], ion_prob=row[1] for K shell.
    let mut total = 0.0;
    for row in fe_factors {
        let mu_ratio = row[0];
        // K: edge[3] * ion_prob[1]
        total += mu_ratio * safe_idx(row, 3) * safe_idx(row, 1);
        // L1: edge[7] * ion_prob[5]
        total += mu_ratio * safe_idx(row, 7) * safe_idx(row, 5);
        // L2: edge[11] * ion_prob[9]
        total += mu_ratio * safe_idx(row, 11) * safe_idx(row, 9);
        // L3: edge[15] * ion_prob[13]
        total += mu_ratio * safe_idx(row, 15) * safe_idx(row, 13);
        // M shells: edge[18] * ion_prob[17], etc.
        total += mu_ratio * safe_idx(row, 18) * safe_idx(row, 17);
        total += mu_ratio * safe_idx(row, 20) * safe_idx(row, 19);
        total += mu_ratio * safe_idx(row, 22) * safe_idx(row, 21);
        total += mu_ratio * safe_idx(row, 24) * safe_idx(row, 23);
        total += mu_ratio * safe_idx(row, 26) * safe_idx(row, 25);
    }
    total
}

/// Calculate Auger energy per photoelectric event (in Joules).
/// Matches Java Crystal.getAugerEnergy().
fn calculate_auger_energy(fe_factors: &[Vec<f64>]) -> f64 {
    let mut auger = 0.0;
    for row in fe_factors {
        let mu_ratio = row[0];
        // K shell only: ion_prob * fl_yield * (1 - edge) * edge
        // Java: K1 = feFactors[i][1] * feFactors[i][2] * (1 - feFactors[i][3])
        // This is K_ion_prob * K_fl_yield * (1 - K_edge_energy)
        // Wait, that doesn't make physical sense with energy...
        // Let me re-read Java:
        //   double K1 = fluorescentEscapeFactors[i][1] * fluorescentEscapeFactors[i][2] *
        //               (1 - fluorescentEscapeFactors[i][3]);
        // [1] = K_ion_prob, [2] = K_edge_energy (in this context), [3] = ???
        // No wait - [2] = K_fl_yield, [3] = K_edge
        // K1 = K_ion_prob * K_fl_yield * (1 - K_edge_energy)
        // That still seems odd... let me look at Auger physics:
        // After K-shell ionization, either fluorescence (prob = fl_yield) or Auger (prob = 1 - fl_yield)
        // Auger energy ≈ edge energy (the ionization energy is released as Auger electron)
        // So: K_auger = ion_prob * (1 - fl_yield) * edge_energy
        // Java: K1 = [1]*[2]*(1-[3]) = ion_prob * fl_yield * (1 - edge)
        // That looks wrong... unless the Java fe_factors layout is different from what I think.

        // Let me check the Java getAugerEnergy more carefully. In Crystal.java:
        // double K1 = fluorescentEscapeFactors[i][1] * fluorescentEscapeFactors[i][2] *
        //             (1 - fluorescentEscapeFactors[i][3]);
        // Then augerEnergy += K1 * muratio;
        // And augerEnergy = augerEnergy * Beam.KEVTOJOULES;

        // In CoefCalcCompute.getFluorescentEscapeFactors (Java):
        // factors[elementCounter][1] = ionisationProb (K shell)
        // factors[elementCounter][2] = edge energy (K shell)  -- NOT fl_yield!
        // factors[elementCounter][3] = fluorescenceYield (K shell)
        // factors[elementCounter][4] = muabs at fluorescence energy

        // Wait, the Java layout is DIFFERENT from the Rust layout!
        // Let me check the Rust calc_fluorescent_escape_factors:
        // row[1] = K ionisation prob
        // row[2] = K fluorescence yield
        // row[3] = K edge energy
        // row[4] = K fl photoelectric coeff

        // Java CoefCalcCompute.getFluorescentEscapeFactors:
        // factors[el][1] = ionisationProb
        // factors[el][2] = edgeEnergy
        // factors[el][3] = fluorescenceYield
        // factors[el][4] = muabs (photoelectric coeff at fluorescence energy)

        // So Java [2] = edge energy, Rust [2] = fl_yield. They're SWAPPED!
        // Java [3] = fl_yield, Rust [3] = edge energy. They're SWAPPED!

        // This means for Java getAugerEnergy:
        // K1 = ion_prob * edge_energy * (1 - fl_yield)
        // = probability of K-shell Auger * edge energy
        // augerEnergy = sum(mu_ratio * ion_prob * edge * (1-fl_yield)) * KEV_TO_JOULES
        // This makes physical sense!

        // For Rust: edge is at [3], fl_yield is at [2]
        // So: K_auger = ion_prob[1] * edge[3] * (1 - fl_yield[2])
        let k_auger = safe_idx(row, 1) * safe_idx(row, 3) * (1.0 - safe_idx(row, 2));
        auger += k_auger * mu_ratio;
    }
    auger * KEV_TO_JOULES
}

/// Calculate total fluorescent energy release per photon absorbed (Joules).
/// Matches Java Crystal.calcFluorescence().
fn calc_fluorescence_energy(fe_factors: &[Vec<f64>]) -> f64 {
    let mut total = 0.0;
    for row in fe_factors {
        let mu_ratio = row[0];
        // Java: K1 = [1]*[2]*[3] = ion_prob * edge * fl_yield
        // Rust: ion_prob[1] * fl_yield[2] * edge[3]
        // Same product, different ordering
        let k = safe_idx(row, 1) * safe_idx(row, 2) * safe_idx(row, 3);
        // L1: ion_prob[5] * fl_yield[6] * edge[7]
        let l1 = safe_idx(row, 5) * safe_idx(row, 6) * safe_idx(row, 7);
        // L2: ion_prob[9] * fl_yield[10] * edge[11]
        let l2 = safe_idx(row, 9) * safe_idx(row, 10) * safe_idx(row, 11);
        // L3: ion_prob[13] * fl_yield[14] * edge[15]
        let l3 = safe_idx(row, 13) * safe_idx(row, 14) * safe_idx(row, 15);
        total += (k + l1 + l2 + l3) * mu_ratio;
    }
    total * KEV_TO_JOULES
}

/// Calculate fluorescent energy per event for each element and shell (Joules).
/// Matches Java Crystal.getFluorescenceEnergyPerEvent().
fn calc_fl_energy_per_event(fe_factors: &[Vec<f64>]) -> Vec<[f64; 4]> {
    let mut result = Vec::with_capacity(fe_factors.len());
    for row in fe_factors {
        let mu_ratio = row[0];
        let k = safe_idx(row, 1) * safe_idx(row, 2) * safe_idx(row, 3) * mu_ratio * KEV_TO_JOULES;
        let l1 = safe_idx(row, 5) * safe_idx(row, 6) * safe_idx(row, 7) * mu_ratio * KEV_TO_JOULES;
        let l2 =
            safe_idx(row, 9) * safe_idx(row, 10) * safe_idx(row, 11) * mu_ratio * KEV_TO_JOULES;
        let l3 =
            safe_idx(row, 13) * safe_idx(row, 14) * safe_idx(row, 15) * mu_ratio * KEV_TO_JOULES;
        result.push([k, l1, l2, l3]);
    }
    result
}

/// Set up PE escape parameters for the current beam energy.
/// Matches Java CrystalPolyhedron.setPEparamsForCurrentBeam().
pub fn setup_pe_escape(
    beam_energy: f64,
    coefcalc: &dyn CoefCalc,
    fe_factors: &[Vec<f64>],
    pix_per_um: f64,
    _cryst_size_um: [f64; 3],
    pe_resolution: Option<i32>,
) -> PeEscape {
    // 1. Calculate binding energy subtraction
    let energy_to_subtract = calculate_pe_energy_subtraction(fe_factors);
    let pe_energy = beam_energy - energy_to_subtract;

    // 2. Calculate auger energy
    let auger_energy = calculate_auger_energy(fe_factors);

    // 3. Calculate FL energy release and proportions
    let fl_energy_release = calc_fluorescence_energy(fe_factors);
    let fl_energy_per_event = calc_fl_energy_per_event(fe_factors);
    let fl_proportion_event: Vec<[f64; 4]> = if fl_energy_release > 0.0 {
        fl_energy_per_event
            .iter()
            .map(|e| {
                [
                    e[0] / fl_energy_release,
                    e[1] / fl_energy_release,
                    e[2] / fl_energy_release,
                    e[3] / fl_energy_release,
                ]
            })
            .collect()
    } else {
        vec![[0.0; 4]; fe_factors.len()]
    };

    // 4. Set Gumbel distribution parameters based on density and PE energy
    let density = coefcalc.density();
    let gumbel_loc = set_gumbel_loc(density, pe_energy);
    let gumbel_scale = set_gumbel_scale(density, pe_energy);

    // 5. Calculate max PE distance
    let max_pe_dist = get_max_pe_distance(pe_energy, &gumbel_loc, &gumbel_scale);

    // 6. Set up distance bins
    let pe_res = pe_resolution.unwrap_or(0);
    let pe_dist_bins = calc_pe_dist_bins(max_pe_dist, pix_per_um, pe_res);
    let bin_interval = max_pe_dist / (pe_dist_bins - 1).max(1) as f64;
    let pe_distances: Vec<f64> = (0..pe_dist_bins).map(|i| i as f64 * bin_interval).collect();

    // 7. Calculate dose deposition profile along PE tracks
    let propn_dose_at_dist =
        calc_propn_dose_deposited(&pe_distances, pe_energy, &gumbel_loc, &gumbel_scale);

    // 8. Set up angular distribution (beam polarization bias)
    let angular_distribution = setup_pe_polarisation(coefcalc, beam_energy, fe_factors);

    // 9. Find voxels reached by PE
    let (relative_vox_xyz, track_bias) =
        find_voxels_reached_by_pe(&pe_distances, pix_per_um, &angular_distribution);

    PeEscape {
        pe_dist_bins,
        pe_distances,
        propn_dose_at_dist,
        relative_vox_xyz,
        track_bias,
        energy_to_subtract,
        auger_energy,
        fl_energy_release,
        fl_proportion_event,
        total_escaped_pe: 0.0,
        total_escaped_fl: 0.0,
    }
}

/// Set up fluorescent escape parameters.
/// Matches Java CrystalPolyhedron.setFLparamsForCurrentBeam() + calcFluorescenceDistribution().
pub fn setup_fl_escape(
    fe_factors: &[Vec<f64>],
    pix_per_um: f64,
    cryst_size_um: [f64; 3],
) -> FlEscape {
    let fl_dist_bins = FL_DIST_BINS;

    // Calculate distance distribution and distances for each element/shell
    let (fl_dist_distribution, fl_distances) =
        calc_fluorescence_distribution(fe_factors, fl_dist_bins, cryst_size_um);

    // Find voxels reached by FL
    let (fl_relative_vox_xyz, num_fl_tracks) =
        find_voxels_reached_by_fl(fe_factors, &fl_distances, fl_dist_bins, pix_per_um);

    FlEscape {
        fl_dist_bins,
        fl_dist_distribution,
        fl_distances,
        fl_relative_vox_xyz,
        num_fl_tracks,
    }
}

// ── Per-voxel escape application ──────────────────────────────────────────────

/// Apply PE escape to a voxel's dose. Redistributes PE dose along random tracks;
/// adds dose to target voxels inside crystal, tracks dose that escapes.
/// Returns the dose that escaped the crystal.
///
/// Matches Java CrystalPolyhedron.addDoseAfterPE().
pub fn apply_pe_escape(
    pe: &PeEscape,
    crystal_size: [usize; 3],
    is_crystal_at: &dyn Fn(usize, usize, usize) -> bool,
    add_dose: &mut dyn FnMut(usize, usize, usize, f64),
    i: usize,
    j: usize,
    k: usize,
    dose_pe: f64,
) -> f64 {
    let mut rng = rand::thread_rng();
    let mut dose_lost = 0.0;

    if pe.track_bias.is_empty() {
        return dose_lost;
    }

    for _ in 0..(PE_ANGLE_RESOLUTION * PE_ANGLE_RESOLUTION) {
        let random_idx = rng.gen_range(0..pe.track_bias.len());
        let random_track = pe.track_bias[random_idx];

        for m in 0..pe.pe_dist_bins {
            if random_track >= pe.relative_vox_xyz[m].len() {
                continue;
            }
            let [rx, ry, rz] = pe.relative_vox_xyz[m][random_track];

            let partial_dose = dose_pe * pe.propn_dose_at_dist[m]
                / (PE_ANGLE_RESOLUTION * PE_ANGLE_RESOLUTION) as f64;

            let ti = (i as f64 + rx).round() as i64;
            let tj = (j as f64 + ry).round() as i64;
            let tk = (k as f64 + rz).round() as i64;

            if ti >= 0
                && tj >= 0
                && tk >= 0
                && (ti as usize) < crystal_size[0]
                && (tj as usize) < crystal_size[1]
                && (tk as usize) < crystal_size[2]
                && is_crystal_at(ti as usize, tj as usize, tk as usize)
            {
                add_dose(ti as usize, tj as usize, tk as usize, partial_dose);
            } else {
                dose_lost += partial_dose;
            }
        }
    }

    dose_lost
}

/// Apply FL escape to a voxel's dose. Redistributes FL dose along isotropic tracks;
/// adds dose to target voxels inside crystal, tracks dose that escapes.
/// Returns the dose that escaped the crystal.
///
/// Matches Java CrystalPolyhedron.addDoseAfterFL().
pub fn apply_fl_escape(
    pe: &PeEscape,
    fl: &FlEscape,
    crystal_size: [usize; 3],
    is_crystal_at: &dyn Fn(usize, usize, usize) -> bool,
    add_dose: &mut dyn FnMut(usize, usize, usize, f64),
    i: usize,
    j: usize,
    k: usize,
    dose_fl: f64,
) -> f64 {
    let mut dose_lost = 0.0;

    let num_elements = pe.fl_proportion_event.len();

    for n in 0..num_elements {
        for l in 0..4 {
            // K, L1, L2, L3
            for m in 0..fl.fl_dist_bins {
                for q in 0..(FL_ANGLE_RESOLUTION * FL_ANGLE_RESOLUTION) {
                    if pe.fl_proportion_event[n][l] == 0.0 {
                        continue;
                    }

                    let fl_partial_dose =
                        dose_fl * pe.fl_proportion_event[n][l] * fl.fl_dist_distribution[n][l][m]
                            / (FL_ANGLE_RESOLUTION * FL_ANGLE_RESOLUTION) as f64;

                    if n >= fl.fl_relative_vox_xyz.len()
                        || l >= fl.fl_relative_vox_xyz[n].len()
                        || m >= fl.fl_relative_vox_xyz[n][l].len()
                        || q >= fl.fl_relative_vox_xyz[n][l][m].len()
                    {
                        continue;
                    }

                    let [rx, ry, rz] = fl.fl_relative_vox_xyz[n][l][m][q];

                    let ti = (i as f64 + rx).round() as i64;
                    let tj = (j as f64 + ry).round() as i64;
                    let tk = (k as f64 + rz).round() as i64;

                    if ti >= 0
                        && tj >= 0
                        && tk >= 0
                        && (ti as usize) < crystal_size[0]
                        && (tj as usize) < crystal_size[1]
                        && (tk as usize) < crystal_size[2]
                        && is_crystal_at(ti as usize, tj as usize, tk as usize)
                    {
                        add_dose(ti as usize, tj as usize, tk as usize, fl_partial_dose);
                    } else {
                        dose_lost += fl_partial_dose;
                    }
                }
            }
        }
    }

    dose_lost
}

// ── Cryo (surrounding) PE escape ──────────────────────────────────────────────

/// Compute cryo grid pixels-per-micron, matching Java's `setCryoPPM()`.
/// Java uses a separate (often coarser) resolution for the surrounding grid.
fn set_cryo_ppm(
    max_pe_distance: f64,
    beam_min_dim: f64,
    cryst_size_um: [f64; 3],
    crystal_pix_per_um: f64,
) -> f64 {
    let min_dim_cryst = cryst_size_um[0].min(cryst_size_um[1]).min(cryst_size_um[2]);

    let mut multiply_factor = 20i32;
    if min_dim_cryst >= beam_min_dim {
        multiply_factor = (multiply_factor as f64 * 1.5) as i32; // Java: *= 1.5 on int
    }

    // MCSim is always false in normal exposure path
    let ideal_ppm = (1.0 / max_pe_distance) * 5.0
        + (1.0 / max_pe_distance) * multiply_factor as f64 * (1.0 / min_dim_cryst);

    if ideal_ppm >= crystal_pix_per_um {
        (ideal_ppm / crystal_pix_per_um).ceil() * crystal_pix_per_um
    } else {
        crystal_pix_per_um / (crystal_pix_per_um / ideal_ppm) as i32 as f64
    }
}

/// Set up cryo (surrounding) PE escape parameters.
pub fn setup_cryo_escape(
    beam_energy: f64,
    coefcalc: &dyn CoefCalc,
    cryo_fe_factors: &[Vec<f64>],
    crystal_pix_per_um: f64,
    cryst_size_um: [f64; 3],
    fl_enabled: bool,
    beam_min_dim: f64,
) -> CryoEscape {
    let cryo_energy_to_subtract = calculate_pe_energy_subtraction(cryo_fe_factors);
    let cryo_auger_energy = calculate_auger_energy(cryo_fe_factors);
    let cryo_fl_energy_release = if fl_enabled {
        calc_fluorescence_energy(cryo_fe_factors)
    } else {
        0.0
    };
    let pe_energy = beam_energy - cryo_energy_to_subtract;

    // Gumbel params for cryo density
    let cryo_density = coefcalc.cryo_density();
    let gumbel_loc = set_gumbel_loc(cryo_density, pe_energy);
    let gumbel_scale = set_gumbel_scale(cryo_density, pe_energy);

    let max_pe_dist = get_max_pe_distance(pe_energy, &gumbel_loc, &gumbel_scale);

    // Cryo grid uses its own PPM, matching Java's setCryoPPM()
    let cryo_ppm = set_cryo_ppm(max_pe_dist, beam_min_dim, cryst_size_um, crystal_pix_per_um);

    // Calculate cryo/crystal density blending
    let average_side = (cryst_size_um[0] + cryst_size_um[1] + cryst_size_um[2]) / 3.0;
    let cryo_and_crystal_density = if average_side >= max_pe_dist {
        0.5
    } else {
        max_pe_dist / (max_pe_dist + average_side)
    };

    // Extra voxels for surrounding (Java uses cryo PPM here)
    let cryo_extra_voxels = (max_pe_dist * cryo_ppm).ceil() as usize;
    let crystal_size_voxels = [
        (cryst_size_um[0] * cryo_ppm).round() as usize + 1,
        (cryst_size_um[1] * cryo_ppm).round() as usize + 1,
        (cryst_size_um[2] * cryo_ppm).round() as usize + 1,
    ];
    let cryo_size_voxels = [
        crystal_size_voxels[0] + 2 * cryo_extra_voxels,
        crystal_size_voxels[1] + 2 * cryo_extra_voxels,
        crystal_size_voxels[2] + 2 * cryo_extra_voxels,
    ];

    // Reuse same PE distance bins as crystal
    let pe_dist_bins = calc_pe_dist_bins(max_pe_dist, cryo_ppm, 0);
    let bin_interval = max_pe_dist / (pe_dist_bins - 1).max(1) as f64;
    let cryo_pe_distances: Vec<f64> = (0..pe_dist_bins).map(|i| i as f64 * bin_interval).collect();

    let propn_dose_at_dist_cryo =
        calc_propn_dose_deposited(&cryo_pe_distances, pe_energy, &gumbel_loc, &gumbel_scale);

    // Angular distribution for cryo
    let cryo_angular_distribution = setup_pe_polarisation(coefcalc, beam_energy, cryo_fe_factors);

    let (relative_vox_xyz_cryo, cryo_track_bias) =
        find_voxels_reached_by_pe(&cryo_pe_distances, cryo_ppm, &cryo_angular_distribution);

    CryoEscape {
        cryo_pe_distances,
        propn_dose_at_dist_cryo,
        relative_vox_xyz_cryo,
        cryo_track_bias,
        cryo_energy_to_subtract,
        cryo_auger_energy,
        cryo_and_crystal_density,
        cryo_extra_voxels,
        cryo_size_voxels,
        cryo_pix_per_um: cryo_ppm,
        cryo_fl_energy_release,
        total_dose_from_surrounding: 0.0,
    }
}

/// Apply cryo PE escape: track PE from a surrounding voxel toward crystal.
/// Returns dose deposited back into crystal.
///
/// Matches Java CrystalPolyhedron.addDoseAfterPECryo().
pub fn apply_cryo_pe_escape(
    cryo: &CryoEscape,
    crystal_size: [usize; 3],
    is_crystal_at: &dyn Fn(usize, usize, usize) -> bool,
    add_dose: &mut dyn FnMut(usize, usize, usize, f64),
    ci: f64,
    cj: f64,
    ck: f64,
    energy_pe: f64,
    energy_to_dose_factor: f64,
) -> f64 {
    let mut rng = rand::thread_rng();
    let mut dose_back = 0.0;
    let pe_dist_bins = cryo.cryo_pe_distances.len();

    if cryo.cryo_track_bias.is_empty() {
        return dose_back;
    }

    for _ in 0..(PE_ANGLE_RESOLUTION * PE_ANGLE_RESOLUTION) {
        let random_idx = rng.gen_range(0..cryo.cryo_track_bias.len());
        let random_track = cryo.cryo_track_bias[random_idx];

        for m in 0..pe_dist_bins {
            if random_track >= cryo.relative_vox_xyz_cryo[m].len() {
                continue;
            }
            let [rx, ry, rz] = cryo.relative_vox_xyz_cryo[m][random_track];

            let partial_energy = energy_pe * cryo.propn_dose_at_dist_cryo[m]
                / (PE_ANGLE_RESOLUTION * PE_ANGLE_RESOLUTION) as f64;
            let partial_dose = (partial_energy / energy_to_dose_factor) * 1e-6; // Energy to Dose in MGy

            let ti = (ci + rx).round() as i64;
            let tj = (cj + ry).round() as i64;
            let tk = (ck + rz).round() as i64;

            if ti >= 0
                && tj >= 0
                && tk >= 0
                && (ti as usize) < crystal_size[0]
                && (tj as usize) < crystal_size[1]
                && (tk as usize) < crystal_size[2]
                && is_crystal_at(ti as usize, tj as usize, tk as usize)
            {
                add_dose(ti as usize, tj as usize, tk as usize, partial_dose);
                dose_back += partial_dose;
            }
        }
    }

    dose_back
}

// ── Gumbel distribution helpers ───────────────────────────────────────────────

/// Gumbel location parameter coefficients based on density and PE energy.
/// Returns [a, b] where mu = a*E^2 + b*E.
/// Matches Java CrystalPolyhedron.setGumbelLoc().
fn set_gumbel_loc(density: f64, pe_energy: f64) -> [f64; 2] {
    let d2 = density * density;
    if pe_energy <= 20.0 {
        [
            0.0105 * d2 - 0.0351 * density + 0.0387,
            0.0245 * d2 - 0.0943 * density + 0.1171,
        ]
    } else if pe_energy <= 50.0 {
        [
            0.009 * d2 - 0.0293 * density + 0.0318,
            0.0459 * d2 - 0.1942 * density + 0.2582,
        ]
    } else {
        [
            0.0072 * d2 - 0.024 * density + 0.0262,
            0.16925 * d2 - 0.05562 * density + 0.6075,
        ]
    }
}

/// Gumbel scale parameter coefficients based on density and PE energy.
/// Returns [a, b] where beta = a*E^2 + b*E.
/// Matches Java CrystalPolyhedron.setGumbelScale().
fn set_gumbel_scale(density: f64, pe_energy: f64) -> [f64; 2] {
    let d2 = density * density;
    if pe_energy <= 20.0 {
        [
            0.0029 * d2 - 0.0081 * density + 0.0076,
            -0.0085 * density + 0.01751,
        ]
    } else if pe_energy <= 50.0 {
        [
            0.0018 * d2 - 0.0053 * density + 0.0051,
            -0.0238 * density + 0.0536,
        ]
    } else {
        [-0.0006 * density + 0.0017, -0.0569 * density + 0.1068]
    }
}

/// Calculate Gumbel parameters for a given beam energy.
/// Returns (mu, beta).
fn get_gumbel_params(pe_energy: f64, loc_coeffs: &[f64; 2], scale_coeffs: &[f64; 2]) -> (f64, f64) {
    let mu = loc_coeffs[0] * pe_energy * pe_energy + loc_coeffs[1] * pe_energy;
    let beta = scale_coeffs[0] * pe_energy * pe_energy + scale_coeffs[1] * pe_energy;
    (mu, beta)
}

/// Gumbel-Left PDF: f(x) = (1/beta) * exp(z - exp(z)) where z = (x-mu)/beta.
/// Note: when beta < 0 (can happen at low PE energies with high density),
/// this produces negative values. Java does not guard against this —
/// the negative signs cancel in the normalization within
/// calc_propn_dose_deposited, yielding correct positive proportions.
fn gumbel_pdf(x: f64, mu: f64, beta: f64) -> f64 {
    if beta == 0.0 {
        return 0.0;
    }
    let z = (x - mu) / beta;
    (z - z.exp()).exp() / beta
}

/// Calculate max PE distance (µm) — distance where Gumbel PDF drops to 0.1% of peak.
/// Matches Java CrystalPolyhedron.getMaxPEDistance().
fn get_max_pe_distance(pe_energy: f64, loc: &[f64; 2], scale: &[f64; 2]) -> f64 {
    let (mu, beta) = get_gumbel_params(pe_energy, loc, scale);
    let modal_height = gumbel_pdf(mu, mu, beta);
    let target_height = modal_height * 0.001;
    let step = if pe_energy < 10.0 {
        pe_energy / 1000.0
    } else {
        pe_energy / 100.0
    };

    let mut max_dist = mu;
    let mut calc_height = modal_height;
    while calc_height >= target_height {
        max_dist += step;
        calc_height = gumbel_pdf(max_dist, mu, beta);
        if max_dist >= mu * 2.0 {
            break;
        }
    }

    max_dist + 1.0
}

/// Calculate number of PE distance bins.
/// Matches Java CrystalPolyhedron.setMaxPEDistance() bin logic.
fn calc_pe_dist_bins(max_pe_dist: f64, pix_per_um: f64, pe_res: i32) -> usize {
    if pe_res >= 2 {
        return pe_res as usize;
    }
    let pixel_size = 1.0 / pix_per_um;
    if pixel_size >= max_pe_dist {
        2
    } else {
        let mut bins = 1 + (max_pe_dist / pixel_size).ceil() as usize;
        if bins < max_pe_dist as usize {
            bins = max_pe_dist as usize;
        }
        bins += 1; // erring on side of caution
        bins += 50; // ramping it up (matches Java)
        bins
    }
}

// ── Energy deposition distribution ────────────────────────────────────────────

/// Load energy coefficients from embedded CSV data.
/// Matches Java ReadEnergyCSV.openCSV().
fn load_energy_coefficients(csv_data: &str, pe_energy: f64) -> [f64; 7] {
    let used_energy = if pe_energy <= 1.0 {
        1.01
    } else if pe_energy >= 100.0 {
        99.99
    } else {
        pe_energy
    };

    let rounded = (used_energy * 100.0).round() / 100.0;
    let mut prev_line: Option<Vec<f64>> = None;

    for line in csv_data.lines() {
        let parts: Vec<f64> = line
            .split(',')
            .filter_map(|s| s.trim().parse().ok())
            .collect();
        if parts.len() < 8 {
            continue;
        }

        if parts[0] > rounded {
            // Choose nearest neighbor (interpolation <= 0.5 → use previous)
            let line_to_use = if let Some(ref prev) = prev_line {
                let interp = (used_energy - prev[0]) / (parts[0] - prev[0]);
                if interp <= 0.5 {
                    prev.clone()
                } else {
                    parts.clone()
                }
            } else {
                parts.clone()
            };

            // Reverse order: coefficients[0]=c0, [1]=c1, ..., [6]=c6
            let mut coeffs = [0.0; 7];
            for i in 0..7 {
                coeffs[i] = line_to_use[7 - i];
            }
            return coeffs;
        }
        prev_line = Some(parts);
    }

    [0.0; 7]
}

/// Calculate energy deposition distribution along a PE track.
/// Matches Java CrystalPolyhedron.calculateEnergyDistn().
fn calculate_energy_distn(distances: &[f64], bins: usize, pe_energy: f64) -> Vec<f64> {
    let low_coeffs = load_energy_coefficients(ENERGY_COEFS_LOW, pe_energy);
    let med_coeffs = load_energy_coefficients(ENERGY_COEFS_MED, pe_energy);
    let high_coeffs = load_energy_coefficients(ENERGY_COEFS_HIGH, pe_energy);

    let max_dist = if bins > 0 && bins - 1 < distances.len() {
        distances[bins - 1]
    } else {
        1.0
    };

    let mut energy = vec![0.0; bins];
    for i in 0..bins {
        let x = if max_dist > 0.0 {
            distances[i] / max_dist
        } else {
            0.0
        };

        if x == 0.0 {
            energy[i] = 0.0;
        } else if x < 2.0 / pe_energy {
            energy[i] = eval_poly6(&low_coeffs, x);
        } else if x > 1.0 - 2.0 / pe_energy {
            energy[i] = eval_poly6(&high_coeffs, x);
        } else {
            energy[i] = eval_poly6(&med_coeffs, x);
        }
    }
    energy
}

/// Evaluate 6th-degree polynomial: c6*x^6 + c5*x^5 + ... + c1*x + c0.
fn eval_poly6(coeffs: &[f64; 7], x: f64) -> f64 {
    coeffs[6] * x.powi(6)
        + coeffs[5] * x.powi(5)
        + coeffs[4] * x.powi(4)
        + coeffs[3] * x.powi(3)
        + coeffs[2] * x.powi(2)
        + coeffs[1] * x
        + coeffs[0]
}

/// Calculate proportion of dose deposited at each distance bin along PE tracks.
/// Combines Gumbel path-length distribution with energy deposition profile.
/// Matches Java CrystalPolyhedron.calcProportionVoxDoseDepositedByDist().
fn calc_propn_dose_deposited(
    distances: &[f64],
    pe_energy_raw: f64,
    loc: &[f64; 2],
    scale: &[f64; 2],
) -> Vec<f64> {
    let pe_dist_bins = distances.len();
    let (mu, beta) = get_gumbel_params(pe_energy_raw, loc, scale);

    // Calculate Gumbel path length distribution
    let mut path_length_distn = vec![0.0; pe_dist_bins];
    for i in 0..pe_dist_bins {
        path_length_distn[i] = gumbel_pdf(distances[i], mu, beta);
    }

    // Calculate energy distribution for each starting position
    let mut total_energy_distn = vec![vec![0.0; pe_dist_bins]; pe_dist_bins];
    for i in 0..pe_dist_bins {
        let sub_bins = pe_dist_bins - i;
        let energy_distn = calculate_energy_distn(distances, sub_bins, pe_energy_raw);
        for j in (1..energy_distn.len()).rev() {
            total_energy_distn[i][j] = energy_distn[j];
        }
    }

    // Integrate path length distribution
    let mut distn_integral = 0.0;
    let mut tot_energy_integral = vec![0.0; pe_dist_bins];
    for l in 0..pe_dist_bins.saturating_sub(1) {
        let width = distances[l + 1] - distances[l];
        let height = (path_length_distn[l + 1] + path_length_distn[l]) / 2.0;
        distn_integral += width * height;

        for i in 0..pe_dist_bins {
            if total_energy_distn[i][l + 1] != 0.0 {
                let energy_height = (total_energy_distn[i][l + 1] + total_energy_distn[i][l]) / 2.0;
                tot_energy_integral[i] += width * energy_height;
            }
        }
    }

    // Calculate distance widths and heights (reversed order)
    let mut distance_widths = vec![0.0; pe_dist_bins];
    let mut distance_heights = vec![0.0; pe_dist_bins];
    let mut path_count: usize = 0;
    for l in (1..pe_dist_bins).rev() {
        distance_widths[path_count] = distances[l] - distances[l - 1];
        distance_heights[path_count] = (path_length_distn[l] + path_length_distn[l - 1]) / 2.0;
        path_count += 1;
    }

    // Combine path length and energy distributions
    let mut propn = vec![0.0; pe_dist_bins];
    for l in (1..pe_dist_bins).rev() {
        for i in 0..pe_dist_bins {
            let mut energy_height = 0.0;
            if total_energy_distn[i][l] != 0.0 {
                energy_height = (total_energy_distn[i][l] + total_energy_distn[i][l - 1]) / 2.0;
            }

            if distn_integral != 0.0 && tot_energy_integral[i] != 0.0 {
                propn[l] += (distance_widths[i] * distance_heights[i] / distn_integral)
                    * (distance_widths[i] * energy_height / tot_energy_integral[i]);
            }
        }
    }

    propn
}

// ── Angular distribution ──────────────────────────────────────────────────────

/// Set up PE angular distribution biased by beam polarization.
/// Matches Java CrystalPolyhedron.setUpPEPolarisation().
fn setup_pe_polarisation(
    _coefcalc: &dyn CoefCalc,
    _beam_energy: f64,
    fe_factors: &[Vec<f64>],
) -> Vec<f64> {
    let beta = 2.0; // s-shell approximation

    // Calculate weight from element absorption
    // Java: weight += feFactors[el][0] * feFactors[el][2]
    // Java [2] = K_ion_prob, so weight = sum(mu_ratio * K_ion_prob)
    // Rust [1] = K_ion_prob
    let mut weight = 0.0;
    for row in fe_factors {
        weight += row[0] * safe_idx(row, 1);
    }

    let num_bins = PE_ANGLE_RES_LIMIT / 2 + 1;
    let mut weighted_avg = vec![0.0; num_bins];

    for i in 0..num_bins {
        let angle = i as f64 * (PE_ANGLE_LIMIT / PE_ANGLE_RES_LIMIT as f64);
        let point = solve_polarisation(angle, 1.0, beta) / solve_polarisation(0.0, 1.0, beta);
        let sum_point = point * weight * DEGREE_OF_POLARISATION;
        weighted_avg[i] =
            (1000.0 * sum_point).round() + 1000.0 * (1.0 - weight * DEGREE_OF_POLARISATION);
    }

    weighted_avg
}

/// Solve polarisation equation for angle.
/// Returns differential cross-section value.
fn solve_polarisation(phi: f64, photo_electric: f64, beta: f64) -> f64 {
    (photo_electric / (4.0 * PI)) * (1.0 + beta * 0.5 * (3.0 * phi.cos().powi(2) - 1.0))
}

// ── Voxel offset pre-calculation ──────────────────────────────────────────────

/// Pre-calculate PE track voxel offsets and angular bias array.
/// Matches Java CrystalPolyhedron.findVoxelsReachedByPE() (non-cryo case).
fn find_voxels_reached_by_pe(
    pe_distances: &[f64],
    pix_per_um: f64,
    angular_distribution: &[f64],
) -> (Vec<Vec<[f64; 3]>>, Vec<usize>) {
    let pe_dist_bins = pe_distances.len();
    let step = 2.0 * PI / PE_ANGLE_RES_LIMIT as f64;
    let max_tracks = PE_ANGLE_RES_LIMIT * PE_ANGLE_RES_LIMIT;

    let mut relative_vox_xyz = vec![vec![[0.0f64; 3]; max_tracks]; pe_dist_bins];
    let mut big_array = Vec::new();
    let mut counter: i32 = -1;

    let mut theta = 0.0;
    while theta < 2.0 * PI {
        let mut phi = 0.0;
        while phi <= PE_ANGLE_LIMIT / 2.0 {
            // Check for duplicate tracks at poles
            let replicate = if (theta == 0.0) || ((theta - PE_ANGLE_LIMIT / 2.0).abs() < 1e-12) {
                phi != 0.0
            } else {
                false
            };

            if !replicate {
                counter += 1;
                let c = counter as usize;

                let x_norm = theta.sin() * phi.cos();
                let y_norm = theta.sin() * phi.sin();
                let z_norm = theta.cos();

                // Calculate angle to beam axis (horizontal goniometer: y-axis)
                let dot_product = y_norm; // horizontal goniometer
                let magnitude = (x_norm * x_norm + y_norm * y_norm + z_norm * z_norm).sqrt();
                let cos_angle = dot_product / magnitude;
                let angle_to_x = cos_angle.acos();

                // Find angular bin
                let place =
                    (angle_to_x * PE_ANGLE_RES_LIMIT as f64 / PE_ANGLE_LIMIT).round() as usize;
                let place = place.min(angular_distribution.len() - 1);

                // Weight track in bias array
                let weight = angular_distribution[place].max(0.0) as usize;
                for _ in 0..weight {
                    big_array.push(c);
                }

                // Pre-calculate voxel offsets for each distance bin
                // No rotation applied (angle=0 for initial setup)
                for m in 0..pe_dist_bins {
                    let r = pe_distances[m] * pix_per_um;
                    if c < max_tracks {
                        relative_vox_xyz[m][c] = [r * x_norm, r * y_norm, r * z_norm];
                    }
                }
            }

            phi += step;
        }
        theta += step;
    }

    // Trim big_array (skip trailing zeros like Java)
    let track_bias = if big_array.is_empty() {
        vec![0]
    } else {
        big_array
    };

    (relative_vox_xyz, track_bias)
}

/// Pre-calculate FL track voxel offsets.
/// Matches Java CrystalPolyhedron.findVoxelsReachedByFL().
fn find_voxels_reached_by_fl(
    fe_factors: &[Vec<f64>],
    fl_distances: &Vec<Vec<Vec<f64>>>,
    fl_dist_bins: usize,
    pix_per_um: f64,
) -> (Vec<Vec<Vec<Vec<[f64; 3]>>>>, usize) {
    let step = PE_ANGLE_LIMIT / FL_ANGLE_RES_LIMIT as f64;
    let max_angle_bins = FL_ANGLE_RES_LIMIT * FL_ANGLE_RES_LIMIT;
    let num_elements = fe_factors.len();

    let mut fl_rel =
        vec![vec![vec![vec![[0.0f64; 3]; max_angle_bins]; fl_dist_bins]; 4]; num_elements];
    let mut counter: i32 = -1;

    let mut phi = 0.0;
    while phi < PE_ANGLE_LIMIT {
        let mut theta = 0.0;
        while theta <= PE_ANGLE_LIMIT / 2.0 {
            let replicate = if (theta == 0.0) || ((theta - PE_ANGLE_LIMIT / 2.0).abs() < 1e-12) {
                phi != 0.0
            } else {
                false
            };

            if !replicate {
                counter += 1;
                let c = counter as usize;

                let x_norm = theta.sin() * phi.cos();
                let y_norm = theta.sin() * phi.sin();
                let z_norm = theta.cos();

                for n in 0..num_elements {
                    for l in 0..4 {
                        for m in 0..fl_dist_bins {
                            if n < fl_distances.len()
                                && l < fl_distances[n].len()
                                && m < fl_distances[n][l].len()
                            {
                                let r = fl_distances[n][l][m] * pix_per_um;
                                if c < max_angle_bins {
                                    fl_rel[n][l][m][c] = [r * x_norm, r * y_norm, r * z_norm];
                                }
                            }
                        }
                    }
                }
            }

            theta += step;
        }
        phi += step;
    }

    let num_tracks = (counter + 1).max(0) as usize;
    (fl_rel, num_tracks)
}

// ── Fluorescence distribution ─────────────────────────────────────────────────

/// Calculate FL distance distribution and distances for each element/shell.
/// Matches Java CrystalPolyhedron.calcFluorescenceDistribution().
///
/// Returns (distribution[element][shell][dist], distances[element][shell][dist]).
fn calc_fluorescence_distribution(
    fe_factors: &[Vec<f64>],
    fl_dist_bins: usize,
    cryst_size_um: [f64; 3],
) -> (Vec<Vec<Vec<f64>>>, Vec<Vec<Vec<f64>>>) {
    let num_elements = fe_factors.len();
    let distance_resolution = fl_dist_bins - 1;

    let mut distribution = vec![vec![vec![0.0; fl_dist_bins]; 4]; num_elements];
    let mut distances = vec![vec![vec![0.0; fl_dist_bins]; 4]; num_elements];

    // Crystal diagonal distance
    let crystal_max_dist =
        (cryst_size_um[0].powi(2) + cryst_size_um[1].powi(2) + cryst_size_um[2].powi(2)).sqrt();

    for i in 0..num_elements {
        let row = &fe_factors[i];
        for j in 0..4 {
            // Check if this shell has fluorescence
            // FL proportion is computed from fl_energy_per_event / total, but here we just
            // check if there's a valid muabs coefficient
            // Java bug compatibility: muabsIndex is hardcoded to 4 (K-shell) for all
            // shells instead of the correct (4*j)+4. The correct formula is commented
            // out in CrystalPolyhedron.java:1239. See docs/java-bugs-analysis.md §Bug 1.
            let muabs_index = 4; // correct: (4 * j) + 4
            let muabs = safe_idx(row, muabs_index);
            // Check if shell is active: ion_prob > 0 and fl_yield > 0
            let ion_idx = match j {
                0 => 1,
                1 => 5,
                2 => 9,
                3 => 13,
                _ => continue,
            };
            let fl_idx = ion_idx + 1;
            let ion_prob = safe_idx(row, ion_idx);
            let fl_yield = safe_idx(row, fl_idx);

            if ion_prob <= 0.0 || fl_yield <= 0.0 {
                continue;
            }

            // Max distance where escape probability = 5%.
            // When muabs <= 0 (e.g. K-shell escapeMuAbs is 0 because beam is below
            // K-edge — see Bug 1 in docs/java-bugs-analysis.md), Java computes
            // -ln(0.05)/0 = Infinity, then clamps to crystal diagonal. The resulting
            // exp(-0 * d) = 1.0 distribution puts 100% in the escape bin.
            let mut max_dist_fl = if muabs <= 0.0 || !muabs.is_finite() {
                f64::INFINITY
            } else {
                -(0.05_f64.ln()) / muabs
            };
            if max_dist_fl > crystal_max_dist {
                max_dist_fl = crystal_max_dist;
            }

            // Populate distances
            for q in 0..=distance_resolution {
                distances[i][j][q] = (max_dist_fl / distance_resolution as f64) * q as f64;
            }

            // Calculate escape probabilities
            let mut running_escape = 0.0;
            for l in (0..fl_dist_bins).rev() {
                if l == fl_dist_bins - 1 {
                    // Escape probability at max distance
                    distribution[i][j][l] = (-muabs * distances[i][j][l]).exp();
                } else {
                    distribution[i][j][l] = (-muabs * distances[i][j][l]).exp() - running_escape;
                }
                running_escape += distribution[i][j][l];
            }
        }
    }

    (distribution, distances)
}

// ── Utility ───────────────────────────────────────────────────────────────────

/// Safe index into a Vec<f64>, returning 0.0 for out-of-bounds.
#[inline]
fn safe_idx(v: &[f64], idx: usize) -> f64 {
    v.get(idx).copied().unwrap_or(0.0)
}
