pub mod cuboid;
#[allow(
    clippy::too_many_arguments,
    clippy::type_complexity,
    clippy::needless_range_loop,
    clippy::explicit_counter_loop,
    clippy::ptr_arg
)]
pub mod escape;
pub mod polyhedron;
pub mod spherical;

pub use cuboid::CrystalCuboid;
pub use polyhedron::CrystalPolyhedron;
pub use spherical::CrystalSpherical;

use crate::beam::Beam;
use crate::coefcalc::CoefCalc;
use crate::constants::*;
use crate::container::Container;
use crate::ddm::DdmModel;
use crate::parser::config::CrystalType;
use crate::wedge::Wedge;

/// Crystal trait: all crystal geometries implement this.
pub trait Crystal: std::fmt::Debug + Send + Sync {
    /// Set up depth finding for a given angle and wedge.
    fn setup_depth_finding(&mut self, ang_rad: f64, wedge: &Wedge);

    /// Find depth of a voxel coordinate from the crystal surface along [0,0,1].
    fn find_depth(&self, vox_coord: &[f64; 3], delta_phi: f64, wedge: &Wedge) -> f64;

    /// Return physical coordinates (µm) of voxel [i,j,k].
    fn get_cryst_coord(&self, i: usize, j: usize, k: usize) -> [f64; 3];

    /// Return number of images in a wedge.
    fn num_images(&self, wedge: &Wedge) -> usize;

    /// Is there crystal material at voxel [i,j,k]?
    fn is_crystal_at(&self, i: usize, j: usize, k: usize) -> bool;

    /// Add dose to voxel.
    fn add_dose(&mut self, i: usize, j: usize, k: usize, dose_increase: f64);

    /// Add fluence to voxel.
    fn add_fluence(&mut self, i: usize, j: usize, k: usize, fluence_increase: f64);

    /// Add elastic scattering to voxel.
    fn add_elastic(&mut self, i: usize, j: usize, k: usize, elastic_increase: f64);

    /// Return crystal description.
    fn crystal_info(&self) -> String;

    /// Return crystal size in voxels [nx, ny, nz].
    fn cryst_size_voxels(&self) -> [usize; 3];

    /// Return crystal size in µm [x, y, z].
    fn cryst_size_um(&self) -> [f64; 3];

    /// Return dose at voxel in MGy.
    fn get_dose(&self, i: usize, j: usize, k: usize) -> f64;

    /// Return fluence at voxel.
    fn get_fluence(&self, i: usize, j: usize, k: usize) -> f64;

    /// Return elastic scattering at voxel.
    fn get_elastic(&self, i: usize, j: usize, k: usize) -> f64;

    /// Crystal resolution in voxels per µm.
    fn crystal_pix_per_um(&self) -> f64;

    /// Mutable access to the CoefCalc.
    fn coefcalc(&self) -> &dyn CoefCalc;
    fn coefcalc_mut(&mut self) -> &mut dyn CoefCalc;

    /// Access to the DDM.
    fn ddm(&self) -> &dyn DdmModel;

    /// Access to the ExposureSummary observer.
    fn exposure_summary(&self) -> &crate::output::ExposureSummary;
    fn exposure_summary_mut(&mut self) -> &mut crate::output::ExposureSummary;

    /// The subprogram mode (RD3D, XFEL, etc.)
    fn subprogram(&self) -> &str;

    /// Whether photoelectron escape is enabled.
    fn photo_electron_escape(&self) -> bool {
        false
    }

    /// Whether fluorescent escape is enabled.
    fn fluorescent_escape(&self) -> bool {
        false
    }

    /// Whether surrounding (cryo) calculation is enabled.
    fn calc_surrounding(&self) -> bool {
        false
    }

    /// PE resolution override (None = use default).
    fn pe_resolution(&self) -> Option<i32> {
        None
    }

    /// Expose the crystal to a beam according to a wedge strategy.
    fn expose(&mut self, beam: &mut dyn Beam, wedge: &Wedge);
}

/// Translate crystal coordinates to beam position.
/// Applies wedge start offset, translation, and rotation by angle.
pub fn translate_crystal_to_position(
    cryst_coords: &[f64; 3],
    wedge_start: &[f64; 3],
    wedge_translation: &[f64; 3],
    angle_cos: f64,
    angle_sin: f64,
) -> [f64; 3] {
    // Apply translations
    let y = cryst_coords[1] + wedge_start[1] + wedge_translation[1];
    let x = cryst_coords[0] + wedge_start[0] + wedge_translation[0];
    let z = cryst_coords[2] + wedge_start[2] + wedge_translation[2];

    // Rotate X-Z plane by angle (rotation about Y axis)
    [
        x * angle_cos + z * angle_sin,
        y,
        -x * angle_sin + z * angle_cos,
    ]
}

/// Default expose implementation for standard RD3D mode.
/// This is the main simulation loop, extracted as a free function
/// so crystal implementations can call it.
pub fn expose_rd3d(
    crystal: &mut dyn Crystal,
    beam: &mut dyn Beam,
    wedge: &Wedge,
    container: &mut dyn Container,
) {
    use std::io::Write;
    // Update coefficients for beam energy
    crystal
        .coefcalc_mut()
        .update_coefficients(beam.photon_energy());

    // Apply container attenuation
    container.calculate_attenuation(beam.photon_energy());
    println!("{}", container.info());
    beam.apply_container_attenuation(container);

    // Generate beam array (only relevant for experimental beam)
    beam.generate_beam_array();

    // Set up escape parameters if enabled
    let pe_enabled = crystal.photo_electron_escape();
    let fl_enabled = crystal.fluorescent_escape();
    let cryo_enabled = crystal.calc_surrounding() && crystal.coefcalc().is_cryo();

    let fe_factors = crystal
        .coefcalc()
        .fluorescent_escape_factors(beam.photon_energy());

    let pe_escape = if pe_enabled {
        Some(escape::setup_pe_escape(
            beam.photon_energy(),
            crystal.coefcalc(),
            &fe_factors,
            crystal.crystal_pix_per_um(),
            crystal.cryst_size_um(),
            crystal.pe_resolution(),
        ))
    } else {
        None
    };

    let fl_escape = if fl_enabled {
        Some(escape::setup_fl_escape(
            &fe_factors,
            crystal.crystal_pix_per_um(),
            crystal.cryst_size_um(),
        ))
    } else {
        None
    };

    let cryo_fe_factors = if cryo_enabled {
        crystal
            .coefcalc_mut()
            .update_cryo_coefficients(beam.photon_energy());
        crystal
            .coefcalc()
            .cryo_fluorescent_escape_factors(beam.photon_energy())
    } else {
        vec![]
    };

    let cryo_escape = if cryo_enabled && pe_enabled {
        let crystal_pe_dist_bins = pe_escape.as_ref().map(|p| p.pe_dist_bins).unwrap_or(2);
        Some(escape::setup_cryo_escape(
            beam.photon_energy(),
            crystal.coefcalc(),
            &cryo_fe_factors,
            crystal.crystal_pix_per_um(),
            crystal.cryst_size_um(),
            fl_enabled,
            beam.beam_minimum_dimension(),
            crystal_pe_dist_bins,
        ))
    } else {
        None
    };

    // Set up angles
    let angles = compute_angles(wedge);
    let crystal_size = crystal.cryst_size_voxels();

    // Notify exposure summary
    crystal
        .exposure_summary_mut()
        .exposure_start(angles.len(), wedge, crystal_size);

    // Progress indicator (matches Java OutputProgressIndicator)
    let image_count = angles.len();
    let mut wedge_progress: usize = 0;
    print!("Exposing wedge: [ 0%");
    let _ = std::io::stdout().flush();

    // Energy sampling (monochromatic for now)
    let energies = vec![beam.photon_energy()];
    let energies_per_angle = energies.len();

    // Mutable escape tracking
    let mut total_escaped_pe = 0.0_f64;
    let mut total_escaped_fl = 0.0_f64;
    let mut _total_dose_from_surrounding = 0.0_f64;

    // Main loop: angles × energies × voxels
    for (n, &angle) in angles.iter().enumerate() {
        for &photon_energy in &energies {
            // Update coefficients for this energy
            crystal.coefcalc_mut().update_coefficients(photon_energy);

            let (esc_pe, esc_fl) = expose_angle(
                crystal,
                beam,
                wedge,
                angle,
                n,
                angles.len(),
                energies_per_angle,
                photon_energy,
                pe_escape.as_ref(),
                fl_escape.as_ref(),
                cryo_escape.as_ref().map(|c| c.cryo_pix_per_um),
            );
            total_escaped_pe += esc_pe;
            total_escaped_fl += esc_fl;
        }

        // Cryo surrounding loop (after crystal exposure for this angle)
        if let (Some(ref pe), Some(ref cryo)) = (&pe_escape, &cryo_escape) {
            let dose_back = expose_cryo_angle(
                crystal,
                beam,
                wedge,
                angle,
                n,
                angles.len(),
                energies_per_angle,
                beam.photon_energy(),
                cryo,
                pe,
            );
            _total_dose_from_surrounding += dose_back;
        }

        // Notify image complete
        let last_angle = if n > 0 { angles[n - 1] } else { 0.0 };
        let vox_vol = crystal.crystal_pix_per_um().powi(-3);
        crystal
            .exposure_summary_mut()
            .image_complete(n, angle, last_angle, vox_vol);

        // Update progress bar
        if image_count > 0 {
            while 100 * (n + 1) / image_count > wedge_progress {
                wedge_progress += 1;
                if wedge_progress.is_multiple_of(4) {
                    print!(".");
                }
                match wedge_progress {
                    20 => print!("20%"),
                    40 => print!("40%"),
                    60 => print!("60%"),
                    80 => print!("80%"),
                    100 => print!("100%"),
                    _ => {}
                }
            }
            let _ = std::io::stdout().flush();
        }
    }

    // Close progress bar
    println!(" ]");

    // Print escape summaries (matches Java)
    if pe_enabled {
        println!(
            "\nEnergy that may escape by Photoelectron Escape: {:.2e} J.\n",
            total_escaped_pe
        );
    }
    if fl_enabled {
        println!(
            "Total energy that may escape by Fluorescent Escape: {:.2e} J.\n",
            total_escaped_fl
        );
    }

    // Summary observations
    let voxel_mass_kg =
        UNIT_CONVERSION * (crystal.crystal_pix_per_um().powi(-3) * crystal.coefcalc().density());

    let size = crystal.cryst_size_voxels();
    for i in 0..size[0] {
        for j in 0..size[1] {
            for k in 0..size[2] {
                if crystal.is_crystal_at(i, j, k) {
                    let dose = crystal.get_dose(i, j, k);
                    crystal.exposure_summary_mut().summary_observation(
                        i,
                        j,
                        k,
                        dose,
                        voxel_mass_kg,
                    );
                }
            }
        }
    }

    crystal.exposure_summary_mut().exposure_complete();
}

/// Compute the angle array for a wedge.
fn compute_angles(wedge: &Wedge) -> Vec<f64> {
    // Normalize angles to [0, 2π)
    let twopi = 2.0 * std::f64::consts::PI;
    let diff = (wedge.start_ang / twopi).floor() as i64;
    let start = wedge.start_ang - twopi * diff as f64;
    let end = wedge.end_ang - twopi * diff as f64;

    if (start - end).abs() < wedge.ang_res {
        // Static exposure: no rotation
        vec![start; STATIC_EXPOSURE]
    } else {
        let sign: f64 = if end < start { -1.0 } else { 1.0 };
        let count = (sign * (end - start) / wedge.ang_res + 1.0) as usize;
        (0..count)
            .map(|i| start + sign * i as f64 * wedge.ang_res)
            .collect()
    }
}

/// Expose one angle of the crystal. Returns (escaped_pe, escaped_fl).
/// `cryo_ppm`: if cryo surrounding is active, the cryo grid pixels-per-micron
/// (may differ from crystal PPM). Used for absorbed energy bug-compat only.
#[allow(clippy::too_many_arguments)]
fn expose_angle(
    crystal: &mut dyn Crystal,
    beam: &dyn Beam,
    wedge: &Wedge,
    angle: f64,
    angle_num: usize,
    angle_count: usize,
    energies_per_angle: usize,
    photon_energy: f64,
    pe_escape: Option<&escape::PeEscape>,
    fl_escape: Option<&escape::FlEscape>,
    cryo_ppm: Option<f64>,
) -> (f64, f64) {
    let crystal_size = crystal.cryst_size_voxels();

    let wedge_start = wedge.start_vector();
    let wedge_translation = wedge.translation_vector(angle);

    let angle_cos = angle.cos();
    let angle_sin = angle.sin();

    crystal.setup_depth_finding(angle, wedge);

    let pix_per_um = crystal.crystal_pix_per_um();
    let cc = crystal.coefcalc();
    let abs_coeff = cc.absorption_coefficient();
    let inelastic_coeff = cc.inelastic_coefficient();
    let elastic_coeff = cc.elastic_coefficient();
    let att_coeff = cc.attenuation_coefficient();
    let density = cc.density();
    // Java bug compat: capture cryo coefficient and PPM upfront
    // (see docs/java-bugs-analysis.md §Bug 2)
    let cryo_abs_coeff_if_active = if cc.is_cryo() {
        Some(cc.cryo_absorption_coefficient())
    } else {
        None
    };

    // Fluence to dose factor (photoelectric)
    let fluence_to_dose_factor = -(-abs_coeff / pix_per_um).exp_m1()
        / (UNIT_CONVERSION * (pix_per_um.powi(-3) * density))
        * GY_TO_MGY;

    // Fluence to dose factor (Compton)
    let fluence_to_dose_factor_compton = -(-inelastic_coeff / pix_per_um).exp_m1()
        / (UNIT_CONVERSION * (pix_per_um.powi(-3) * density))
        * GY_TO_MGY;

    // Fluence to elastic factor
    let fluence_to_elastic_factor =
        -(-elastic_coeff / pix_per_um).exp_m1() / (photon_energy * KEV_TO_JOULES);

    // Beam attenuation factor
    let beam_attenuation_factor =
        pix_per_um.powi(-2) * wedge.total_sec() / (angle_count * energies_per_angle) as f64;

    let beam_attenuation_exp_factor = -att_coeff;

    let beam_energy_j = photon_energy * KEV_TO_JOULES;

    // Compton electron energy calculation constants
    let electron_mass_kg: f64 = 9.109_383_56e-31;
    let c_squared: f64 = 3.0e8 * 3.0e8;
    let mc_squared = electron_mass_kg * c_squared;

    let mut total_escaped_pe = 0.0_f64;
    let mut total_escaped_fl = 0.0_f64;

    // Build a snapshot of crystal occupancy for escape closures.
    // We need to pass is_crystal_at and add_dose as closures but both
    // borrow crystal. Instead, we collect voxel data first, then apply.
    // Actually, we can use the fact that the escape functions need
    // shared access while we have exclusive access.
    // We'll use a two-pass approach for escape: collect doses first, apply later.

    for i in 0..crystal_size[0] {
        for j in 0..crystal_size[1] {
            for k in 0..crystal_size[2] {
                if !crystal.is_crystal_at(i, j, k) {
                    continue;
                }

                let cryst_coords = crystal.get_cryst_coord(i, j, k);
                let translated = translate_crystal_to_position(
                    &cryst_coords,
                    &wedge_start,
                    &wedge_translation,
                    angle_cos,
                    angle_sin,
                );

                let unattenuated_intensity =
                    beam.beam_intensity(translated[0], translated[1], wedge.off_axis_um);

                if unattenuated_intensity <= 0.0 {
                    continue;
                }

                let depth = crystal.find_depth(&translated, angle, wedge);

                // Attenuated fluence (J)
                let vox_fluence = unattenuated_intensity
                    * beam_attenuation_factor
                    * (depth * beam_attenuation_exp_factor).exp();

                // Compton electron energy
                let compton_ratio = mc_squared / (2.0 * beam_energy_j + mc_squared);
                let compton_electron_energy = beam_energy_j * (1.0 - compton_ratio.sqrt());

                let num_photons = vox_fluence / beam_energy_j;
                let compton_fluence = num_photons * compton_electron_energy;

                // Elastic yield
                let elastic_yield = fluence_to_elastic_factor * vox_fluence;

                // Dose from photoelectric absorption
                let vox_dose = fluence_to_dose_factor * vox_fluence;

                // Dose from Compton
                let vox_dose_compton = fluence_to_dose_factor_compton * compton_fluence;

                if vox_dose > 0.0 {
                    let total_dose_before = crystal.get_dose(i, j, k);
                    // Match Java: interpolatedVoxelDose = getDose(after) + voxImageDose/2
                    // Java's DWD loop runs after addDose, so getDose(after) = total_before + PE + Compton.
                    // interpolatedVoxelDose = (total_before + PE + Compton) + PE/2
                    //                      = total_before + 1.5*PE + Compton
                    let rde = crystal
                        .ddm()
                        .calc_decay(total_dose_before + 1.5 * vox_dose + vox_dose_compton);
                    // Java bug compatibility: when cryo surrounding + PE escape are both
                    // active, Crystal.java:1236 overwrites energyPerFluence with the cryo
                    // absorption coefficient AND cryo PPM, and that leaked value is used at
                    // line 1361 for crystal absorbed energy.
                    // See docs/java-bugs-analysis.md §Bug 2.
                    let energy_per_fluence = if let (Some(cryo_coeff), Some(_), Some(cppm)) =
                        (cryo_abs_coeff_if_active, &pe_escape, cryo_ppm)
                    {
                        // Java scope-leak: uses cryo abs coeff AND cryo PPM
                        -(-cryo_coeff / cppm).exp_m1()
                    } else {
                        // correct: crystal absorption coefficient and crystal PPM
                        -(-abs_coeff / pix_per_um).exp_m1()
                    };
                    let energy_absorbed =
                        energy_per_fluence * vox_fluence + energy_per_fluence * compton_fluence;

                    // Notify exposure summary.
                    // added_dose includes both PE and Compton, matching Java's
                    // addedDose = totalVoxelDose - voxImageDoseLast (which includes Compton).
                    crystal.exposure_summary_mut().exposure_observation(
                        angle_num,
                        i,
                        j,
                        k,
                        vox_dose + vox_dose_compton,
                        total_dose_before,
                        vox_fluence,
                        rde,
                        energy_absorbed,
                        elastic_yield,
                        angle_count as f64,
                    );

                    crystal.add_fluence(i, j, k, vox_fluence);

                    // Apply escape logic based on enabled flags
                    match (pe_escape, fl_escape) {
                        (Some(pe), Some(fl)) => {
                            // Case 1: Both PE and FL escape
                            let tot_fl_energy = pe.fl_energy_release * num_photons;
                            let vox_fl_dose = fluence_to_dose_factor * tot_fl_energy;
                            let tot_auger_dose =
                                pe.auger_energy * num_photons * fluence_to_dose_factor;

                            // FL escape
                            if vox_fl_dose > 0.0 {
                                let fl_lost = apply_escape_fl(
                                    crystal,
                                    pe,
                                    fl,
                                    crystal_size,
                                    i,
                                    j,
                                    k,
                                    vox_fl_dose,
                                );
                                total_escaped_fl += fl_lost;
                            }

                            // Auger stays in voxel
                            if tot_auger_dose > 0.0 {
                                crystal.add_dose(i, j, k, tot_auger_dose);
                            }

                            // PE dose = total - FL - Auger
                            let dose_pe = vox_dose - vox_fl_dose - tot_auger_dose;
                            if dose_pe > 0.0 {
                                let pe_lost =
                                    apply_escape_pe(crystal, pe, crystal_size, i, j, k, dose_pe);
                                total_escaped_pe += pe_lost;
                            }
                        }
                        (Some(pe), None) => {
                            // Case 2: PE escape only
                            let tot_auger_dose =
                                pe.auger_energy * num_photons * fluence_to_dose_factor;

                            // PE dose = total - binding energy fraction
                            let binding_frac = if photon_energy > 0.0 {
                                pe.energy_to_subtract / photon_energy
                            } else {
                                0.0
                            };
                            let dose_pe = vox_dose - binding_frac * vox_dose;

                            if dose_pe > 0.0 {
                                let pe_lost =
                                    apply_escape_pe(crystal, pe, crystal_size, i, j, k, dose_pe);
                                total_escaped_pe += pe_lost;
                            }

                            // Auger stays in voxel
                            if tot_auger_dose > 0.0 {
                                crystal.add_dose(i, j, k, tot_auger_dose);
                            }
                        }
                        (None, Some(_fl)) => {
                            // Case 3: FL escape only (need pe params for fl_proportion)
                            // Without PE setup, FL escape requires fl_energy_release data
                            // which comes from PE setup. Fall through to no-escape.
                            crystal.add_dose(i, j, k, vox_dose);
                        }
                        (None, None) => {
                            // Case 4: No escape
                            crystal.add_dose(i, j, k, vox_dose);
                        }
                    }

                    // Compton dose and elastic always added
                    crystal.add_dose(i, j, k, vox_dose_compton);
                    crystal.add_elastic(i, j, k, elastic_yield);
                } else if vox_dose < 0.0 {
                    panic!("Negative dose encountered - this should never happen");
                }
            }
        }
    }

    (total_escaped_pe, total_escaped_fl)
}

/// Helper: apply PE escape for a single voxel. Works around borrow issues
/// by doing escape calculations directly.
fn apply_escape_pe(
    crystal: &mut dyn Crystal,
    pe: &escape::PeEscape,
    crystal_size: [usize; 3],
    i: usize,
    j: usize,
    k: usize,
    dose_pe: f64,
) -> f64 {
    use rand::Rng;

    let mut rng = rand::thread_rng();
    let mut dose_lost = 0.0;

    if pe.track_bias.is_empty() {
        return dose_lost;
    }

    for _ in 0..(escape::pe_angle_resolution() * escape::pe_angle_resolution()) {
        let random_idx = rng.gen_range(0..pe.track_bias.len());
        let random_track = pe.track_bias[random_idx];

        for m in 0..pe.pe_dist_bins {
            if random_track >= pe.relative_vox_xyz[m].len() {
                continue;
            }
            let [rx, ry, rz] = pe.relative_vox_xyz[m][random_track];
            let partial_dose = dose_pe * pe.propn_dose_at_dist[m]
                / (escape::pe_angle_resolution() * escape::pe_angle_resolution()) as f64;

            let ti = (i as f64 + rx).round() as i64;
            let tj = (j as f64 + ry).round() as i64;
            let tk = (k as f64 + rz).round() as i64;

            if ti >= 0
                && tj >= 0
                && tk >= 0
                && (ti as usize) < crystal_size[0]
                && (tj as usize) < crystal_size[1]
                && (tk as usize) < crystal_size[2]
                && crystal.is_crystal_at(ti as usize, tj as usize, tk as usize)
            {
                crystal.add_dose(ti as usize, tj as usize, tk as usize, partial_dose);
            } else {
                dose_lost += partial_dose;
            }
        }
    }

    dose_lost
}

/// Helper: apply FL escape for a single voxel.
#[allow(clippy::too_many_arguments)]
fn apply_escape_fl(
    crystal: &mut dyn Crystal,
    pe: &escape::PeEscape,
    fl: &escape::FlEscape,
    crystal_size: [usize; 3],
    i: usize,
    j: usize,
    k: usize,
    dose_fl: f64,
) -> f64 {
    let mut dose_lost = 0.0;
    let num_elements = pe.fl_proportion_event.len();

    for n in 0..num_elements {
        for l in 0..4 {
            for m in 0..fl.fl_dist_bins {
                for q in 0..(escape::FL_ANGLE_RESOLUTION_PUB * escape::FL_ANGLE_RESOLUTION_PUB) {
                    if pe.fl_proportion_event[n][l] == 0.0 {
                        continue;
                    }
                    if n >= fl.fl_relative_vox_xyz.len()
                        || l >= fl.fl_relative_vox_xyz[n].len()
                        || m >= fl.fl_relative_vox_xyz[n][l].len()
                        || q >= fl.fl_relative_vox_xyz[n][l][m].len()
                    {
                        continue;
                    }

                    let fl_partial_dose =
                        dose_fl * pe.fl_proportion_event[n][l] * fl.fl_dist_distribution[n][l][m]
                            / (escape::FL_ANGLE_RESOLUTION_PUB * escape::FL_ANGLE_RESOLUTION_PUB)
                                as f64;

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
                        && crystal.is_crystal_at(ti as usize, tj as usize, tk as usize)
                    {
                        crystal.add_dose(ti as usize, tj as usize, tk as usize, fl_partial_dose);
                    } else {
                        dose_lost += fl_partial_dose;
                    }
                }
            }
        }
    }

    dose_lost
}

/// Expose surrounding (cryo) voxels for one angle and calculate PE dose
/// that enters the crystal from surrounding material.
/// Returns total dose deposited back into crystal from surrounding.
/// Matches Java Crystal.java lines 1220-1340.
#[allow(clippy::too_many_arguments)]
fn expose_cryo_angle(
    crystal: &mut dyn Crystal,
    beam: &dyn Beam,
    wedge: &Wedge,
    angle: f64,
    _angle_num: usize,
    angle_count: usize,
    _energies_per_angle: usize,
    photon_energy: f64,
    cryo: &escape::CryoEscape,
    _pe: &escape::PeEscape,
) -> f64 {
    use rand::Rng;

    let crystal_size = crystal.cryst_size_voxels();
    let crystal_ppm = crystal.crystal_pix_per_um();
    let cryo_ppm = cryo.cryo_pix_per_um;
    let cc = crystal.coefcalc();
    let cryo_abs_coeff = cc.cryo_absorption_coefficient();
    let att_coeff = cc.attenuation_coefficient();
    let density = cc.density();

    if cryo_abs_coeff <= 0.0 {
        return 0.0;
    }

    // Java line 1236: energyPerFluence uses cryo abs coeff / cryo PPM
    let energy_per_fluence = -(-cryo_abs_coeff / cryo_ppm).exp_m1();
    // Java line 1226: fluenceToDoseFactor uses cryo abs coeff / cryo PPM
    // but crystal PPM and crystal density for mass conversion
    let _fluence_to_dose_factor = -(-cryo_abs_coeff / cryo_ppm).exp_m1()
        / (UNIT_CONVERSION * (crystal_ppm.powi(-3) * density))
        * GY_TO_MGY;
    // Java line 1241: energyToDoseFactor uses crystal PPM (dose = energy/mass)
    let energy_to_dose_factor = UNIT_CONVERSION * (crystal_ppm.powi(-3) * density);
    // Java line 1244: beamAttenuationFactor uses cryo PPM
    let beam_attenuation_factor = cryo_ppm.powi(-2) * wedge.total_sec() / angle_count as f64;
    let beam_attenuation_exp_factor = -att_coeff;

    // ppmRatio for converting cryo voxels to crystal voxels (Java line 1252)
    let ppm_ratio = crystal_ppm / cryo_ppm;

    // Crystal bounding box (for depth clamping, Java lines 1280-1290)
    let cryst_size_um = crystal.cryst_size_um();
    let coord_0 = crystal.get_cryst_coord(0, 0, 0);
    let min_dims = coord_0;
    let max_dims = [
        coord_0[0] + cryst_size_um[0],
        coord_0[1] + cryst_size_um[1],
        coord_0[2] + cryst_size_um[2],
    ];

    let wedge_start = wedge.start_vector();
    let wedge_translation = wedge.translation_vector(angle);
    let angle_cos = angle.cos();
    let angle_sin = angle.sin();

    let extra = cryo.cryo_extra_voxels;
    let mut total_dose_back = 0.0;
    let mut rng = rand::thread_rng();

    let beam_energy_j = photon_energy * KEV_TO_JOULES;

    crystal.setup_depth_finding(angle, wedge);

    for ci in 0..cryo.cryo_size_voxels[0] {
        for cj in 0..cryo.cryo_size_voxels[1] {
            for ck in 0..cryo.cryo_size_voxels[2] {
                // Java line 1259: convert cryo voxel to crystal voxel coords via ppmRatio
                let i_cryst = (ci as f64 - extra as f64) * ppm_ratio;
                let j_cryst = (cj as f64 - extra as f64) * ppm_ratio;
                let k_cryst = (ck as f64 - extra as f64) * ppm_ratio;

                // Java line 1263: round to crystal voxel indices
                let ii = i_cryst.round() as i64;
                let jj = j_cryst.round() as i64;
                let kk = k_cryst.round() as i64;

                // Java line 1267: only process voxels OUTSIDE the crystal
                let inside = ii >= 0
                    && jj >= 0
                    && kk >= 0
                    && (ii as usize) < crystal_size[0]
                    && (jj as usize) < crystal_size[1]
                    && (kk as usize) < crystal_size[2]
                    && crystal.is_crystal_at(ii as usize, jj as usize, kk as usize);
                if inside {
                    continue;
                }

                // Java line 1268: get physical coords from cryo grid
                // Cryo coord = crystal_origin + (ci - extra) / cryo_ppm
                let x_um = coord_0[0] + (ci as f64 - extra as f64) / cryo_ppm;
                let y_um = coord_0[1] + (cj as f64 - extra as f64) / cryo_ppm;
                let z_um = coord_0[2] + (ck as f64 - extra as f64) / cryo_ppm;
                let coords = [x_um, y_um, z_um];
                let translated = translate_crystal_to_position(
                    &coords,
                    &wedge_start,
                    &wedge_translation,
                    angle_cos,
                    angle_sin,
                );

                let intensity =
                    beam.beam_intensity(translated[0], translated[1], wedge.off_axis_um);
                if intensity <= 0.0 {
                    continue;
                }

                // Java lines 1280-1290: clamp depth coordinates to crystal bounding box
                let mut depth_coords = translated;
                for dim in 0..3 {
                    if depth_coords[dim] < min_dims[dim] {
                        depth_coords[dim] = min_dims[dim];
                    } else if depth_coords[dim] > max_dims[dim] {
                        depth_coords[dim] = max_dims[dim];
                    }
                }

                // Java line 1292: depth uses clamped coordinates
                let depth = crystal.find_depth(&depth_coords, angle, wedge);

                // Java line 1294: attenuated fluence
                let vox_fluence = intensity
                    * beam_attenuation_factor
                    * (depth * beam_attenuation_exp_factor).exp();

                // Java line 1298-1299
                let _num_photons = vox_fluence / beam_energy_j;

                // Java line 1301: cryo energy absorbed
                let cryo_vox_energy = energy_per_fluence * vox_fluence;

                if cryo_vox_energy <= 0.0 {
                    continue;
                }

                // Java line 1325: PE energy using binding fraction
                // (when FL escape enabled, uses binding fraction approach)
                let binding_frac = if photon_energy > 0.0 {
                    cryo.cryo_energy_to_subtract / photon_energy
                } else {
                    0.0
                };
                let energy_pe = cryo_vox_energy - binding_frac * cryo_vox_energy;

                if energy_pe <= 0.0 {
                    continue;
                }

                // Java line 1329: addDoseAfterPECryo(iCryst, jCryst, kCryst, energyPE, energyToDoseFactor)
                // Track PE from cryo voxel, depositing dose into crystal voxels.
                // Base coords (i_cryst, j_cryst, k_cryst) and relative offsets are both
                // in crystal voxel space (see Java line 1935).
                if cryo.cryo_track_bias.is_empty() {
                    continue;
                }

                let pe_dist_bins = cryo.cryo_pe_distances.len();
                for _ in 0..(escape::pe_angle_resolution() * escape::pe_angle_resolution()) {
                    let random_idx = rng.gen_range(0..cryo.cryo_track_bias.len());
                    let random_track = cryo.cryo_track_bias[random_idx];

                    for m in 0..pe_dist_bins {
                        if random_track >= cryo.relative_vox_xyz_cryo[m].len() {
                            continue;
                        }
                        let [rx, ry, rz] = cryo.relative_vox_xyz_cryo[m][random_track];

                        let partial_energy = energy_pe * cryo.propn_dose_at_dist_cryo[m]
                            / (escape::pe_angle_resolution() * escape::pe_angle_resolution())
                                as f64;
                        let partial_dose = (partial_energy / energy_to_dose_factor) * 1e-6;

                        // i_cryst + rx are both in crystal voxel space
                        let ti = (i_cryst + rx).round() as i64;
                        let tj = (j_cryst + ry).round() as i64;
                        let tk = (k_cryst + rz).round() as i64;

                        if ti >= 0
                            && tj >= 0
                            && tk >= 0
                            && (ti as usize) < crystal_size[0]
                            && (tj as usize) < crystal_size[1]
                            && (tk as usize) < crystal_size[2]
                            && crystal.is_crystal_at(ti as usize, tj as usize, tk as usize)
                        {
                            crystal.add_dose(ti as usize, tj as usize, tk as usize, partial_dose);
                            total_dose_back += partial_dose;
                        }
                    }
                }
            }
        }
    }

    total_dose_back
}

/// Create a crystal from parsed configuration.
pub fn create_crystal(
    config: &crate::parser::config::CrystalConfig,
) -> Result<Box<dyn Crystal>, String> {
    let crystal_type = config.crystal_type.unwrap_or(CrystalType::Cuboid);

    match crystal_type {
        CrystalType::Cuboid => Ok(Box::new(CrystalCuboid::from_config(config)?)),
        CrystalType::Polyhedron => Ok(Box::new(CrystalPolyhedron::from_config(config)?)),
        CrystalType::Cylinder => Ok(Box::new(polyhedron::crystal_cylinder_from_config(config)?)),
        CrystalType::SphericalNew => Ok(Box::new(polyhedron::crystal_spherical_new_from_config(
            config,
        )?)),
        // Match Java: "Type Spherical" maps to CrystalSphericalNew (icosphere mesh),
        // same as "Type SphericalNew". The analytic CrystalSpherical is still available
        // in code but not selectable via input file, matching Java's behavior.
        CrystalType::Spherical => Ok(Box::new(polyhedron::crystal_spherical_new_from_config(
            config,
        )?)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::config::{CoefCalcType, CrystalConfig};

    const CRYSTAL_RESOLUTION_MARKER: f64 = 0.533743110;

    fn default_crystal_config(crystal_type: CrystalType) -> CrystalConfig {
        CrystalConfig {
            crystal_type: Some(crystal_type),
            coefcalc: Some(CoefCalcType::Average),
            dim_x: Some(10.0),
            dim_y: Some(10.0),
            dim_z: Some(10.0),
            pixels_per_micron: Some(CRYSTAL_RESOLUTION_MARKER),
            angle_p: Some(0.0),
            angle_l: Some(0.0),
            cell_a: Some(100.0),
            cell_b: Some(100.0),
            cell_c: Some(100.0),
            ..Default::default()
        }
    }

    #[test]
    fn create_crystal_cuboid() {
        let c = create_crystal(&default_crystal_config(CrystalType::Cuboid)).unwrap();
        assert!(
            (c.crystal_pix_per_um() - CRYSTAL_RESOLUTION_MARKER).abs() < 1e-6,
            "Resolution mismatch"
        );
    }

    #[test]
    fn create_crystal_spherical() {
        let c = create_crystal(&default_crystal_config(CrystalType::Spherical)).unwrap();
        assert!(
            (c.crystal_pix_per_um() - CRYSTAL_RESOLUTION_MARKER).abs() < 1e-6,
            "Resolution mismatch"
        );
    }
}
