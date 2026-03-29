pub mod cuboid;
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

    // Main loop: angles × energies × voxels
    for (n, &angle) in angles.iter().enumerate() {
        for &photon_energy in &energies {
            // Update coefficients for this energy
            crystal.coefcalc_mut().update_coefficients(photon_energy);

            expose_angle(
                crystal,
                beam,
                wedge,
                angle,
                n,
                angles.len(),
                energies_per_angle,
                photon_energy,
            );
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

/// Expose one angle of the crystal.
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
) {
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
                    let rde = crystal.ddm().calc_decay(total_dose_before);
                    let energy_per_fluence = -(-abs_coeff / pix_per_um).exp_m1();
                    let energy_absorbed =
                        energy_per_fluence * vox_fluence + energy_per_fluence * compton_fluence;

                    // Notify exposure summary
                    crystal.exposure_summary_mut().exposure_observation(
                        angle_num,
                        i,
                        j,
                        k,
                        vox_dose,
                        total_dose_before,
                        vox_fluence,
                        rde,
                        energy_absorbed,
                        elastic_yield,
                        angle_count as f64,
                    );

                    crystal.add_fluence(i, j, k, vox_fluence);

                    // No escape: add full dose to voxel
                    crystal.add_dose(i, j, k, vox_dose);
                    crystal.add_dose(i, j, k, vox_dose_compton);
                    crystal.add_elastic(i, j, k, elastic_yield);
                } else if vox_dose < 0.0 {
                    panic!("Negative dose encountered - this should never happen");
                }
            }
        }
    }
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
        CrystalType::Spherical => Ok(Box::new(CrystalSpherical::from_config(config)?)),
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
