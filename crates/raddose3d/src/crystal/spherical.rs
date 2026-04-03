use crate::beam::Beam;
use crate::coefcalc::{self, CoefCalc};
use crate::container::{self, Container, ContainerTransparent};
use crate::ddm::{self, DdmModel};
use crate::output::ExposureSummary;
use crate::parser::config::CrystalConfig;
use crate::wedge::Wedge;

/// Analytic sphere crystal: occupancy by radius check, depth by analytic formula.
/// This matches Java's CrystalSpherical (no polyhedron mesh).
#[derive(Debug)]
#[allow(dead_code)] // phase-6 fields (photo_electron_escape, fluorescent_escape)
pub struct CrystalSpherical {
    diameter: f64,
    pix_per_um: f64,
    cryst_size_voxels: [usize; 3],
    cryst_coord: Vec<[f64; 3]>,
    dose: Vec<f64>,
    fluence: Vec<f64>,
    elastic: Vec<f64>,
    /// Pre-computed occupancy (true = crystal material present).
    crystocc: Vec<bool>,
    coefcalc: Box<dyn CoefCalc>,
    ddm: Box<dyn DdmModel>,
    container: Box<dyn Container>,
    exposure_summary: ExposureSummary,
    subprogram: String,
    photo_electron_escape: bool,
    fluorescent_escape: bool,
}

impl CrystalSpherical {
    const DEFAULT_RESOLUTION: f64 = 0.5;

    pub fn from_config(config: &CrystalConfig) -> Result<Self, String> {
        let diameter = config
            .dim_x
            .ok_or("Spherical crystal requires DimX (diameter)")?;
        let radius = diameter / 2.0;

        let pix_per_um = config.pixels_per_micron.unwrap_or_else(|| {
            let default_res = 10.0 / diameter;
            default_res.min(Self::DEFAULT_RESOLUTION)
        });

        let n = (diameter * pix_per_um).round() as usize + 1;
        let total_voxels = n * n * n;

        let mut cryst_coord = vec![[0.0f64; 3]; total_voxels];
        let mut crystocc = vec![false; total_voxels];

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let idx = i * n * n + j * n + k;
                    let x = -radius + i as f64 / pix_per_um;
                    let y = -radius + j as f64 / pix_per_um;
                    let z = -radius + k as f64 / pix_per_um;
                    cryst_coord[idx] = [x, y, z];
                    // Occupancy: voxel is inside if distance from origin <= radius
                    let dist = (x * x + y * y + z * z).sqrt();
                    crystocc[idx] = dist <= radius;
                }
            }
        }

        let coefcalc = coefcalc::create_coefcalc(config)?;
        let ddm = ddm::create_ddm(
            config.ddm,
            config.gamma_param,
            config.b0_param,
            config.beta_param,
        );

        let subprogram = config.program.as_deref().unwrap_or("RD3D").to_uppercase();
        match subprogram.as_str() {
            "RD3D" | "" => {}
            "MONTECARLO" | "GOS" | "XFEL" | "EMSP" | "EMED" => {
                return Err(format!(
                    "Subprogram '{subprogram}' is not supported for Spherical crystals; use Type Cuboid"
                ));
            }
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

        let container = container::create_container(config);

        Ok(CrystalSpherical {
            diameter,
            pix_per_um,
            cryst_size_voxels: [n, n, n],
            cryst_coord,
            dose: vec![0.0; total_voxels],
            fluence: vec![0.0; total_voxels],
            elastic: vec![0.0; total_voxels],
            crystocc,
            coefcalc,
            ddm,
            container,
            exposure_summary: ExposureSummary::new(),
            subprogram,
            photo_electron_escape,
            fluorescent_escape,
        })
    }

    #[inline]
    fn voxel_idx(&self, i: usize, j: usize, k: usize) -> usize {
        let [_, ny, nz] = self.cryst_size_voxels;
        i * ny * nz + j * nz + k
    }
}

impl super::Crystal for CrystalSpherical {
    /// No-op: analytic sphere doesn't need depth finding setup.
    fn setup_depth_finding(&mut self, _ang_rad: f64, _wedge: &Wedge) {}

    fn find_depth(&self, vox_coord: &[f64; 3], delta_phi: f64, wedge: &Wedge) -> f64 {
        let delta = delta_phi - wedge.start_ang;

        // Map voxel coordinates relative to beam position
        let mapped_x = vox_coord[0] - (wedge.start_x + wedge.trans_x * delta);
        let mapped_y = vox_coord[1] - (wedge.start_y + wedge.trans_y * delta);
        let mapped_z = vox_coord[2];

        let radius = self.diameter / 2.0;
        let r_sq = radius * radius;
        let xy_sq = mapped_x * mapped_x + mapped_y * mapped_y;

        if xy_sq > r_sq {
            // Outside sphere in X-Y plane, no depth
            return 0.0;
        }

        let surface_z = (r_sq - xy_sq).sqrt();
        surface_z - mapped_z
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
        self.crystocc[self.voxel_idx(i, j, k)]
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
            "Spherical crystal of diameter {:.0} um at a resolution of {:.2} microns per voxel edge.",
            self.diameter,
            1.0 / self.pix_per_um
        )
    }

    fn cryst_size_voxels(&self) -> [usize; 3] {
        self.cryst_size_voxels
    }

    fn cryst_size_um(&self) -> [f64; 3] {
        [self.diameter, self.diameter, self.diameter]
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

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::Crystal;
    use crate::parser::config::{CrystalConfig, CrystalType, WedgeConfig};
    use crate::wedge::Wedge;

    fn make_spherical_config(diameter: f64) -> CrystalConfig {
        CrystalConfig {
            crystal_type: Some(CrystalType::Spherical),
            dim_x: Some(diameter),
            // Minimal unit cell for CoefCalcFromParams
            cell_a: Some(78.0),
            cell_b: Some(78.0),
            cell_c: Some(37.0),
            num_residues: Some(60),
            ..Default::default()
        }
    }

    fn default_wedge() -> Wedge {
        Wedge::from_config(&WedgeConfig {
            start_ang: 0.0,
            end_ang: 0.0,
            ..Default::default()
        })
    }

    #[test]
    fn test_spherical_occupancy() {
        let diameter = 100.0f64;
        let config = make_spherical_config(diameter);
        let crystal = CrystalSpherical::from_config(&config).expect("Spherical should construct");

        // Center voxel should be inside
        let [nx, ny, nz] = crystal.cryst_size_voxels;
        let ci = nx / 2;
        let cj = ny / 2;
        let ck = nz / 2;
        assert!(
            crystal.is_crystal_at(ci, cj, ck),
            "Center voxel should be inside sphere"
        );

        // Corner voxel should be outside (distance = sqrt(3) * radius > radius)
        assert!(
            !crystal.is_crystal_at(0, 0, 0),
            "Corner voxel (0,0,0) should be outside sphere"
        );
    }

    #[test]
    fn test_spherical_depth_at_center() {
        let diameter = 100.0f64;
        let config = make_spherical_config(diameter);
        let crystal = CrystalSpherical::from_config(&config).expect("Spherical should construct");

        let wedge = default_wedge();
        // At center (0, 0, 0), mapped_x = 0, mapped_y = 0
        // surface_z = radius = 50.0, mapped_z = 0
        // depth = 50.0 - 0.0 = 50.0
        let depth = crystal.find_depth(&[0.0, 0.0, 0.0], 0.0, &wedge);
        assert!(
            (depth - diameter / 2.0).abs() < 0.01,
            "Depth at center should be radius ({:.1}), got {:.4}",
            diameter / 2.0,
            depth
        );
    }

    #[test]
    fn test_spherical_depth_outside() {
        let diameter = 100.0f64;
        let config = make_spherical_config(diameter);
        let crystal = CrystalSpherical::from_config(&config).expect("Spherical should construct");

        let wedge = default_wedge();
        // Far outside sphere in X
        let depth = crystal.find_depth(&[200.0, 0.0, 0.0], 0.0, &wedge);
        assert_eq!(depth, 0.0, "Depth outside sphere should be 0");
    }
}
