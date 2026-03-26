use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

/// Prints a simple ASCII progress bar during simulation.
///
/// Port of Java `OutputProgressIndicator`.
/// Implements both `Output` and `ExposeObserver`.
#[derive(Debug)]
pub struct OutputProgressIndicator {
    writer: DebugWriter,
    image_count: usize,
    wedge_progress: usize,
}

impl OutputProgressIndicator {
    const DOT_INTERVAL: usize = 4;

    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputProgressIndicator {
            writer: DebugWriter(writer),
            image_count: 0,
            wedge_progress: 0,
        }
    }
}

impl super::Output for OutputProgressIndicator {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(
        &mut self,
        _wedge: &Wedge,
        _summary: &ExposureSummary,
        _crystal: Option<&dyn Crystal>,
    ) {
    }

    fn close(&mut self, _crystal: Option<&dyn Crystal>) {
        let _ = self.writer.flush();
    }
}

impl super::ExposeObserver for OutputProgressIndicator {
    fn register(&mut self, _crystal: &dyn Crystal) {}

    fn exposure_start(&mut self, wedge_images: usize, _wedge: &Wedge, _crystal_size: [usize; 3]) {
        let _ = write!(self.writer, "Exposing wedge: [ 0%");
        let _ = self.writer.flush();
        self.image_count = wedge_images;
        self.wedge_progress = 0;
    }

    #[allow(clippy::too_many_arguments)]
    fn exposure_observation(
        &mut self,
        _wedge_image: usize,
        _i: usize,
        _j: usize,
        _k: usize,
        _added_dose: f64,
        _total_dose: f64,
        _fluence: f64,
        _rde: f64,
        _absorbed_energy: f64,
        _elastic: f64,
        _angle_count: f64,
    ) {
    }

    fn image_complete(&mut self, image: usize, _ang_rad: f64, _last_angle: f64, _vox_vol: f64) {
        if self.image_count == 0 {
            return;
        }
        // Advance progress to match (image+1) / image_count * 100
        loop {
            let next_progress = 100 * (image + 1) / self.image_count;
            if next_progress <= self.wedge_progress {
                break;
            }
            self.wedge_progress += 1;

            if self.wedge_progress.is_multiple_of(Self::DOT_INTERVAL) {
                let _ = write!(self.writer, ".");
            }
            match self.wedge_progress {
                20 => {
                    let _ = write!(self.writer, "20%");
                }
                40 => {
                    let _ = write!(self.writer, "40%");
                }
                60 => {
                    let _ = write!(self.writer, "60%");
                }
                80 => {
                    let _ = write!(self.writer, "80%");
                }
                100 => {
                    let _ = write!(self.writer, "100%");
                }
                _ => {}
            }
        }
        let _ = self.writer.flush();
    }

    fn summary_observation(
        &mut self,
        _i: usize,
        _j: usize,
        _k: usize,
        _total_dose: f64,
        _voxel_mass_kg: f64,
    ) {
    }

    fn exposure_complete(&mut self) {
        let _ = writeln!(self.writer, " ]");
        let _ = self.writer.flush();
    }
}

/// Estimates simulation progress by reporting crystal voxel count and wedge slice count.
///
/// Port of Java `OutputProgressEstimate`.
#[derive(Debug)]
pub struct OutputProgressEstimate {
    writer: Option<DebugWriter>,
    crystal_voxels: u64,
    sum_wedge_slices: u64,
    crystal_voxel_list: Vec<u64>,
    sum_wedge_slice_list: Vec<u64>,
}

impl OutputProgressEstimate {
    pub fn new(writer: Option<Box<dyn Write + Send + Sync>>) -> Self {
        OutputProgressEstimate {
            writer: writer.map(DebugWriter),
            crystal_voxels: 0,
            sum_wedge_slices: 0,
            crystal_voxel_list: Vec::new(),
            sum_wedge_slice_list: Vec::new(),
        }
    }

    fn record_estimate(&mut self, voxels: u64, slices: u64) {
        if let Some(ref mut w) = self.writer {
            let _ = writeln!(
                w,
                "Progress Estimate Information: CrystalVoxels={} WedgeSlices={}",
                voxels, slices
            );
        }
        self.crystal_voxel_list.push(voxels);
        self.sum_wedge_slice_list.push(slices);
    }

    pub fn crystal_voxel_list(&self) -> &[u64] {
        &self.crystal_voxel_list
    }

    pub fn wedge_slice_list(&self) -> &[u64] {
        &self.sum_wedge_slice_list
    }
}

impl super::Output for OutputProgressEstimate {
    fn publish_crystal(&mut self, crystal: &dyn Crystal) {
        if self.crystal_voxels > 0 || self.sum_wedge_slices > 0 {
            let v = self.crystal_voxels;
            let s = self.sum_wedge_slices;
            self.record_estimate(v, s);
        }
        let size = crystal.cryst_size_voxels();
        self.crystal_voxels = (size[0] * size[1] * size[2]) as u64;
        self.sum_wedge_slices = 0;
    }

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(
        &mut self,
        wedge: &Wedge,
        _summary: &ExposureSummary,
        _crystal: Option<&dyn Crystal>,
    ) {
        let diff = (wedge.start_ang - wedge.end_ang).abs();
        if diff < wedge.ang_res {
            self.sum_wedge_slices += crate::constants::STATIC_EXPOSURE as u64;
        } else {
            self.sum_wedge_slices += (diff / wedge.ang_res + 1.0) as u64;
        }
    }

    fn close(&mut self, _crystal: Option<&dyn Crystal>) {
        let v = self.crystal_voxels;
        let s = self.sum_wedge_slices;
        self.record_estimate(v, s);
        if let Some(ref mut w) = self.writer {
            let _ = w.flush();
        }
    }
}
