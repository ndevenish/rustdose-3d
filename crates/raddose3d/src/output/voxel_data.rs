use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

/// Outputs per-image, per-voxel dose values as CSV.
///
/// Port of Java `OutputVoxelDose`.
/// Row 0: voxel coordinate labels (x_y_z,...).
/// Rows 1..N: cumulative dose for each occupied voxel at each image.
#[derive(Debug)]
pub struct OutputVoxelDose {
    writer: DebugWriter,
    wedge: Option<Wedge>,
}

impl OutputVoxelDose {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputVoxelDose {
            writer: DebugWriter(writer),
            wedge: None,
        }
    }
}

impl super::Output for OutputVoxelDose {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(&mut self, wedge: &Wedge, _summary: &ExposureSummary) {
        self.wedge = Some(wedge.clone());
    }

    fn close(&mut self, crystal: Option<&dyn Crystal>) {
        let crystal = match crystal {
            Some(c) => c,
            None => {
                let _ = writeln!(
                    self.writer,
                    "OutputVoxelDose: No crystal object has been seen."
                );
                let _ = self.writer.flush();
                return;
            }
        };

        let size = crystal.cryst_size_voxels();

        // Header row: coordinate labels for occupied voxels
        for i in 0..size[0] {
            for j in 0..size[1] {
                for k in 0..size[2] {
                    if crystal.is_crystal_at(i, j, k) {
                        let coord = crystal.get_cryst_coord(i, j, k);
                        let _ = write!(
                            self.writer,
                            "{}_{}_{},",
                            coord[0] as f32, coord[1] as f32, coord[2] as f32
                        );
                    }
                }
            }
        }
        let _ = writeln!(self.writer);

        // Per-image dose rows from ExposureSummary snapshot
        let summary = crystal.exposure_summary();
        let vox_dose_image = &summary.vox_dose_image;

        let num_images = match &self.wedge {
            Some(w) => crystal.num_images(w),
            None => vox_dose_image.len(),
        };

        for l in 0..num_images {
            if l >= vox_dose_image.len() {
                break;
            }
            let snap = &vox_dose_image[l];
            for i in 0..size[0] {
                for j in 0..size[1] {
                    for k in 0..size[2] {
                        if crystal.is_crystal_at(i, j, k) {
                            let idx = i * size[1] * size[2] + j * size[2] + k;
                            let dose = if idx < snap.len() {
                                snap[idx] as f32
                            } else {
                                0.0
                            };
                            let _ = write!(self.writer, "{},", dose);
                        }
                    }
                }
            }
            let _ = writeln!(self.writer);
        }

        self.wedge = None;
        let _ = self.writer.flush();
    }
}

/// Outputs per-image, per-voxel fluence values as CSV.
///
/// Port of Java `OutputVoxelFluences`.
#[derive(Debug)]
pub struct OutputVoxelFluences {
    writer: DebugWriter,
    wedge: Option<Wedge>,
}

impl OutputVoxelFluences {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputVoxelFluences {
            writer: DebugWriter(writer),
            wedge: None,
        }
    }
}

impl super::Output for OutputVoxelFluences {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(&mut self, wedge: &Wedge, _summary: &ExposureSummary) {
        self.wedge = Some(wedge.clone());
    }

    fn close(&mut self, crystal: Option<&dyn Crystal>) {
        let crystal = match crystal {
            Some(c) => c,
            None => {
                let _ = writeln!(
                    self.writer,
                    "OutputVoxelFluences: No crystal object has been seen."
                );
                let _ = self.writer.flush();
                return;
            }
        };

        let size = crystal.cryst_size_voxels();

        // Header row: coordinate labels for occupied voxels
        for i in 0..size[0] {
            for j in 0..size[1] {
                for k in 0..size[2] {
                    if crystal.is_crystal_at(i, j, k) {
                        let coord = crystal.get_cryst_coord(i, j, k);
                        let _ = write!(
                            self.writer,
                            "{}_{}_{},",
                            coord[0] as f32, coord[1] as f32, coord[2] as f32
                        );
                    }
                }
            }
        }
        let _ = writeln!(self.writer);

        let summary = crystal.exposure_summary();
        let vox_fluence_image = &summary.vox_fluence_image;

        let num_images = match &self.wedge {
            Some(w) => crystal.num_images(w),
            None => vox_fluence_image.len(),
        };

        for l in 0..num_images {
            if l >= vox_fluence_image.len() {
                break;
            }
            let snap = &vox_fluence_image[l];
            for i in 0..size[0] {
                for j in 0..size[1] {
                    for k in 0..size[2] {
                        if crystal.is_crystal_at(i, j, k) {
                            let idx = i * size[1] * size[2] + j * size[2] + k;
                            let flu = if idx < snap.len() {
                                snap[idx] as f32
                            } else {
                                0.0
                            };
                            let _ = write!(self.writer, "{},", flu);
                        }
                    }
                }
            }
            let _ = writeln!(self.writer);
        }

        self.wedge = None;
        let _ = self.writer.flush();
    }
}
