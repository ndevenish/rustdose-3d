use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

/// Outputs per-image DWD, angle, volume, and RDE values at 5 resolution shells.
///
/// Port of Java `OutputDWDs`.
#[derive(Debug)]
pub struct OutputDWDs {
    writer: DebugWriter,
}

impl OutputDWDs {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        let mut w = DebugWriter(writer);
        let _ = writeln!(
            w,
            "RADDOSE Image Number, DWD Angle, DWD, Vol, 1A RDE, 2A RDE, 3A RDE, 4A RDE, max res RDE"
        );
        OutputDWDs { writer: w }
    }
}

impl super::Output for OutputDWDs {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(&mut self, _wedge: &Wedge, summary: &ExposureSummary) {
        let image_dwd = summary.image_dwds();
        let angle_dwd = summary.angle_dwd_array();
        let image_vol = summary.image_vol_array();
        let image_rde = summary.image_rde_array();

        for i in 0..angle_dwd.len() {
            let image = i + 1;
            let angle_deg = angle_dwd[i].to_degrees();
            let dwd = if i < image_dwd.len() {
                image_dwd[i]
            } else {
                0.0
            };
            let vol = if i < image_vol.len() {
                image_vol[i]
            } else {
                0.0
            };
            let rde = if i < image_rde.len() {
                image_rde[i]
            } else {
                [0.0; 5]
            };
            // Java order: image, angle, DWD, vol, 1A(idx1), 2A(idx2), 3A(idx3), 4A(idx4), max_res(idx0)
            let _ = writeln!(
                self.writer,
                "{},{},{},{},{},{},{},{},{}",
                image, angle_deg, dwd, vol, rde[1], rde[2], rde[3], rde[4], rde[0]
            );
        }
    }

    fn close(&mut self, _crystal: Option<&dyn Crystal>) {
        let _ = self.writer.flush();
    }
}
