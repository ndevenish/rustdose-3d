use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

/// Outputs per-image angle, fluence-weighted average RDE, and minimum RDE as CSV.
///
/// Port of Java `OutputRDECSV`.
/// Header: `Image Number, Angle, Avg RDE, min RDE`
#[derive(Debug)]
pub struct OutputRDECSV {
    writer: DebugWriter,
}

impl OutputRDECSV {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        let mut w = DebugWriter(writer);
        let _ = writeln!(w, "Image Number, Angle, Avg RDE, min RDE");
        OutputRDECSV { writer: w }
    }
}

impl super::Output for OutputRDECSV {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(&mut self, _wedge: &Wedge, summary: &ExposureSummary) {
        let weighted = summary.fluence_weighted_rde_array();
        let min_arr = summary.min_rde_array();

        for i in 0..weighted.len() {
            let image = i + 1;
            let angle_deg = weighted[i][0].to_degrees();
            let avg_rde = weighted[i][1];
            let min_rde = if i < min_arr.len() {
                min_arr[i][1]
            } else {
                0.0
            };
            let _ = writeln!(
                self.writer,
                "{},{},{},{}",
                image, angle_deg, avg_rde, min_rde
            );
        }
    }

    fn close(&mut self, _crystal: Option<&dyn Crystal>) {
        let _ = self.writer.flush();
    }
}
