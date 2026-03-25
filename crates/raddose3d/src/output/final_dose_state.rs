use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

/// Outputs the final per-voxel dose, fluence, and elastic scattering as CSV.
///
/// Port of Java `OutputFinalDoseStateCSV`.
/// Format: x,y,z,dose,fluence,elastic  (one line per occupied voxel)
#[derive(Debug)]
pub struct OutputFinalDoseStateCSV {
    writer: DebugWriter,
}

impl OutputFinalDoseStateCSV {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputFinalDoseStateCSV {
            writer: DebugWriter(writer),
        }
    }

    fn write_crystal(&mut self, crystal: &dyn Crystal) {
        let size = crystal.cryst_size_voxels();
        for i in 0..size[0] {
            for j in 0..size[1] {
                for k in 0..size[2] {
                    if crystal.is_crystal_at(i, j, k) {
                        let coord = crystal.get_cryst_coord(i, j, k);
                        let dose = crystal.get_dose(i, j, k) as f32;
                        let fluence = crystal.get_fluence(i, j, k) as f32;
                        let elastic = crystal.get_elastic(i, j, k) as f32;

                        let dose_str = if dose <= f32::MIN_POSITIVE {
                            "0".to_string()
                        } else {
                            dose.to_string()
                        };
                        let fluence_str = if fluence <= f32::MIN_POSITIVE {
                            "0".to_string()
                        } else {
                            fluence.to_string()
                        };
                        let elastic_str = if elastic <= f32::MIN_POSITIVE {
                            "0".to_string()
                        } else {
                            elastic.to_string()
                        };

                        let _ = writeln!(
                            self.writer,
                            "{},{},{},{},{},{}",
                            coord[0] as f32,
                            coord[1] as f32,
                            coord[2] as f32,
                            dose_str,
                            fluence_str,
                            elastic_str,
                        );
                    }
                }
            }
        }
    }
}

impl super::Output for OutputFinalDoseStateCSV {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(&mut self, _wedge: &Wedge, _summary: &ExposureSummary) {}

    fn close(&mut self, crystal: Option<&dyn Crystal>) {
        match crystal {
            Some(c) => self.write_crystal(c),
            None => {
                let _ = writeln!(
                    self.writer,
                    "OutputFinalDoseStateCSV: No crystal object has been seen."
                );
            }
        }
        let _ = self.writer.flush();
    }
}

/// Generates an R script for 3D dose visualisation using the `rgl` library.
///
/// Port of Java `OutputFinalDoseStateR`.
#[derive(Debug)]
pub struct OutputFinalDoseStateR {
    writer: DebugWriter,
}

impl OutputFinalDoseStateR {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputFinalDoseStateR {
            writer: DebugWriter(writer),
        }
    }

    fn write_r_script(&mut self, crystal: &dyn Crystal) {
        let csize = crystal.cryst_size_voxels();

        let _ = writeln!(self.writer, "# Crystal dose state visualization using R");
        let _ = writeln!(self.writer, "# http://www.r-project.org/");
        let _ = writeln!(self.writer, "#");
        let _ = writeln!(
            self.writer,
            "# Crystal size: {} x {} x {} voxels",
            csize[0], csize[1], csize[2]
        );
        let _ = writeln!(self.writer, "#\n");

        let _ = writeln!(self.writer, "contourlevels <- c(0.1, 20, 30) # MGy");
        let _ = writeln!(
            self.writer,
            "contourcolours <- c('lightblue', 'darkblue', 'red')"
        );
        let _ = writeln!(self.writer, "contouropacity <- c(0.2, 0.5, 1)\n");

        let _ = writeln!(self.writer, "require(\"rgl\")");
        let _ = writeln!(self.writer, "require(\"misc3d\")\n");

        let _ = writeln!(self.writer, "# Three dimensional dose array (MGy)");
        let _ = writeln!(
            self.writer,
            "dose <- array(0, c({}, {}, {}))",
            csize[0], csize[1], csize[2]
        );

        for k in 0..csize[2] {
            let _ = write!(self.writer, "dose[,,{}]<-c(", k + 1);
            let mut first = true;
            for j in 0..csize[1] {
                for i in 0..csize[0] {
                    if !first {
                        let _ = write!(self.writer, ",");
                    }
                    first = false;
                    let dose = crystal.get_dose(i, j, k) as f32;
                    if dose <= f32::MIN_POSITIVE {
                        let _ = write!(self.writer, "0");
                    } else {
                        let _ = write!(self.writer, "{}", dose);
                    }
                }
            }
            let _ = writeln!(self.writer, ")");
        }

        let _ = writeln!(
            self.writer,
            "contour3d(dose, level=contourlevels, color=contourcolours, alpha=contouropacity)"
        );
        let _ = writeln!(self.writer, "# axes3d()");
        let _ = writeln!(
            self.writer,
            "wire3d(translate3d(scale3d(cube3d(),{}/2,{}/2,{}/2),{}/2,{}/2,{}/2),col = 'grey')",
            csize[0], csize[1], csize[2], csize[0], csize[1], csize[2]
        );
        for (theta, phi) in [(0, 0), (0, 45), (45, 45), (90, 45), (90, 90)] {
            let _ = writeln!(
                self.writer,
                "rgl.viewpoint( theta = {}, phi = {})",
                theta, phi
            );
            let _ = writeln!(
                self.writer,
                "rgl.snapshot( \"plot_{}_{}.png\", fmt = \"png\", top = TRUE)",
                theta, phi
            );
            let _ = writeln!(self.writer, "Sys.sleep(1)");
        }
        let _ = writeln!(self.writer, "print(\"Plots Saved\")");
    }
}

impl super::Output for OutputFinalDoseStateR {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(&mut self, _wedge: &Wedge, _summary: &ExposureSummary) {}

    fn close(&mut self, crystal: Option<&dyn Crystal>) {
        match crystal {
            Some(c) => self.write_r_script(c),
            None => {
                let _ = writeln!(
                    self.writer,
                    "# OutputFinalDoseStateR: No crystal object has been seen."
                );
            }
        }
        let _ = self.writer.flush();
    }
}
