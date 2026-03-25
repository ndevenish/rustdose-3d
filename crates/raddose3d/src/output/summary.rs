use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::ExposureSummary;
use crate::wedge::Wedge;

/// Newtype wrapper around Write that implements Debug.
struct DebugWriter(Box<dyn Write + Send + Sync>);

impl std::fmt::Debug for DebugWriter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Writer")
    }
}

impl std::ops::Deref for DebugWriter {
    type Target = dyn Write;
    fn deref(&self) -> &Self::Target {
        &*self.0
    }
}

impl std::ops::DerefMut for DebugWriter {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut *self.0
    }
}

/// OutputSummaryText: produces human-readable dose summary.
#[derive(Debug)]
pub struct OutputSummaryText {
    writer: DebugWriter,
    wedge_num: usize,
}

impl OutputSummaryText {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputSummaryText {
            writer: DebugWriter(writer),
            wedge_num: 0,
        }
    }
}

impl super::Output for OutputSummaryText {
    fn publish_crystal(&mut self, crystal: &dyn Crystal) {
        let _ = writeln!(self.writer, "{}", crystal.crystal_info());
    }

    fn publish_beam(&mut self, beam: &dyn Beam) {
        let _ = writeln!(self.writer, "{}", beam.description());
    }

    fn publish_wedge(&mut self, wedge: &Wedge, summary: &ExposureSummary) {
        self.wedge_num += 1;
        let _ = writeln!(self.writer, "Wedge {}:", self.wedge_num);
        let _ = writeln!(self.writer, "{}", wedge.description());

        let abs_en_threshold = 0.5;

        let _ = writeln!(
            self.writer,
            "{:<42}: {:.6} MGy",
            "Average Diffraction Weighted Dose",
            summary.avg_diffracted_dose()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.6} MGy",
            "Last Diffraction Weighted Dose",
            summary.last_dwd()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.2e} photons",
            "Elastic Yield",
            summary.wedge_elastic()
        );

        let dwd = summary.avg_diffracted_dose();
        if dwd > 0.0 {
            let _ = writeln!(
                self.writer,
                "{:<42}: {:.2e} photons/MGy",
                "Diffraction Efficiency (Elastic Yield/DWD)",
                summary.wedge_elastic() / dwd
            );
        }

        let _ = writeln!(
            self.writer,
            "{:<42}: {:.6} MGy",
            "Average Dose (Whole Crystal)",
            summary.avg_dose_whole_crystal()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.6} MGy",
            "Average Dose (Exposed Region)",
            summary.avg_dose_exposed_region()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.6} MGy",
            "Max Dose",
            summary.max_dose()
        );

        let _ = writeln!(
            self.writer,
            "Average Dose ({:.1} % of total absorbed energy threshold ({:.2} MGy)): {:.6} MGy",
            abs_en_threshold * 100.0,
            summary.abs_dose_threshold(abs_en_threshold),
            summary.avg_dose_threshold(abs_en_threshold)
        );

        let _ = writeln!(
            self.writer,
            "{:<42}: {:.2}",
            "Dose Contrast (Max/Threshold Av.)",
            summary.dose_contrast(abs_en_threshold)
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.1}%",
            "Used Volume",
            summary.used_volume_fraction()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.2e} J",
            "Absorbed Energy (this Wedge)",
            summary.abs_energy_total()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.1} 1/g",
            "Dose Inefficiency (Max Dose/mJ Absorbed)",
            summary.dose_inefficiency()
        );
    }

    fn close(&mut self) {
        let _ = self.writer.flush();
    }
}

/// OutputSummaryCSV: produces machine-readable CSV dose summary.
#[derive(Debug)]
pub struct OutputSummaryCSV {
    writer: DebugWriter,
    header_written: bool,
    wedge_num: usize,
}

impl OutputSummaryCSV {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputSummaryCSV {
            writer: DebugWriter(writer),
            header_written: false,
            wedge_num: 0,
        }
    }
}

impl super::Output for OutputSummaryCSV {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(&mut self, _wedge: &Wedge, summary: &ExposureSummary) {
        if !self.header_written {
            let _ = writeln!(
                self.writer,
                "Wedge,AvgDWD,AvgDoseWhole,AvgDoseExposed,MaxDose,ElasticYield,AbsEnergy,UsedVol%"
            );
            self.header_written = true;
        }
        self.wedge_num += 1;
        let _ = writeln!(
            self.writer,
            "{},{},{},{},{},{},{},{}",
            self.wedge_num,
            summary.avg_diffracted_dose(),
            summary.avg_dose_whole_crystal(),
            summary.avg_dose_exposed_region(),
            summary.max_dose(),
            summary.wedge_elastic(),
            summary.abs_energy_total(),
            summary.used_volume_fraction(),
        );
    }

    fn close(&mut self) {
        let _ = self.writer.flush();
    }
}
