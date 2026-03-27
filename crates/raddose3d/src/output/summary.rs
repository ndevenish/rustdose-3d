use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

const ABS_EN_THRESHOLD: f64 = 0.95;

// ── OutputSummaryText ────────────────────────────────────────────────────────

/// Outputs a human-readable dose summary, matching Java `OutputSummaryText`.
#[derive(Debug)]
pub struct OutputSummaryText {
    writer: DebugWriter,
    wedge_num: usize,
    /// Coefficient description captured at publish_crystal time.
    coefcalc_desc: Option<String>,
}

impl OutputSummaryText {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputSummaryText {
            writer: DebugWriter(writer),
            wedge_num: 0,
            coefcalc_desc: None,
        }
    }
}

impl super::Output for OutputSummaryText {
    fn publish_crystal(&mut self, crystal: &dyn Crystal) {
        let _ = writeln!(self.writer, "{}", crystal.crystal_info());
        let _ = writeln!(self.writer, "{}.", crystal.ddm().name());
    }

    fn publish_beam(&mut self, beam: &dyn Beam) {
        let _ = writeln!(self.writer, "{}", beam.description());
    }

    fn publish_wedge(
        &mut self,
        wedge: &Wedge,
        summary: &ExposureSummary,
        crystal: Option<&dyn Crystal>,
    ) {
        // Capture coefcalc description now (after expose, so coefficients are non-zero).
        if let Some(c) = crystal {
            self.coefcalc_desc = Some(c.coefcalc().description());
        }
        self.wedge_num += 1;
        let _ = writeln!(self.writer, "Wedge {}:", self.wedge_num);
        let _ = write!(self.writer, "{}", wedge.description());

        if let Some(ref desc) = self.coefcalc_desc {
            let _ = write!(self.writer, "{}", desc);
        }
        let _ = writeln!(self.writer);

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

        let threshold_label = format!(
            "Average Dose ({:.1} % of total absorbed energy threshold ({:.2} MGy))",
            ABS_EN_THRESHOLD * 100.0,
            summary.abs_dose_threshold(ABS_EN_THRESHOLD),
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.6} MGy",
            threshold_label,
            summary.avg_dose_threshold(ABS_EN_THRESHOLD)
        );

        let _ = writeln!(
            self.writer,
            "{:<42}: {:.2}",
            "Dose Contrast (Max/Threshold Av.)",
            summary.dose_contrast(ABS_EN_THRESHOLD)
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.1}%",
            "Used Volume",
            summary.used_volume_fraction()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.2e} J.",
            "Absorbed Energy (this Wedge)",
            summary.abs_energy_total()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.1} 1/g",
            "Dose Inefficiency (Max Dose/mJ Absorbed)",
            summary.dose_inefficiency()
        );
        let _ = writeln!(
            self.writer,
            "{:<42}: {:.1} 1/g",
            "Dose Inefficiency PE (Max Dose/mJ Deposited)",
            summary.dose_inefficiency_pe()
        );

        // Final dose histogram from ExposureSummary
        let normed = summary.dose_hist_normalised();
        let bins = normed.len(); // 11
        let _ = writeln!(self.writer, "Final Dose Histogram:");
        // Bin 1: underflow (< 0.1 MGy)
        let _ = writeln!(
            self.writer,
            "Bin  1,  0.0 to {:4.1} MGy: {:4.1} % ",
            ExposureSummary::dose_hist_break(1),
            normed[0] * 100.0
        );
        #[allow(clippy::needless_range_loop)]
        for i in 1..bins - 1 {
            let _ = writeln!(
                self.writer,
                "Bin {:2}, {:4.1} to {:4.1} MGy: {:4.1} % ",
                i + 1,
                ExposureSummary::dose_hist_break(i),
                ExposureSummary::dose_hist_break(i + 1),
                normed[i] * 100.0
            );
        }
        // Overflow bin
        let _ = writeln!(
            self.writer,
            "Bin {:2}, {:4.1} MGy upwards: {:4.1} %",
            bins,
            ExposureSummary::dose_hist_break(bins - 1),
            normed[bins - 1] * 100.0
        );
    }

    fn close(&mut self, _crystal: Option<&dyn Crystal>) {
        let _ = self.writer.flush();
    }
}

// ── OutputSummaryCSV ─────────────────────────────────────────────────────────

/// Outputs a machine-readable CSV dose summary, matching Java `OutputSummaryCSV`.
///
/// Header:
/// `Wedge Number, Average DWD, Last DWD, Elastic Yield (wedge),
///  Diffraction Efficiency, AD-WC, AD-ExpRegion, Max Dose,
///  Dose Threshold, Abs En Threshold, TAD, Dose Contrast,
///  Used Volume, Wedge Absorbed Energy, Dose Inefficiency, Dose Inefficiency PE`
#[derive(Debug)]
pub struct OutputSummaryCSV {
    writer: DebugWriter,
    wedge_num: usize,
}

impl OutputSummaryCSV {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        let mut w = DebugWriter(writer);
        let _ = writeln!(
            w,
            "Wedge Number, Average DWD, Last DWD,Elastic Yield (wedge), \
             Diffraction Efficiency, AD-WC, AD-ExpRegion, Max Dose, \
             Dose Threshold, Abs En Threshold, TAD, Dose Contrast, \
             Used Volume, Wedge Absorbed Energy, Dose Inefficiency, Dose Inefficiency PE"
        );
        OutputSummaryCSV {
            writer: w,
            wedge_num: 0,
        }
    }
}

impl super::Output for OutputSummaryCSV {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {}

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(
        &mut self,
        _wedge: &Wedge,
        summary: &ExposureSummary,
        _crystal: Option<&dyn Crystal>,
    ) {
        self.wedge_num += 1;
        let dwd = summary.avg_diffracted_dose();
        let diff_eff = if dwd > 0.0 {
            summary.wedge_elastic() / dwd
        } else {
            0.0
        };
        let _ = writeln!(
            self.writer,
            "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}",
            self.wedge_num,
            dwd,
            summary.last_dwd(),
            summary.wedge_elastic(),
            diff_eff,
            summary.avg_dose_whole_crystal(),
            summary.avg_dose_exposed_region(),
            summary.max_dose(),
            summary.abs_dose_threshold(ABS_EN_THRESHOLD),
            ABS_EN_THRESHOLD * 100.0,
            summary.avg_dose_threshold(ABS_EN_THRESHOLD),
            summary.dose_contrast(ABS_EN_THRESHOLD),
            summary.used_volume_fraction(),
            summary.abs_energy_total(),
            summary.dose_inefficiency(),
            summary.dose_inefficiency_pe(),
        );
    }

    fn close(&mut self, _crystal: Option<&dyn Crystal>) {
        let _ = self.writer.flush();
    }
}
