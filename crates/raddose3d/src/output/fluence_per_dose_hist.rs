use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

/// Histogram accumulating weighted observations in logarithmically-spaced bins.
///
/// Port of Java `Histogram`.
struct Histogram {
    min_value: f64,
    max_value: f64,
    bucket_step: f64,
    bucket_count: usize,
    values: Vec<f64>,
    observation_weight_sum: f64,
}

impl Histogram {
    fn new(range_min: f64, range_max: f64, buckets: usize) -> Self {
        let bucket_step = (range_max - range_min) / buckets as f64;
        let bucket_count = buckets + 2;
        Histogram {
            min_value: range_min,
            max_value: range_max,
            bucket_step,
            bucket_count,
            values: vec![0.0; bucket_count],
            observation_weight_sum: 0.0,
        }
    }

    fn add_value(&mut self, position: f64, weight: f64) {
        let bucket = if position < self.min_value {
            0
        } else if position >= self.max_value {
            self.bucket_count - 1
        } else {
            1 + ((position - self.min_value) / self.bucket_step) as usize
        };
        self.values[bucket] += weight;
        self.observation_weight_sum += weight;
    }

    fn reset(&mut self) {
        for v in &mut self.values {
            *v = 0.0;
        }
        self.observation_weight_sum = 0.0;
    }

    /// Bin boundary values: -∞, min, min+step, …, max.
    fn histogram_breaks(&self) -> Vec<f64> {
        let mut breaks = vec![f64::NEG_INFINITY];
        for i in 0..self.bucket_count {
            breaks.push(self.min_value + i as f64 * self.bucket_step);
        }
        breaks
    }

    fn weight_histogram(&self) -> &[f64] {
        &self.values
    }
}

/// Outputs a per-image elastic-yield-weighted dose histogram as CSV.
///
/// Port of Java `OutputFluencePerDoseHistCSV`.
#[derive(Debug)]
pub struct OutputFluencePerDoseHistCSV {
    writer: DebugWriter,
    wedge_counter: usize,
    /// DDM-weighted elastic yield sum (numerator for avg diff dose).
    diff_intensity_num: f64,
    /// Weight sum (denominator).
    diff_intensity_denom: f64,
    hist: HistogramDebug,
}

/// Wrapper to make Histogram work with #[derive(Debug)].
struct HistogramDebug(Histogram);
impl std::fmt::Debug for HistogramDebug {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Histogram")
    }
}

impl OutputFluencePerDoseHistCSV {
    const DEFAULT_BINS: usize = 199;
    const DEFAULT_MIN: f64 = 0.1;
    const DEFAULT_MAX: f64 = 100.0;

    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        Self::with_range(
            writer,
            Self::DEFAULT_MIN,
            Self::DEFAULT_MAX,
            Self::DEFAULT_BINS,
        )
    }

    pub fn with_range(
        writer: Box<dyn Write + Send + Sync>,
        hist_min: f64,
        hist_max: f64,
        bins: usize,
    ) -> Self {
        OutputFluencePerDoseHistCSV {
            writer: DebugWriter(writer),
            wedge_counter: 1,
            diff_intensity_num: 0.0,
            diff_intensity_denom: 0.0,
            hist: HistogramDebug(Histogram::new(hist_min, hist_max, bins)),
        }
    }

    fn write_header(&mut self) {
        let breaks = self.hist.0.histogram_breaks();
        let _ = write!(
            self.writer,
            "Wedge Number, Angular position, Average Diffracted Dose, < {:e}",
            breaks[1]
        );
        for b in breaks.iter().take(breaks.len() - 1).skip(1) {
            let _ = write!(self.writer, ", {:e}", b);
        }
        let _ = writeln!(self.writer, ", > {:e}", breaks[breaks.len() - 1]);
    }
}

impl super::Output for OutputFluencePerDoseHistCSV {
    fn publish_crystal(&mut self, _crystal: &dyn Crystal) {
        self.diff_intensity_num = 0.0;
        self.diff_intensity_denom = 0.0;
        self.hist.0.reset();
        self.write_header();
    }

    fn publish_beam(&mut self, _beam: &dyn Beam) {}

    fn publish_wedge(
        &mut self,
        _wedge: &Wedge,
        _summary: &ExposureSummary,
        _crystal: Option<&dyn Crystal>,
    ) {
        self.wedge_counter += 1;
        self.hist.0.reset();
        self.diff_intensity_num = 0.0;
        self.diff_intensity_denom = 0.0;
    }

    fn close(&mut self, _crystal: Option<&dyn Crystal>) {
        let _ = self.writer.flush();
    }
}

impl super::ExposeObserver for OutputFluencePerDoseHistCSV {
    fn register(&mut self, _crystal: &dyn Crystal) {}

    fn exposure_start(&mut self, _wedge_images: usize, _wedge: &Wedge, _crystal_size: [usize; 3]) {}

    #[allow(clippy::too_many_arguments)]
    fn exposure_observation(
        &mut self,
        _wedge_image: usize,
        _i: usize,
        _j: usize,
        _k: usize,
        added_dose: f64,
        total_dose: f64,
        _fluence: f64,
        rde: f64,
        _absorbed_energy: f64,
        elastic: f64,
        _angle_count: f64,
    ) {
        let half_dose = total_dose + added_dose / 2.0;
        let decay = rde;
        let weighted_elastic = elastic * decay;

        self.hist.0.add_value(half_dose, weighted_elastic);
        self.diff_intensity_num += half_dose * weighted_elastic;
        self.diff_intensity_denom += weighted_elastic;
    }

    fn image_complete(&mut self, _image: usize, ang_rad: f64, _last_angle: f64, _vox_vol: f64) {
        let avg_diff_dose = if self.diff_intensity_denom == 0.0 {
            0.0
        } else {
            self.diff_intensity_num / self.diff_intensity_denom
        };

        let _ = write!(
            self.writer,
            "{}, {:e}, {:e}",
            self.wedge_counter,
            ang_rad.to_degrees(),
            avg_diff_dose
        );

        for &val in self.hist.0.weight_histogram() {
            let _ = write!(self.writer, ", {:e}", val);
        }
        let _ = writeln!(self.writer);

        self.hist.0.reset();
        self.diff_intensity_num = 0.0;
        self.diff_intensity_denom = 0.0;
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

    fn exposure_complete(&mut self) {}
}
