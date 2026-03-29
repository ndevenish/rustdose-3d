pub mod dwds;
pub mod exposure_summary;
pub mod final_dose_state;
pub mod fluence_per_dose_hist;
pub mod html_report;
pub mod progress;
pub mod rde_csv;
pub mod summary;
pub mod voxel_data;

pub use dwds::OutputDWDs;
pub use exposure_summary::ExposureSummary;
pub use final_dose_state::{OutputFinalDoseStateCSV, OutputFinalDoseStateR};
pub use fluence_per_dose_hist::OutputFluencePerDoseHistCSV;
pub use html_report::OutputDoseStateHTML;
pub use progress::{OutputProgressEstimate, OutputProgressIndicator};
pub use rde_csv::OutputRDECSV;
pub use summary::{OutputSummaryCSV, OutputSummaryText};
pub use voxel_data::{OutputVoxelDose, OutputVoxelFluences};

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::wedge::Wedge;

/// Newtype wrapper around `Write` that implements `Debug`.
/// Used by all output modules that hold a boxed writer.
pub(crate) struct DebugWriter(pub(crate) Box<dyn std::io::Write + Send + Sync>);

impl std::fmt::Debug for DebugWriter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Writer")
    }
}

impl std::ops::Deref for DebugWriter {
    type Target = dyn std::io::Write;
    fn deref(&self) -> &Self::Target {
        &*self.0
    }
}

impl std::ops::DerefMut for DebugWriter {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut *self.0
    }
}

/// Output observer trait: receives notifications about experiment events.
pub trait Output: std::fmt::Debug + Send + Sync {
    /// A new crystal has been set up.
    fn publish_crystal(&mut self, crystal: &dyn Crystal);
    /// A new beam has been set up.
    fn publish_beam(&mut self, beam: &dyn Beam);
    /// A wedge exposure has completed.
    /// `crystal` is the crystal after exposure (coefficients already updated).
    fn publish_wedge(
        &mut self,
        wedge: &Wedge,
        summary: &ExposureSummary,
        crystal: Option<&dyn Crystal>,
    );
    /// Clean up and close.
    /// `crystal` is the last seen crystal (still alive at this point), if any.
    /// Default implementation ignores the crystal and just flushes.
    fn close(&mut self, crystal: Option<&dyn Crystal>) {
        let _ = crystal;
    }
}

/// Extended output that also receives warnings and references.
pub trait ExperimentNotices: Output {
    fn raise_warning(&mut self, warning: &str);
    fn add_reference(&mut self, reference: &str);
}

/// ExposeObserver: observes individual voxel exposures within a crystal.
pub trait ExposeObserver: std::fmt::Debug + Send + Sync {
    /// Called when registering onto a crystal.
    fn register(&mut self, crystal: &dyn Crystal);

    /// Called before exposure begins.
    fn exposure_start(&mut self, wedge_images: usize, wedge: &Wedge, crystal_size: [usize; 3]);

    /// Single voxel exposure event.
    #[allow(clippy::too_many_arguments)]
    fn exposure_observation(
        &mut self,
        wedge_image: usize,
        i: usize,
        j: usize,
        k: usize,
        added_dose: f64,
        total_dose: f64,
        fluence: f64,
        rde: f64,
        absorbed_energy: f64,
        elastic: f64,
        angle_count: f64,
    );

    /// Completion of exposure at a particular angle.
    fn image_complete(&mut self, image: usize, ang_rad: f64, last_angle: f64, vox_vol: f64);

    /// Voxel summary after wedge exposure.
    fn summary_observation(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        total_dose: f64,
        voxel_mass_kg: f64,
    );

    /// Called at end of wedge exposure.
    fn exposure_complete(&mut self);
}
