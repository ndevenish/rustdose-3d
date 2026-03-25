pub mod exposure_summary;
pub mod summary;

pub use exposure_summary::ExposureSummary;
pub use summary::{OutputSummaryCSV, OutputSummaryText};

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::wedge::Wedge;

/// Output observer trait: receives notifications about experiment events.
pub trait Output: std::fmt::Debug + Send + Sync {
    /// A new crystal has been set up.
    fn publish_crystal(&mut self, crystal: &dyn Crystal);
    /// A new beam has been set up.
    fn publish_beam(&mut self, beam: &dyn Beam);
    /// A wedge exposure has completed.
    fn publish_wedge(&mut self, wedge: &Wedge, summary: &ExposureSummary);
    /// Clean up and close.
    fn close(&mut self);
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
