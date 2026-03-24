use crate::beam::{self, Beam};
use crate::crystal::{self, Crystal};
use crate::output::Output;
use crate::wedge::Wedge;
use crate::parser::config::Config;

/// Experiment is the central coordinating class.
/// It receives Crystal, Beam, and Wedge objects, initiates exposure,
/// and notifies subscribed Output modules.
#[derive(Debug)]
pub struct Experiment {
    current_crystal: Option<Box<dyn Crystal>>,
    current_beam: Option<Box<dyn Beam>>,
    observers: Vec<Box<dyn Output>>,
}

impl Experiment {
    pub fn new() -> Self {
        Experiment {
            current_crystal: None,
            current_beam: None,
            observers: Vec::new(),
        }
    }

    /// Add an output observer.
    pub fn add_observer(&mut self, observer: Box<dyn Output>) {
        self.observers.push(observer);
    }

    /// Set the crystal and notify observers.
    pub fn set_crystal(&mut self, crystal: Box<dyn Crystal>) {
        self.current_crystal = Some(crystal);
        if let Some(ref c) = self.current_crystal {
            for o in &mut self.observers {
                o.publish_crystal(&**c);
            }
        }
    }

    /// Set the beam and notify observers.
    pub fn set_beam(&mut self, beam: Box<dyn Beam>) {
        self.current_beam = Some(beam);
        if let Some(ref b) = self.current_beam {
            for o in &mut self.observers {
                o.publish_beam(&**b);
            }
        }
    }

    /// Expose a wedge and notify observers.
    pub fn expose_wedge(&mut self, wedge: &Wedge) {
        if let (Some(ref mut crystal), Some(ref mut beam)) =
            (&mut self.current_crystal, &mut self.current_beam)
        {
            crystal.expose(&mut **beam, wedge);
            let summary = crystal.exposure_summary();
            for o in &mut self.observers {
                o.publish_wedge(wedge, summary);
            }
        }
    }

    /// Close the experiment.
    pub fn close(&mut self) {
        for o in &mut self.observers {
            o.close();
        }
        self.observers.clear();
        self.current_beam = None;
        self.current_crystal = None;
    }

    /// Access the crystal (if set).
    pub fn crystal(&self) -> Option<&dyn Crystal> {
        self.current_crystal.as_deref()
    }

    /// Run a full simulation from a parsed Config.
    pub fn run_from_config(config: &Config) -> Result<Self, String> {
        let mut experiment = Experiment::new();

        // Create default output observers
        let summary_text = crate::output::OutputSummaryText::new(Box::new(std::io::stdout()));
        experiment.add_observer(Box::new(summary_text));

        // Process each crystal/beam/wedge
        for crystal_config in &config.crystals {
            let crystal = crystal::create_crystal(crystal_config)?;
            experiment.set_crystal(crystal);
        }

        for beam_config in &config.beams {
            let beam = beam::create_beam(beam_config)?;
            experiment.set_beam(beam);
        }

        for wedge_config in &config.wedges {
            let wedge = Wedge::from_config(wedge_config);
            experiment.expose_wedge(&wedge);
        }

        experiment.close();
        Ok(experiment)
    }
}
