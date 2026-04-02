use crate::beam::{self, Beam};
use crate::crystal::{self, Crystal};
use crate::output::Output;
use crate::parser::config::Config;
use crate::wedge::Wedge;

/// Experiment is the central coordinating class.
/// It receives Crystal, Beam, and Wedge objects, initiates exposure,
/// and notifies subscribed Output modules.
#[derive(Debug)]
pub struct Experiment {
    current_crystal: Option<Box<dyn Crystal>>,
    current_beam: Option<Box<dyn Beam>>,
    observers: Vec<Box<dyn Output>>,
}

impl Default for Experiment {
    fn default() -> Self {
        Self::new()
    }
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
                o.publish_wedge(wedge, summary, Some(crystal.as_ref()));
            }
        }
    }

    /// Close the experiment.
    pub fn close(&mut self) {
        let crystal_ref: Option<&dyn Crystal> = self.current_crystal.as_deref();
        for o in &mut self.observers {
            o.close(crystal_ref);
        }
        self.observers.clear();
        self.current_beam = None;
        self.current_crystal = None;
    }

    /// Access the crystal (if set).
    pub fn crystal(&self) -> Option<&dyn Crystal> {
        self.current_crystal.as_deref()
    }
}

impl Experiment {
    /// Run a full simulation from a parsed Config.
    pub fn run_from_config(config: &Config) -> Result<Self, String> {
        let mut experiment = Experiment::new();

        // Create default output observers
        let summary_text = crate::output::OutputSummaryText::new(Box::new(std::io::stdout()));
        experiment.add_observer(Box::new(summary_text));

        // Process items in declaration order so each beam applies to the
        // wedges that follow it (not all wedges using the last beam).
        use crate::parser::config::ConfigItem;
        for item in &config.items {
            match item {
                ConfigItem::Crystal(c) => {
                    let crystal = crystal::create_crystal(c)?;
                    experiment.set_crystal(crystal);
                }
                ConfigItem::Beam(b) => {
                    let beam = beam::create_beam(b)?;
                    experiment.set_beam(beam);
                }
                ConfigItem::Wedge(w) => {
                    let wedge = Wedge::from_config(w);
                    experiment.expose_wedge(&wedge);
                }
            }
        }

        experiment.close();
        Ok(experiment)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};

    /// Event types recorded by the test observer.
    #[derive(Debug, Clone, PartialEq)]
    enum Event {
        Crystal,
        Beam,
        Wedge,
        Close,
    }

    /// A mock Output observer that records events in order.
    #[derive(Debug)]
    struct TestObserver {
        events: Arc<Mutex<Vec<Event>>>,
    }

    impl TestObserver {
        fn new(events: Arc<Mutex<Vec<Event>>>) -> Self {
            TestObserver { events }
        }
    }

    impl crate::output::Output for TestObserver {
        fn publish_crystal(&mut self, _crystal: &dyn Crystal) {
            self.events.lock().unwrap().push(Event::Crystal);
        }
        fn publish_beam(&mut self, _beam: &dyn Beam) {
            self.events.lock().unwrap().push(Event::Beam);
        }
        fn publish_wedge(
            &mut self,
            _wedge: &Wedge,
            _summary: &crate::output::ExposureSummary,
            _crystal: Option<&dyn Crystal>,
        ) {
            self.events.lock().unwrap().push(Event::Wedge);
        }
        fn close(&mut self, _crystal: Option<&dyn Crystal>) {
            self.events.lock().unwrap().push(Event::Close);
        }
    }

    /// Helper to build a minimal crystal config for testing.
    fn make_test_crystal() -> Box<dyn Crystal> {
        use crate::parser::config::{CoefCalcType, CrystalConfig, CrystalType};
        let config = CrystalConfig {
            crystal_type: Some(CrystalType::Cuboid),
            coefcalc: Some(CoefCalcType::Average),
            dim_x: Some(10.0),
            dim_y: Some(10.0),
            dim_z: Some(10.0),
            pixels_per_micron: Some(0.5),
            angle_p: Some(0.0),
            angle_l: Some(0.0),
            cell_a: Some(100.0),
            cell_b: Some(100.0),
            cell_c: Some(100.0),
            ..Default::default()
        };
        crystal::create_crystal(&config).unwrap()
    }

    fn make_test_beam() -> Box<dyn Beam> {
        use crate::parser::config::{BeamConfig, BeamType, Collimation};
        let config = BeamConfig {
            beam_type: Some(BeamType::Tophat),
            energy: Some(12.0),
            flux: Some(1e10),
            collimation: Some(Collimation::Rectangular { h: 10.0, v: 10.0 }),
            ..Default::default()
        };
        beam::create_beam(&config).unwrap()
    }

    #[test]
    fn experiment_close_notifies_observer() {
        let events = Arc::new(Mutex::new(Vec::new()));
        let mut exp = Experiment::new();
        exp.add_observer(Box::new(TestObserver::new(events.clone())));

        let crystal = make_test_crystal();
        exp.set_crystal(crystal);
        // close without beam or wedge
        exp.close();

        let log = events.lock().unwrap();
        assert_eq!(log[0], Event::Crystal);
        assert_eq!(*log.last().unwrap(), Event::Close);
        // No beam or wedge events since we didn't set those
        assert!(!log.contains(&Event::Beam));
        assert!(!log.contains(&Event::Wedge));
    }

    #[test]
    fn experiment_multiple_observers_receive_events() {
        let events1 = Arc::new(Mutex::new(Vec::new()));
        let events2 = Arc::new(Mutex::new(Vec::new()));
        let events3 = Arc::new(Mutex::new(Vec::new()));

        let mut exp = Experiment::new();

        // Observer 1 subscribes first
        exp.add_observer(Box::new(TestObserver::new(events1.clone())));
        let crystal = make_test_crystal();
        exp.set_crystal(crystal);

        // Observer 2 subscribes after crystal
        exp.add_observer(Box::new(TestObserver::new(events2.clone())));

        // Set beam
        let beam = make_test_beam();
        exp.set_beam(beam);

        // Observer 3 subscribes late
        exp.add_observer(Box::new(TestObserver::new(events3.clone())));

        exp.close();

        // Observer 1 should have seen Crystal, Beam, Close
        let log1 = events1.lock().unwrap();
        assert!(log1.contains(&Event::Crystal));
        assert!(log1.contains(&Event::Beam));
        assert!(log1.contains(&Event::Close));

        // Observer 2 should not have seen Crystal (subscribed after)
        let log2 = events2.lock().unwrap();
        assert!(!log2.contains(&Event::Crystal));
        assert!(log2.contains(&Event::Beam));
        assert!(log2.contains(&Event::Close));

        // Observer 3 should only have seen Close
        let log3 = events3.lock().unwrap();
        assert!(!log3.contains(&Event::Crystal));
        assert!(!log3.contains(&Event::Beam));
        assert!(log3.contains(&Event::Close));
    }
}
