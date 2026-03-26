pub mod beam;
pub mod coefcalc;
pub mod container;
pub mod crystal;
pub mod ddm;
pub mod element;
pub mod energy_distribution;
pub mod experiment;
pub mod output;
pub mod residue;
pub mod simulation;
pub mod wedge;
pub mod writer;

pub use raddose3d_parser as parser;

/// Physical constants used throughout RADDOSE-3D.
pub mod constants {
    /// Elementary charge in coulombs.
    pub const ELEMENTARY_CHARGE: f64 = 1.602_176_565e-19;
    /// Conversion factor from keV to Joules.
    pub const KEV_TO_JOULES: f64 = 1000.0 * ELEMENTARY_CHARGE;
    /// Atomic mass unit in grams.
    pub const ATOMIC_MASS_UNIT: f64 = 1.66e-24;
    /// Avogadro's number.
    pub const AVOGADRO_NUM: f64 = 6.022e23;
    /// Water concentration in mM.
    pub const WATER_CONCENTRATION: f64 = 51666.0;
    /// Angstroms cubed to mL conversion.
    pub const ANGSTROMS_TO_ML: f64 = 1e-24;
    /// Number of static exposure steps when crystal is exposed without rotation.
    pub const STATIC_EXPOSURE: usize = 100;
    /// Conversion factor from Gy to MGy.
    pub const GY_TO_MGY: f64 = 1e-6;
    /// Unit conversion to get voxel mass in kg (1e-15 = 1e-18/1e-3).
    pub const UNIT_CONVERSION: f64 = 1e-15;
}
