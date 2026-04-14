// Many variables in this module are intermediate simulation-state trackers
// that will be connected to outputs once the full MicroED pipeline is wired up.
#![allow(dead_code, unused_variables, unused_assignments)]
/// MicroED (Micro-Crystal Electron Diffraction) simulation engine.
///
/// Ports Java MicroED.java (~4847 lines).
use std::f64::consts::PI;

use crate::beam::Beam;
use crate::coefcalc::CoefCalc;
use crate::wedge::Wedge;

/// keV to Joules conversion factor.
pub const KEV_TO_JOULES: f64 = 1.602_176_565e-16;

/// Speed of light (m/s).
pub const C: f64 = 299_792_458.0;
/// Electron rest mass (kg).
pub const M_E: f64 = 9.109_383_56e-31;

/// Number of regions for dose tracking.
pub const NUM_REGIONS: usize = 10;

/// Number of Monte Carlo electrons to simulate.
pub const NUM_MONTE_CARLO_ELECTRONS: u64 = 10_000;

/// Broad energy sweep array (keV).
const ENERGY_ARRAY: &[f64] = &[
    10.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0,
    1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
];

/// Broad thickness sweep array (nm).
const THICKNESS_ARRAY: &[f64] = &[
    10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0,
    1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
];

// ── Main struct ───────────────────────────────────────────────────────────────

/// MicroED simulation state.
///
/// Simulates electron-beam damage in micro-crystal electron diffraction
/// using a stopping-power approach with 5 nm slices through the crystal.
pub struct MicroEdSimulation {
    // Crystal geometry
    pub cryst_coord: Vec<[f64; 3]>,
    pub cryst_pix_per_um: f64,
    pub cryst_size_voxels: [usize; 3],

    /// Crystal type string (e.g. "CUBOID", "CYLINDER", "SPHERICAL").
    pub crystal_type: String,
    /// Whether the geometry has been rotated for the current wedge angle.
    pub rotated: bool,

    // Crystal dimensions
    /// X dimension in µm.
    pub x_dimension: f64,
    /// Y dimension in µm.
    pub y_dimension: f64,
    /// Z dimension in µm.
    pub z_dimension: f64,
    /// Sample thickness in nm (= z_dimension * 1000).
    pub sample_thickness: f64,
    /// Crystal surface area in Å².
    pub crystal_surface_area: f64,
    /// Crystal volume in dm³.
    pub crystal_volume: f64,

    // Dose results
    pub dose_output: f64,
    pub gos_dose_output: f64,
    pub monte_carlo_dose: f64,
    pub monte_carlo_image_dose: f64,

    // Per-slice tracking
    /// Number of slices (= sample_thickness / slice_thickness).
    pub number_slices: usize,
    /// Slice thickness in nm.
    pub slice_thickness: f64,
    /// Per-slice dose values (MGy), length = number_slices.
    pub all_doses: Vec<f64>,

    // Electron counting metrics
    pub number_elastic: f64,
    pub number_inelastic: f64,
    pub number_single_elastic: f64,
    pub number_productive: f64,
    pub num_electrons: f64,
    pub num_simulated_electrons: u64,

    // Voxel dose / charge arrays
    voxel_dose: Vec<Vec<Vec<f64>>>,
    voxel_charge: Vec<Vec<Vec<f64>>>,
    max_x: usize,
    max_y: usize,
    max_z: usize,

    // Region-based dose tracking (NUM_REGIONS concentric shells)
    region_dose: Vec<f64>,
    region_volume: Vec<f64>,

    // Monte Carlo misc counters
    monte_carlo_tot_elastic_count: f64,
    monte_carlo_single_elastic_count: f64,
    monte_carlo_fse_escape: f64,
    monte_carlo_fse_entry: f64,
    monte_carlo_fl_escape: f64,
    monte_carlo_auger_escape: f64,
    monte_carlo_auger_entry: f64,
    monte_carlo_productive: f64,
    monte_carlo_unproductive: f64,
    monte_carlo_productive_image: f64,
    monte_carlo_productive_solvent: f64,
    monte_carlo_unproductive_micro_ed: f64,
    monte_carlo_unproductive_micro_ed_solvent: f64,
    monte_carlo_runtime: f64,
    monte_carlo_charge: f64,
    monte_carlo_charge_density: f64,
    monte_carlo_electrons_exited: f64,
    monte_carlo_electrons_entered: f64,
    monte_carlo_gos_dose: f64,
    monte_carlo_gos_escape: f64,
    energy_lost_gos: f64,
    gos_image_dose: f64,
    gos_surrounding_elastic: f64,
    gos_surrounding_elastic_image: f64,
    image_entry: f64,
    image_sec_deposited: f64,

    // FSE / Auger bookkeeping
    extra_fl_escape: f64,
    extra_auger_escape: f64,
    new_monte_carlo_fse_escape: f64,
    fse_sum: f64,
    fse_count: i32,
    low_en_dose: f64,

    // Energy partition sums
    elastic_energy_tot: f64,
    displacement_energy: f64,
    tot_fse_energy: f64,
    tot_auger_energy: f64,
    tot_shell_energy: f64,
    tot_plasmon_energy: f64,
    tot_breakdown_energy: f64,

    // Other stats
    avg_dist: f64,
    avg_voxel_dose: f64,
    avg_region_dose: f64,
    num_auger: i32,
    num_fl: i32,
    num_fse_from_surr: i32,
    num_fse_from_sample: i32,
    hit: i32,

    /// Stopping power from ESTAR (Selenium lookup — not used in Rust port).
    stopping_power_estar: f64,

    /// GOS flag (always true in Java's non-MC path).
    gos: bool,

    /// Whether the scattered medium is solvent.
    scattered_solvent: bool,

    avg_w: f64,
    w_count: f64,
    avg_shell: f64,
    inel_angle_sum: f64,
    inel_angle_count: f64,
    number_not_inelastic_equ: f64,
    number_not_inelastic_ratio: f64,

    sim_number: u64,
}

impl MicroEdSimulation {
    /// Construct a new `MicroEdSimulation` from crystal geometry data.
    ///
    /// Dimensions are derived directly from the vertex bounding box:
    /// - `vertices` – mesh vertices in µm (used only for dimension calculation)
    /// - `cryst_coord`    – flat voxel centre coordinates in µm
    /// - `cryst_pix_per_um` – voxels per µm
    /// - `cryst_size_voxels` – grid dimensions \[nx, ny, nz\]
    /// - `crystal_type` – one of `"CUBOID"`, `"CYLINDER"`, `"SPHERICAL"`, etc.
    pub fn new(
        vertices: Vec<[f64; 3]>,
        cryst_coord: Vec<[f64; 3]>,
        cryst_pix_per_um: f64,
        cryst_size_voxels: [usize; 3],
        crystal_type: String,
    ) -> Self {
        let x_min_max = min_max_vertices(0, &vertices);
        let y_min_max = min_max_vertices(1, &vertices);
        let z_min_max = min_max_vertices(2, &vertices);

        // Java: XDimension = 1000 * (xMinMax[1] - xMinMax[0])  (nm)
        // We store in µm to match the Rust convention throughout.
        let x_dim_um = x_min_max[1] - x_min_max[0]; // µm
        let y_dim_um = y_min_max[1] - y_min_max[0]; // µm
        let z_dim_um = z_min_max[1] - z_min_max[0]; // µm

        // Java stores XDimension etc. in nm, so multiply by 1000.
        // In this Rust port we keep them in µm and convert to nm only where needed.
        let x_dimension = x_dim_um; // µm
        let y_dimension = y_dim_um; // µm
        let mut z_dimension = z_dim_um; // µm

        // Java: if CYLINDER, ZDimension -= 0.001 (nm → subtract 1e-6 µm)
        if crystal_type.eq_ignore_ascii_case("CYLINDER") {
            z_dimension -= 1e-6; // 0.001 nm in µm
        }

        // sample_thickness is stored in nm (as in Java).
        let sample_thickness = z_dimension * 1000.0; // nm

        // Surface area in Å² — Java: XDimension(nm) * YDimension(nm) * 1E02
        // Convert: (x_µm * 1000 nm/µm) * (y_µm * 1000 nm/µm) * 100 A²/nm² = x_µm * y_µm * 1e8
        let mut crystal_surface_area = x_dimension * y_dimension * 1.0e8; // Å²
        if crystal_type.eq_ignore_ascii_case("CYLINDER") {
            crystal_surface_area =
                PI * (x_dimension * 1000.0 / 2.0) * (y_dimension * 1000.0 / 2.0) * 1.0e2;
        }

        // Crystal volume in dm³
        // Java: crystalSurfaceArea(Å²) * sampleThickness(nm) * 10 * 1E-27
        // = (A² * nm * 10 * 1E-27) → 1 Å = 1E-1 nm, so Å³ = 1E-3 nm³; ...
        let crystal_volume = if crystal_type.eq_ignore_ascii_case("SPHERICAL") {
            // Java: (4/3)*PI*(X/2)*(Y/2)*(Z/2) * 1E-24 with dimensions in nm
            let xnm = x_dimension * 1000.0;
            let ynm = y_dimension * 1000.0;
            let znm = z_dimension * 1000.0;
            (4.0 / 3.0) * PI * (xnm / 2.0) * (ynm / 2.0) * (znm / 2.0) * 1.0e-24
        } else {
            crystal_surface_area * (sample_thickness * 10.0) * 1.0e-27
        };

        // Initialise voxel dose/charge arrays.
        let max_pixel = get_max_pixel_coordinates(&vertices, cryst_pix_per_um);
        let nx = max_pixel[0];
        let ny = max_pixel[1];
        let nz = max_pixel[2];
        let voxel_dose = vec![vec![vec![0.0_f64; nz]; ny]; nx];
        let voxel_charge = vec![vec![vec![0.0_f64; nz]; ny]; nx];

        let mut region_volume = vec![0.0_f64; NUM_REGIONS];
        populate_region_volumes(
            x_dimension * 1000.0,
            y_dimension * 1000.0,
            z_dimension * 1000.0,
            &mut region_volume,
        );

        MicroEdSimulation {
            cryst_coord,
            cryst_pix_per_um,
            cryst_size_voxels,
            crystal_type,
            rotated: false,
            x_dimension,
            y_dimension,
            z_dimension,
            sample_thickness,
            crystal_surface_area,
            crystal_volume,
            dose_output: 0.0,
            gos_dose_output: 0.0,
            monte_carlo_dose: 0.0,
            monte_carlo_image_dose: 0.0,
            number_slices: 0,
            slice_thickness: 5.0, // nm
            all_doses: Vec::new(),
            number_elastic: 0.0,
            number_inelastic: 0.0,
            number_single_elastic: 0.0,
            number_productive: 0.0,
            num_electrons: 0.0,
            num_simulated_electrons: 0,
            voxel_dose,
            voxel_charge,
            max_x: nx,
            max_y: ny,
            max_z: nz,
            region_dose: vec![0.0; NUM_REGIONS],
            region_volume,
            monte_carlo_tot_elastic_count: 0.0,
            monte_carlo_single_elastic_count: 0.0,
            monte_carlo_fse_escape: 0.0,
            monte_carlo_fse_entry: 0.0,
            monte_carlo_fl_escape: 0.0,
            monte_carlo_auger_escape: 0.0,
            monte_carlo_auger_entry: 0.0,
            monte_carlo_productive: 0.0,
            monte_carlo_unproductive: 0.0,
            monte_carlo_productive_image: 0.0,
            monte_carlo_productive_solvent: 0.0,
            monte_carlo_unproductive_micro_ed: 0.0,
            monte_carlo_unproductive_micro_ed_solvent: 0.0,
            monte_carlo_runtime: 0.0,
            monte_carlo_charge: 0.0,
            monte_carlo_charge_density: 0.0,
            monte_carlo_electrons_exited: 0.0,
            monte_carlo_electrons_entered: 0.0,
            monte_carlo_gos_dose: 0.0,
            monte_carlo_gos_escape: 0.0,
            energy_lost_gos: 0.0,
            gos_image_dose: 0.0,
            gos_surrounding_elastic: 0.0,
            gos_surrounding_elastic_image: 0.0,
            image_entry: 0.0,
            image_sec_deposited: 0.0,
            extra_fl_escape: 0.0,
            extra_auger_escape: 0.0,
            new_monte_carlo_fse_escape: 0.0,
            fse_sum: 0.0,
            fse_count: 0,
            low_en_dose: 0.0,
            elastic_energy_tot: 0.0,
            displacement_energy: 0.0,
            tot_fse_energy: 0.0,
            tot_auger_energy: 0.0,
            tot_shell_energy: 0.0,
            tot_plasmon_energy: 0.0,
            tot_breakdown_energy: 0.0,
            avg_dist: 0.0,
            avg_voxel_dose: 0.0,
            avg_region_dose: 0.0,
            num_auger: 0,
            num_fl: 0,
            num_fse_from_surr: 0,
            num_fse_from_sample: 0,
            hit: 0,
            stopping_power_estar: 0.0,
            gos: true,
            scattered_solvent: false,
            avg_w: 0.0,
            w_count: 0.0,
            avg_shell: 0.0,
            inel_angle_sum: 0.0,
            inel_angle_count: 0.0,
            number_not_inelastic_equ: 0.0,
            number_not_inelastic_ratio: 0.0,
            sim_number: 0,
        }
    }

    // ── Public entry point ────────────────────────────────────────────────────

    /// Main MicroED calculation entry point.
    ///
    /// Mirrors Java `CalculateEM(Beam, Wedge, CoefCalc)`.
    pub fn calculate_em(&mut self, beam: &dyn Beam, wedge: &Wedge, coef_calc: &mut dyn CoefCalc) {
        let csda_distance = self.get_csda_range(coef_calc, beam);

        println!("The density is: {:.2e}", coef_calc.density());

        let _wavelength = self.get_wavelength(beam);
        let _res_rough = get_resolution_rough(_wavelength);
        let _max_res = get_max_res(_wavelength);

        let _dose1 = self.em_let_way(beam, wedge, coef_calc);

        let _dose2 = self.em_equation_way(beam, wedge, coef_calc, true);

        let dose3_results = self.em_stopping_power_way(beam, wedge, coef_calc);
        let dose3 = dose3_results[0];

        println!(
            "\nThe Dose in the exposed area by stopping power: {:.4e} MGy",
            dose3
        );
        println!("The CSDA range is: {:.2} nm", csda_distance);

        // Monte Carlo is commented out in Java — stub only.
        // self.start_monte_carlo(coef_calc, beam);

        println!("\nNumber elastic events: {:.2}", self.number_elastic);
        println!(
            "Number single elastic events: {:.2}",
            self.number_single_elastic
        );
        println!("Number Inelastic events: {:.2}", self.number_inelastic);
        println!("Number productive events: {:.2}", self.number_productive);

        let information_coef = if dose3 != 0.0 {
            self.number_productive / dose3
        } else {
            0.0
        };
        println!(
            "Information Coefficient (productive/MGy): {:.2}",
            information_coef
        );

        let optimal_en = self.get_optimal_energy(beam, wedge, coef_calc);
        println!(
            "\nThe optimal accelerating voltage is: {:.0} kV",
            optimal_en
        );
        let optimal_t = self.get_optimal_thickness(beam, wedge, coef_calc);
        println!("The optimal thickness is: {:.0} nm\n", optimal_t);

        // Store primary dose result.
        self.dose_output = dose3;

        // Store slice doses (indices 1..=N in result array).
        self.all_doses = dose3_results[1..].to_vec();

        // Write summary CSV — matches Java's WriterFile("outputMicroED.CSV").
        // Format: one header + one data row, all fields as %f (6 decimal places).
        let summary_path = "outputMicroED.CSV";
        let beam_en = beam.photon_energy();
        match Self::write_summary_csv(
            summary_path,
            beam_en,
            dose3,
            self.number_elastic,
            self.number_single_elastic,
            self.number_inelastic,
            self.number_productive,
            information_coef,
            optimal_en,
            optimal_t,
        ) {
            Ok(()) => {}
            Err(e) => eprintln!("Warning: could not write {}: {}", summary_path, e),
        }

        // Write slice data CSV — matches Java's WriterFile("outputMicroED_slices.csv").
        let slice_path = "outputMicroED_slices.csv";
        match Self::write_slice_csv(slice_path, &self.all_doses, self.slice_thickness) {
            Ok(()) => println!("Slice data written to {}", slice_path),
            Err(e) => eprintln!("Warning: could not write {}: {}", slice_path, e),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn write_summary_csv(
        path: &str,
        beam_en: f64,
        dose: f64,
        elastic: f64,
        single_elastic: f64,
        inelastic: f64,
        productive: f64,
        info_coef: f64,
        best_en: f64,
        best_t: f64,
    ) -> std::io::Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(path)?;
        writeln!(
            f,
            "Beam_en, Dose, Elastic, Single_elastic, Inelastic, Productive, Info_coef, Best_en, Best_t"
        )?;
        writeln!(
            f,
            " {beam_en:.6}, {dose:.6}, {elastic:.6}, {single_elastic:.6}, {inelastic:.6}, {productive:.6}, {info_coef:.6}, {best_en:.6}, {best_t:.6}"
        )?;
        Ok(())
    }

    fn write_slice_csv(path: &str, doses: &[f64], slice_thickness: f64) -> std::io::Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(path)?;
        writeln!(f, "Slice Number,Slice Depth (nm),Slice Dose (MGy)")?;
        for (i, &sd) in doses.iter().enumerate() {
            let depth = i as f64 * slice_thickness;
            writeln!(f, "{},{:.4},{:.8e}", i + 1, depth, sd)?;
        }
        Ok(())
    }

    // ── CSDA range ────────────────────────────────────────────────────────────

    /// Compute the Continuous-Slowing-Down Approximation (CSDA) range in nm.
    ///
    /// Mirrors Java `getCSDArange(CoefCalc, Beam)`.
    pub fn get_csda_range(&self, coef_calc: &mut dyn CoefCalc, beam: &dyn Beam) -> f64 {
        let mut en = beam.photon_energy(); // keV
        let divisions = 100;
        let energy_step = en / divisions as f64;
        let mut distance = 0.0;
        while en > 0.05 {
            let stopping_power = coef_calc.stopping_power(en, false);
            if stopping_power > 0.0 {
                distance += energy_step / stopping_power;
            }
            en -= energy_step;
        }
        distance
    }

    // ── Primary dose calculation: stopping-power way ──────────────────────────

    /// Compute dose using the stopping-power slice method.
    ///
    /// Returns a `Vec<f64>` where index 0 is the total (mean) dose in MGy,
    /// and indices 1..=N are the per-slice doses.
    ///
    /// Mirrors Java `EMStoppingPowerWay(Beam, Wedge, CoefCalc)`.
    pub fn em_stopping_power_way(
        &mut self,
        beam: &dyn Beam,
        _wedge: &Wedge,
        coef_calc: &mut dyn CoefCalc,
    ) -> Vec<f64> {
        let exposure = beam.exposure(); // e/Å²
        let exposed_area = self.exposed_area(beam); // µm²
        let electron_number = exposure * (exposed_area * 1.0e8); // electrons

        let exposed_volume = exposed_area * (self.sample_thickness / 1000.0) * 1.0e-15; // dm³
        let exposed_mass = (coef_calc.density() * 1000.0 * exposed_volume) / 1000.0; // kg

        self.number_slices = (self.sample_thickness / self.slice_thickness) as usize;
        if self.number_slices == 0 {
            self.number_slices = 1;
        }

        let exposed_mass_per_slice = exposed_mass / self.number_slices as f64;

        let mut avg_energy = beam.photon_energy(); // keV
        let mut all_doses = vec![0.0_f64; self.number_slices];
        let mut total_dose = 0.0_f64;
        let mut end_loop = false;

        for slice_dose in all_doses.iter_mut().take(self.number_slices) {
            let stopping_power = coef_calc.stopping_power(avg_energy, false); // keV/nm
            let mut energy_per_el =
                stopping_power * (self.sample_thickness / self.number_slices as f64); // keV
            avg_energy -= energy_per_el;

            if avg_energy < 0.05 {
                end_loop = true;
                energy_per_el += avg_energy; // deposit remaining energy
            }

            let energy_deposited = electron_number * energy_per_el * KEV_TO_JOULES; // J
            *slice_dose = if exposed_mass_per_slice > 0.0 {
                (energy_deposited / exposed_mass_per_slice) / 1.0e6 // MGy
            } else {
                0.0
            };
            total_dose += *slice_dose;

            if end_loop {
                break;
            }
        }

        // result[0] = mean dose; result[1..] = per-slice doses
        let mut result = vec![0.0_f64; self.number_slices + 1];
        result[0] = total_dose / self.number_slices as f64;
        result[1..].copy_from_slice(&all_doses);
        result
    }

    // ── LET-way dose ──────────────────────────────────────────────────────────

    /// LET-based dose from hardcoded lookup table (100/200/300 keV values).
    ///
    /// Mirrors Java `EMLETWay(Beam, Wedge, CoefCalc)`.
    pub fn em_let_way(
        &mut self,
        beam: &dyn Beam,
        _wedge: &Wedge,
        _coef_calc: &mut dyn CoefCalc,
    ) -> f64 {
        let exposure = beam.exposure(); // e/Å²
        let beam_energy = beam.photon_energy(); // keV

        let base_dose = if (beam_energy - 100.0).abs() < f64::EPSILON {
            6.6
        } else if (beam_energy - 200.0).abs() < f64::EPSILON {
            4.5
        } else if (beam_energy - 300.0).abs() < f64::EPSILON {
            3.7
        } else {
            0.0
        };

        base_dose * exposure
    }

    // ── Equation-way dose ─────────────────────────────────────────────────────

    /// Equation-based dose calculation.
    ///
    /// Mirrors Java `EMEquationWay(Beam, Wedge, CoefCalc, boolean)`.
    pub fn em_equation_way(
        &mut self,
        beam: &dyn Beam,
        _wedge: &Wedge,
        coef_calc: &mut dyn CoefCalc,
        _use_inel_equ: bool,
    ) -> f64 {
        let exposure = beam.exposure();
        let exposed_area = self.exposed_area(beam); // µm²
        let electron_number = exposure * (exposed_area * 1.0e8);
        self.num_electrons = electron_number;

        let avg_energy = beam.photon_energy();

        self.number_slices = (self.sample_thickness / self.slice_thickness) as usize;
        if self.number_slices == 0 {
            self.number_slices = 1;
        }

        // Elastic
        let elastic_lambda = coef_calc.electron_elastic_mfpl(avg_energy, false); // nm
        println!("Elastic Lambda: {}", elastic_lambda);
        let elastic_prob = if elastic_lambda > 0.0 {
            1.0 - (-self.sample_thickness / elastic_lambda).exp()
        } else {
            0.0
        };
        self.number_elastic = elastic_prob * electron_number;

        self.number_single_elastic = if elastic_lambda > 0.0 {
            let ratio = self.sample_thickness / elastic_lambda;
            electron_number * (-ratio).exp() * ratio
        } else {
            0.0
        };

        // GOS inelastic
        let gos_inelastic_lambda = coef_calc.gos_inel(false, avg_energy); // nm
        println!("GOS Inelastic Lambda: {}", gos_inelastic_lambda);
        let inel_prob = if gos_inelastic_lambda > 0.0 {
            1.0 - (-self.sample_thickness / gos_inelastic_lambda).exp()
        } else {
            0.0
        };
        self.number_inelastic = inel_prob * electron_number;
        self.number_productive = self.number_single_elastic * (1.0 - inel_prob);

        0.0 // Java also returns 0.0 here
    }

    // ── Optimal energy / thickness ────────────────────────────────────────────

    /// Find the accelerating voltage (keV) that maximises the information
    /// coefficient (productive / dose).
    ///
    /// Mirrors Java `getOptimalEnergy(Beam, Wedge, CoefCalc)`.
    pub fn get_optimal_energy(
        &mut self,
        beam: &dyn Beam,
        _wedge: &Wedge,
        coef_calc: &mut dyn CoefCalc,
    ) -> f64 {
        let exposure = beam.exposure();
        let exposed_area = self.exposed_area(beam);
        let electron_number = exposure * (exposed_area * 1.0e8);
        self.num_electrons = electron_number;

        let exposed_volume = exposed_area * (self.sample_thickness / 1000.0) * 1.0e-15;
        let exposed_mass = (coef_calc.density() * 1000.0 * exposed_volume) / 1000.0;

        let mut info_coefs = vec![0.0_f64; ENERGY_ARRAY.len()];
        let mut exit_i = 0usize;

        for (i, &this_energy) in ENERGY_ARRAY.iter().enumerate() {
            exit_i = i;
            let ic = self.get_info_coef(
                coef_calc,
                this_energy,
                electron_number,
                exposed_mass,
                self.sample_thickness,
            );
            info_coefs[i] = ic;
            if i > 1 && ic < info_coefs[i - 1] {
                break;
            }
        }

        // Fine-grain search in the neighbourhood of exit_i
        let step = 10.0_f64;
        let start_en = if exit_i >= 2 {
            ENERGY_ARRAY[exit_i - 2]
        } else {
            ENERGY_ARRAY[0]
        };
        let end_en = ENERGY_ARRAY[exit_i];
        let length_loop = ((end_en - start_en) / step) as usize + 1;
        let mut prev_info_coef = 0.0_f64;
        let mut this_energy = start_en;

        for i in 0..length_loop {
            this_energy = start_en + i as f64 * step;
            let ic = self.get_info_coef(
                coef_calc,
                this_energy,
                electron_number,
                exposed_mass,
                self.sample_thickness,
            );
            if i > 0 && ic < prev_info_coef {
                break;
            }
            prev_info_coef = ic;
        }

        this_energy - step
    }

    /// Find the sample thickness (nm) that maximises the information coefficient.
    ///
    /// Mirrors Java `getOptimalThickness(Beam, Wedge, CoefCalc)`.
    pub fn get_optimal_thickness(
        &mut self,
        beam: &dyn Beam,
        _wedge: &Wedge,
        coef_calc: &mut dyn CoefCalc,
    ) -> f64 {
        let avg_energy = beam.photon_energy();
        let exposure = beam.exposure();
        let exposed_area = self.exposed_area(beam);
        let electron_number = exposure * (exposed_area * 1.0e8);
        self.num_electrons = electron_number;

        let mut info_coefs = vec![0.0_f64; THICKNESS_ARRAY.len()];
        let mut exit_i = 0usize;

        for (i, &this_thickness) in THICKNESS_ARRAY.iter().enumerate() {
            exit_i = i;
            let exposed_volume = exposed_area * (this_thickness / 1000.0) * 1.0e-15;
            let exposed_mass = (coef_calc.density() * 1000.0 * exposed_volume) / 1000.0;
            let ic = self.get_info_coef(
                coef_calc,
                avg_energy,
                electron_number,
                exposed_mass,
                this_thickness,
            );
            info_coefs[i] = ic;
            if i > 1 && ic < info_coefs[i - 1] {
                break;
            }
        }

        let step = 5.0_f64;
        let start_t = if exit_i >= 2 {
            THICKNESS_ARRAY[exit_i - 2]
        } else {
            THICKNESS_ARRAY[0]
        };
        let end_t = THICKNESS_ARRAY[exit_i];
        let length_loop = ((end_t - start_t) / step) as usize + 1;
        let mut prev_info_coef = 0.0_f64;
        let mut this_thickness = start_t;

        for i in 0..length_loop {
            this_thickness = start_t + i as f64 * step;
            let exposed_volume = exposed_area * (this_thickness / 1000.0) * 1.0e-15;
            let exposed_mass = (coef_calc.density() * 1000.0 * exposed_volume) / 1000.0;
            let ic = self.get_info_coef(
                coef_calc,
                avg_energy,
                electron_number,
                exposed_mass,
                this_thickness,
            );
            if i > 0 && ic < prev_info_coef {
                break;
            }
            prev_info_coef = ic;
        }

        this_thickness - step
    }

    // ── Exposed dimension helpers ─────────────────────────────────────────────

    /// Returns the exposed X extent: min(beam_x, crystal_x) in µm.
    ///
    /// Mirrors Java `getExposedX(Beam)`.
    pub fn get_exposed_x(&self, beam: &dyn Beam) -> f64 {
        let beam_x = beam.beam_size_x(); // µm
        let crystal_x = self.x_dimension; // µm
        if crystal_x > beam_x {
            beam_x
        } else {
            crystal_x
        }
    }

    /// Returns the exposed Y extent: min(beam_y, crystal_y) in µm.
    ///
    /// Mirrors Java `getExposedY(Beam)`.
    pub fn get_exposed_y(&self, beam: &dyn Beam) -> f64 {
        let beam_y = beam.beam_size_y(); // µm
        let crystal_y = self.y_dimension; // µm
        if crystal_y > beam_y {
            beam_y
        } else {
            crystal_y
        }
    }

    // ── Electron wavelength / resolution ──────────────────────────────────────

    /// Relativistic de Broglie wavelength in Å.
    ///
    /// Mirrors Java `getWavelength(Beam)`.
    pub fn get_wavelength(&self, beam: &dyn Beam) -> f64 {
        let h = 6.626_070_040e-34_f64;
        let c = 299_792_458.0_f64;
        let m0 = 9.109_383_356e-31_f64;
        let v0 = beam.photon_energy() * KEV_TO_JOULES;
        let lambda = (h * c) / (v0 * v0 + 2.0 * v0 * m0 * c * c).sqrt(); // m
        lambda * 1.0e10 // Å
    }

    // ── Monte Carlo stub ──────────────────────────────────────────────────────

    /// Monte Carlo simulation stub.
    ///
    /// The full Monte Carlo path is almost entirely commented out in the Java
    /// source; this stub mirrors the structure without the implementation.
    pub fn start_monte_carlo(&mut self, _coef_calc: &dyn CoefCalc, _beam: &dyn Beam) {
        // Monte Carlo simulation is commented out in Java MicroED.java.
        // Stub only.
    }

    /// Process Monte Carlo dose results (stub).
    ///
    /// Mirrors Java `processMonteCarloDose(Beam, CoefCalc)`.
    pub fn process_monte_carlo_dose(
        &self,
        _beam: &dyn Beam,
        _coef_calc: &dyn CoefCalc,
    ) -> [f64; 2] {
        [self.monte_carlo_dose, self.monte_carlo_image_dose]
    }

    /// Convert voxel energy (keV per simulated electron) to dose (MGy).
    ///
    /// Mirrors Java `convertVoxEnergyToDose(...)`.
    pub fn convert_vox_energy_to_dose(
        &self,
        energy: f64,
        beam: &dyn Beam,
        coef_calc: &dyn CoefCalc,
    ) -> f64 {
        if self.num_simulated_electrons == 0 {
            return 0.0;
        }
        let electron_number = beam.exposure() * (beam.beam_area() * 1.0e8);
        let total_j =
            (energy * (electron_number / self.num_simulated_electrons as f64)) * KEV_TO_JOULES;
        let voxel_side = 1.0 / self.cryst_pix_per_um / 1.0e4; // cm
        let voxel_volume = voxel_side.powi(3); // cm³
        let voxel_mass = (coef_calc.density() * voxel_volume) / 1000.0; // kg
        if voxel_mass > 0.0 {
            (total_j / voxel_mass) / 1.0e6
        } else {
            0.0
        }
    }

    /// Convert region energy (keV per simulated electron) to dose (MGy).
    ///
    /// Mirrors Java `convertRegionEnergyToDose(...)`.
    pub fn convert_region_energy_to_dose(
        &self,
        energy: f64,
        index: usize,
        beam: &dyn Beam,
        coef_calc: &dyn CoefCalc,
    ) -> f64 {
        if self.num_simulated_electrons == 0 {
            return 0.0;
        }
        let electron_number = beam.exposure() * (beam.beam_area() * 1.0e8);
        let total_j =
            (energy * (electron_number / self.num_simulated_electrons as f64)) * KEV_TO_JOULES;
        let volume = self.region_volume[index] / 1.0e21; // cm³
        let region_mass = (coef_calc.density() * volume) / 1000.0; // kg
        if region_mass > 0.0 {
            (total_j / region_mass) / 1.0e6
        } else {
            0.0
        }
    }

    // ── Utility methods ───────────────────────────────────────────────────────

    /// Compute the information coefficient (productive events / dose).
    ///
    /// Mirrors Java `getInfoCoef(CoefCalc, double, double, double, double)`.
    fn get_info_coef(
        &self,
        coef_calc: &mut dyn CoefCalc,
        test_energy: f64,
        electron_number: f64,
        exposed_mass: f64,
        this_thickness: f64,
    ) -> f64 {
        let stopping_power = coef_calc.stopping_power(test_energy, false);
        let energy_per_el = stopping_power * this_thickness;
        let energy_deposited = electron_number * energy_per_el * KEV_TO_JOULES;
        let dose = if exposed_mass > 0.0 {
            (energy_deposited / exposed_mass) / 1.0e6
        } else {
            0.0
        };

        let elastic_lambda = coef_calc.electron_elastic_mfpl(test_energy, false);
        let inel_lambda = coef_calc.gos_inel(false, test_energy);

        let inel_prob = if inel_lambda > 0.0 {
            1.0 - (-this_thickness / inel_lambda).exp()
        } else {
            0.0
        };
        let number_single_elastic = if elastic_lambda > 0.0 {
            let ratio = this_thickness / elastic_lambda;
            electron_number * (-ratio).exp() * ratio
        } else {
            0.0
        };
        let productive_el = number_single_elastic * (1.0 - inel_prob);

        if dose > 0.0 {
            productive_el / dose
        } else {
            0.0
        }
    }

    /// Compute exposed area in µm² (rectangular or circular beam).
    fn exposed_area(&self, beam: &dyn Beam) -> f64 {
        let ex = self.get_exposed_x(beam);
        let ey = self.get_exposed_y(beam);
        if beam.is_circular() {
            PI * (ex / 2.0) * (ey / 2.0)
        } else {
            ex * ey
        }
    }

    // ── GOS angle / energy helpers (used in full MC — stubs) ─────────────────

    /// Wk to Wak conversion.
    pub fn wk_to_wak(&self, e: f64, wk: f64, uk: f64) -> f64 {
        if e * 1000.0 > 3.0 * wk - 2.0 * uk {
            wk
        } else {
            (e * 1000.0 + 2.0 * uk) / 3.0
        }
    }

    /// Q_ak parameter.
    pub fn get_qak(&self, e: f64, wk: f64, uk: f64) -> f64 {
        if e * 1000.0 > 3.0 * wk - 2.0 * uk {
            uk
        } else {
            uk * (e * 1000.0 / (3.0 * wk - 2.0 * uk))
        }
    }

    /// Distant energy-loss sampling.
    pub fn get_energy_loss_distant(&self, wdis: f64, uk: f64) -> f64 {
        let rnd: f64 = rand_unit();
        wdis - (rnd * (wdis - uk).powi(2)).sqrt()
    }

    /// Close-collision parameter `a`.
    pub fn get_close_a(&self, e: f64) -> f64 {
        let csq = C * C;
        let vo = e * KEV_TO_JOULES;
        let denom = vo + M_E * csq;
        (vo / denom).powi(2)
    }

    /// GOS primary scattering angle (theta) for a long-range interaction.
    pub fn get_gos_primary_theta_long(
        &self,
        e_kev: f64,
        q: f64,
        wa_kev: f64,
        _previous_theta: f64,
    ) -> f64 {
        let csq = C * C;
        let e = e_kev * KEV_TO_JOULES;
        let wak = (wa_kev / 1000.0) * KEV_TO_JOULES;
        let num = e * (e + 2.0 * M_E * csq) + (e - wak) * (e - wak + 2.0 * M_E * csq)
            - q * (q + 2.0 * M_E * csq);
        let den =
            2.0 * (e * (e + 2.0 * M_E * csq) * (e - wak) * (e - wak + 2.0 * M_E * csq)).sqrt();
        if den.abs() > 0.0 {
            (num / den).clamp(-1.0, 1.0).acos()
        } else {
            0.0
        }
    }

    /// Random k for GOS sampling.
    pub fn get_random_k(&self, e: f64, qk: f64) -> f64 {
        let wcc = 0.0_f64; // Wcc constant
        let kc = wcc.max(qk) / (e * 1000.0 + qk);
        let a = self.get_close_a(e);
        let rnd: f64 = rand_unit();
        let zeta = rnd * (1.0 + 5.0 * a * kc / 2.0);
        if zeta < 1.0 {
            kc / (1.0 - zeta * (1.0 - 2.0 * kc))
        } else {
            kc + (zeta - 1.0) * (1.0 - 2.0 * kc) / (5.0 * a * kc)
        }
    }

    // ── FSE cross-section (used in full MC) ───────────────────────────────────

    /// Gryzinski FSE cross-section (nm²/sr per electron).
    pub fn get_fse_x_section(&self, energy: f64) -> f64 {
        // Simplified Gryzinski-style estimate; full version uses tables in MC.
        let e_ev = energy * 1000.0;
        let phi = (e_ev / 2.0) / (1.0 + e_ev / 2.0);
        phi * 6.56e-14 / (e_ev * e_ev)
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Returns min and max vertex coordinate along `dimension` (0=x, 1=y, 2=z).
///
/// Mirrors Java `minMaxVertices(int, double[][])`.
pub fn min_max_vertices(dimension: usize, vertices: &[[f64; 3]]) -> [f64; 2] {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for v in vertices {
        if v[dimension] < min {
            min = v[dimension];
        }
        if v[dimension] > max {
            max = v[dimension];
        }
    }
    [min, max]
}

/// Compute maximum pixel coordinate extents from vertex bounding box.
///
/// Mirrors Java `getMaxPixelCoordinates()`.
fn get_max_pixel_coordinates(vertices: &[[f64; 3]], cryst_pix_per_um: f64) -> [usize; 3] {
    let xmm = min_max_vertices(0, vertices);
    let ymm = min_max_vertices(1, vertices);
    let zmm = min_max_vertices(2, vertices);

    let nx = ((xmm[1] - xmm[0]) * cryst_pix_per_um).round() as usize + 1;
    let ny = ((ymm[1] - ymm[0]) * cryst_pix_per_um).round() as usize + 1;
    let nz = ((zmm[1] - zmm[0]) * cryst_pix_per_um).round() as usize + 1;
    [nx.max(1), ny.max(1), nz.max(1)]
}

/// Populate concentric-shell region volumes.
///
/// Mirrors Java `populateRegionVolumes()`.
fn populate_region_volumes(xnm: f64, ynm: f64, znm: f64, region_volume: &mut [f64]) {
    let total_volume = xnm * ynm * znm;
    let mut sum_volume = 0.0_f64;
    for (i, rv) in region_volume.iter_mut().enumerate().take(NUM_REGIONS) {
        let inner_x = xnm - (i + 1) as f64 * (xnm / NUM_REGIONS as f64);
        let inner_y = ynm - (i + 1) as f64 * (ynm / NUM_REGIONS as f64);
        let inner_volume = inner_x.max(0.0) * inner_y.max(0.0) * znm;
        *rv = total_volume - (inner_volume + sum_volume);
        sum_volume += *rv;
    }
}

/// Rough maximum resolution.
///
/// Mirrors Java `getResolutionRough(double wavelength)`.
pub fn get_resolution_rough(wavelength: f64) -> f64 {
    let a = 0.01_f64; // radians
    let n = 1.0_f64;
    wavelength / (2.0 * n * a.sin())
}

/// Maximum resolution from aberration limit.
///
/// Mirrors Java `getMaxRes(double wavelength)`.
pub fn get_max_res(wavelength: f64) -> f64 {
    let cs = 1.0e7_f64; // Å
    (cs * wavelength.powi(3) / 6.0).powf(0.25)
}

/// Very basic PRNG shim — returns a value in [0, 1).
///
/// In the Java Monte Carlo, `Math.random()` is used extensively.
/// The full MC is not implemented; this stub is used by the few non-MC
/// helpers that perform random sampling.
fn rand_unit() -> f64 {
    // Simple LCG — not for scientific use; replace with proper RNG if MC is enabled.
    use std::sync::atomic::{AtomicU64, Ordering};
    static STATE: AtomicU64 = AtomicU64::new(0x853c_49e6_748f_ea9b);
    let mut x = STATE.load(Ordering::Relaxed);
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    STATE.store(x, Ordering::Relaxed);
    (x as f64) / (u64::MAX as f64)
}
