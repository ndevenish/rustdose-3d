use serde::{Deserialize, Serialize};

/// Top-level configuration parsed from a RADDOSE-3D input file.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Config {
    pub crystals: Vec<CrystalConfig>,
    pub beams: Vec<BeamConfig>,
    pub wedges: Vec<WedgeConfig>,
}

/// Absorption coefficient calculation method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CoefCalcType {
    Average,    // 1: DUMMY / AVERAGE
    Default,    // 2: DEFAULT / RDJAVA / RD3D
    RdFortran,  // 3: RDV2 / RDV3
    Pdb,        // 4: EXP
    Saxs,       // 5: SAXS
    Sequence,   // 6: SEQUENCE
    SaxsSeq,    // 7: SAXSSEQ
    SmallMole,  // 8: SMALLMOLE
    Cif,        // 9: EXPSM / CIF
}

/// Diffraction decay model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DdmType {
    Simple,
    Linear,
    Leal,
    Bfactor,
}

/// Container material type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ContainerMaterialType {
    None,
    Mixture,
    Elemental,
}

/// Collimation shape.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Collimation {
    Rectangular { h: f64, v: f64 },
    Circular { h: f64, v: f64 },
    Horizontal { h: f64 },
    Vertical { v: f64 },
    Unspecified,
}

/// Element-count pair (e.g., "S 2700" or "Si 1 O 2").
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElementCount {
    pub symbol: String,
    pub count: f64,
}

/// Crystal block configuration.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CrystalConfig {
    pub crystal_type: Option<String>,
    pub coefcalc: Option<CoefCalcType>,
    pub ddm: Option<DdmType>,
    pub container_material: Option<ContainerMaterialType>,

    // Dimensions
    pub dim_x: Option<f64>,
    pub dim_y: Option<f64>,
    pub dim_z: Option<f64>,
    pub pixels_per_micron: Option<f64>,

    // Angles
    pub angle_p: Option<f64>,
    pub angle_l: Option<f64>,

    // Unit cell
    pub cell_a: Option<f64>,
    pub cell_b: Option<f64>,
    pub cell_c: Option<f64>,
    pub cell_alpha: Option<f64>,
    pub cell_beta: Option<f64>,
    pub cell_gamma: Option<f64>,

    // Composition
    pub num_monomers: Option<i32>,
    pub num_residues: Option<i32>,
    pub num_rna: Option<i32>,
    pub num_dna: Option<i32>,
    pub num_carb: Option<i32>,
    pub heavy_protein_atoms: Vec<ElementCount>,
    pub small_mole_atoms: Vec<ElementCount>,
    pub heavy_solution_conc: Vec<ElementCount>,
    pub solvent_fraction: Option<f64>,

    // PDB / CIF / Sequence
    pub pdb: Option<String>,
    pub cif: Option<String>,
    pub seq_file: Option<String>,
    pub protein_conc: Option<f64>,

    // Wireframe / model
    pub wireframe_type: Option<String>,
    pub model_file: Option<String>,

    // Escape calculations
    pub calculate_pe_escape: Option<String>,
    pub calculate_fl_escape: Option<String>,
    pub fl_resolution: Option<i32>,
    pub pe_resolution: Option<i32>,

    // Container
    pub container_thickness: Option<f64>,
    pub container_density: Option<f64>,
    pub container_mixture: Option<String>,
    pub container_elements: Vec<ElementCount>,

    // DDM parameters
    pub gamma_param: Option<f64>,
    pub b0_param: Option<f64>,
    pub beta_param: Option<f64>,

    // Surrounding / oil
    pub surrounding_heavy_conc: Vec<ElementCount>,
    pub surrounding_thickness: Option<[f64; 3]>,
    pub oil_based: Option<String>,
    pub calc_surrounding: Option<String>,
    pub oil_elements: Vec<ElementCount>,
    pub oil_density: Option<f64>,

    // Goniometer
    pub goniometer_axis: Option<f64>,
    pub polarisation_direction: Option<f64>,

    // Simulation
    pub sim_electrons: Option<i64>,
    pub program: Option<String>,
    pub runs: Option<i32>,
}

/// Beam block configuration.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BeamConfig {
    pub beam_type: Option<String>,
    pub flux: Option<f64>,
    pub fwhm_x: Option<f64>,
    pub fwhm_y: Option<f64>,
    pub energy: Option<f64>,
    pub collimation: Option<Collimation>,
    pub file: Option<String>,
    pub pixel_size_x: Option<f64>,
    pub pixel_size_y: Option<f64>,
    pub exposure: Option<f64>,
    pub semi_angle: Option<f64>,
    pub aperture_radius: Option<f64>,
    pub image_x: Option<f64>,
    pub image_y: Option<f64>,
    pub pulse_energy: Option<f64>,
    pub energy_fwhm: Option<f64>,
}

/// Wedge block configuration.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WedgeConfig {
    pub start_ang: f64,
    pub end_ang: f64,
    pub exposure_time: Option<f64>,
    pub angular_resolution: Option<f64>,
    pub start_offset_x: Option<f64>,
    pub start_offset_y: Option<f64>,
    pub start_offset_z: Option<f64>,
    pub translate_x: Option<f64>,
    pub translate_y: Option<f64>,
    pub translate_z: Option<f64>,
    pub rot_ax_beam_offset: Option<f64>,
    pub max_resolution: Option<f64>,
}
