pub mod config;

use pest::Parser;
use pest_derive::Parser;

pub use config::*;

#[derive(Parser)]
#[grammar = "grammar.pest"]
struct InputParser;

/// Parse a RADDOSE-3D input file string into a Config.
pub fn parse(input: &str) -> Result<Config, ParseError> {
    let pairs = InputParser::parse(Rule::configfile, input)
        .map_err(|e| ParseError::Grammar(e.to_string()))?;

    let mut config = Config::default();

    for pair in pairs {
        if pair.as_rule() == Rule::configfile {
            for inner in pair.into_inner() {
                match inner.as_rule() {
                    Rule::crystal => {
                        let c = parse_crystal(inner)?;
                        config.items.push(ConfigItem::Crystal(Box::new(c.clone())));
                        config.crystals.push(c);
                    }
                    Rule::beam => {
                        let b = parse_beam(inner)?;
                        config.items.push(ConfigItem::Beam(Box::new(b.clone())));
                        config.beams.push(b);
                    }
                    Rule::wedge => {
                        let w = parse_wedge(inner)?;
                        config.items.push(ConfigItem::Wedge(Box::new(w.clone())));
                        config.wedges.push(w);
                    }
                    _ => {}
                }
            }
        }
    }

    Ok(config)
}

#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("Grammar error: {0}")]
    Grammar(String),
    #[error("Invalid number: {0}")]
    InvalidNumber(String),
}

type Pairs<'a> = pest::iterators::Pair<'a, Rule>;

fn parse_float(pair: &Pairs<'_>) -> Result<f64, ParseError> {
    pair.as_str()
        .parse::<f64>()
        .map_err(|e| ParseError::InvalidNumber(format!("{}: '{}'", e, pair.as_str())))
}

fn parse_int(pair: &Pairs<'_>) -> Result<i32, ParseError> {
    let f = parse_float(pair)?;
    Ok(f as i32)
}

fn parse_long(pair: &Pairs<'_>) -> Result<i64, ParseError> {
    let f = parse_float(pair)?;
    Ok(f as i64)
}

fn parse_element_counts(pair: Pairs<'_>) -> Result<Vec<ElementCount>, ParseError> {
    let mut result = Vec::new();
    let mut inner = pair.into_inner();
    while let Some(elem) = inner.next() {
        if elem.as_rule() == Rule::element {
            let symbol = elem.as_str().to_string();
            let count_tok = inner
                .next()
                .ok_or_else(|| ParseError::Grammar(format!("element '{}' has no count", symbol)))?;
            let count = parse_float(&count_tok)?;
            result.push(ElementCount { symbol, count });
        }
    }
    Ok(result)
}

fn collect_floats(pair: &Pairs<'_>) -> Result<Vec<f64>, ParseError> {
    pair.clone()
        .into_inner()
        .filter(|p| p.as_rule() == Rule::float_val)
        .map(|p| parse_float(&p))
        .collect()
}

fn parse_crystal(pair: Pairs<'_>) -> Result<CrystalConfig, ParseError> {
    let mut c = CrystalConfig::default();

    for line in pair.into_inner() {
        match line.as_rule() {
            Rule::crystal_type => {
                let val = line.into_inner().next().unwrap();
                c.crystal_type = Some(match val.as_str().to_lowercase().as_str() {
                    "cuboid" => CrystalType::Cuboid,
                    "polyhedron" | "obj" => CrystalType::Polyhedron,
                    "cylinder" => CrystalType::Cylinder,
                    "sphericalnew" => CrystalType::SphericalNew,
                    "spherical" => CrystalType::Spherical,
                    other => {
                        return Err(ParseError::Grammar(format!(
                            "Unknown crystal type: '{other}'"
                        )))
                    }
                });
            }
            Rule::crystal_ddm => {
                let kw = line.into_inner().next().unwrap();
                c.ddm = Some(match kw.as_str().to_lowercase().as_str() {
                    "simple" => DdmType::Simple,
                    "linear" => DdmType::Linear,
                    "leal" => DdmType::Leal,
                    "bfactor" => DdmType::Bfactor,
                    _ => DdmType::Simple,
                });
            }
            Rule::crystal_coefcalc => {
                let kw = line.into_inner().next().unwrap();
                c.coefcalc = Some(match kw.as_str().to_lowercase().as_str() {
                    "dummy" | "average" => CoefCalcType::Average,
                    "default" | "rdjava" | "rd3d" | "rd3" => CoefCalcType::Default,
                    "rdv2" | "rdv3" => CoefCalcType::RdFortran,
                    "exp" => CoefCalcType::Pdb,
                    "saxs" => CoefCalcType::Saxs,
                    "sequence" => CoefCalcType::Sequence,
                    "saxsseq" => CoefCalcType::SaxsSeq,
                    "smallmole" => CoefCalcType::SmallMole,
                    "expsm" => CoefCalcType::Cif,
                    "microed" => CoefCalcType::MicroED,
                    _ => CoefCalcType::Default,
                });
            }
            Rule::crystal_dim => {
                let floats = collect_floats(&line)?;
                if let Some(&x) = floats.first() {
                    c.dim_x = Some(x);
                }
                if let Some(&y) = floats.get(1) {
                    c.dim_y = Some(y);
                }
                if let Some(&z) = floats.get(2) {
                    c.dim_z = Some(z);
                }
            }
            Rule::crystal_ppm => {
                let f = line.into_inner().next().unwrap();
                c.pixels_per_micron = Some(parse_float(&f)?);
            }
            Rule::crystal_angle_p => {
                let f = line.into_inner().next().unwrap();
                c.angle_p = Some(parse_float(&f)?);
            }
            Rule::crystal_angle_l => {
                let f = line.into_inner().next().unwrap();
                c.angle_l = Some(parse_float(&f)?);
            }
            Rule::crystal_decay_param => {
                let floats = collect_floats(&line)?;
                c.gamma_param = floats.first().copied();
                c.b0_param = floats.get(1).copied();
                c.beta_param = floats.get(2).copied();
            }
            Rule::unitcell => {
                let floats = collect_floats(&line)?;
                if let Some(&a) = floats.first() {
                    c.cell_a = Some(a);
                }
                if let Some(&b) = floats.get(1) {
                    c.cell_b = Some(b);
                }
                if let Some(&cc) = floats.get(2) {
                    c.cell_c = Some(cc);
                }
                if let Some(&al) = floats.get(3) {
                    c.cell_alpha = Some(al);
                }
                if let Some(&be) = floats.get(4) {
                    c.cell_beta = Some(be);
                }
                if let Some(&ga) = floats.get(5) {
                    c.cell_gamma = Some(ga);
                }
            }
            Rule::num_monomers => {
                let f = line.into_inner().next().unwrap();
                c.num_monomers = Some(parse_int(&f)?);
            }
            Rule::num_residues => {
                let f = line.into_inner().next().unwrap();
                c.num_residues = Some(parse_int(&f)?);
            }
            Rule::num_rna => {
                let f = line.into_inner().next().unwrap();
                c.num_rna = Some(parse_int(&f)?);
            }
            Rule::num_dna => {
                let f = line.into_inner().next().unwrap();
                c.num_dna = Some(parse_int(&f)?);
            }
            Rule::num_carb => {
                let f = line.into_inner().next().unwrap();
                c.num_carb = Some(parse_int(&f)?);
            }
            Rule::heavy_protein_atoms => {
                c.heavy_protein_atoms = parse_element_counts(line)?;
            }
            Rule::small_mole_atoms => {
                c.small_mole_atoms = parse_element_counts(line)?;
            }
            Rule::heavy_solution_conc => {
                c.heavy_solution_conc = parse_element_counts(line)?;
            }
            Rule::solvent_fraction => {
                let f = line.into_inner().next().unwrap();
                c.solvent_fraction = Some(parse_float(&f)?);
            }
            Rule::pdb => {
                let val = line.into_inner().next().unwrap();
                c.pdb = Some(val.as_str().to_string());
            }
            Rule::cif => {
                let val = line.into_inner().next().unwrap();
                c.cif = Some(val.as_str().to_string());
            }
            Rule::wireframe_type => {
                let val = line.into_inner().next().unwrap();
                c.wireframe_type = Some(val.as_str().to_string());
            }
            Rule::model_file => {
                let val = line.into_inner().next().unwrap();
                c.model_file = Some(val.as_str().to_string());
            }
            Rule::calculate_pe_escape => {
                let val = line.into_inner().next().unwrap();
                c.calculate_pe_escape = Some(val.as_str().to_string());
            }
            Rule::calculate_fl_escape => {
                let val = line.into_inner().next().unwrap();
                c.calculate_fl_escape = Some(val.as_str().to_string());
            }
            Rule::fl_resolution => {
                let f = line.into_inner().next().unwrap();
                c.fl_resolution = Some(parse_int(&f)?);
            }
            Rule::pe_resolution => {
                let f = line.into_inner().next().unwrap();
                c.pe_resolution = Some(parse_int(&f)?);
            }
            Rule::protein_concentration => {
                let f = line.into_inner().next().unwrap();
                c.protein_conc = Some(parse_float(&f)?);
            }
            Rule::container_material_elements => {
                c.container_elements = parse_element_counts(line)?;
            }
            Rule::container_thickness => {
                let f = line.into_inner().next().unwrap();
                c.container_thickness = Some(parse_float(&f)?);
            }
            Rule::container_density => {
                let f = line.into_inner().next().unwrap();
                c.container_density = Some(parse_float(&f)?);
            }
            Rule::container_material_mixture => {
                let val = line.into_inner().next().unwrap();
                c.container_mixture = Some(val.as_str().to_string());
            }
            Rule::crystal_container_material => {
                let kw = line.into_inner().next().unwrap();
                c.container_material = Some(match kw.as_str().to_lowercase().as_str() {
                    "none" => ContainerMaterialType::None,
                    "mixture" => ContainerMaterialType::Mixture,
                    "elemental" => ContainerMaterialType::Elemental,
                    _ => ContainerMaterialType::None,
                });
            }
            Rule::sequence_file => {
                let val = line.into_inner().next().unwrap();
                c.seq_file = Some(val.as_str().to_string());
            }
            Rule::surrounding_heavy_conc => {
                c.surrounding_heavy_conc = parse_element_counts(line)?;
            }
            Rule::surrounding_thickness => {
                let floats = collect_floats(&line)?;
                if floats.len() >= 3 {
                    c.surrounding_thickness = Some([floats[0], floats[1], floats[2]]);
                }
            }
            Rule::oil_based => {
                let val = line.into_inner().next().unwrap();
                c.oil_based = Some(val.as_str().to_string());
            }
            Rule::calc_surrounding => {
                let val = line.into_inner().next().unwrap();
                c.calc_surrounding = Some(val.as_str().to_string());
            }
            Rule::oil_elements => {
                c.oil_elements = parse_element_counts(line)?;
            }
            Rule::oil_density => {
                let f = line.into_inner().next().unwrap();
                c.oil_density = Some(parse_float(&f)?);
            }
            Rule::goniometer_axis => {
                let f = line.into_inner().next().unwrap();
                c.goniometer_axis = Some(parse_float(&f)?);
            }
            Rule::polarisation_direction => {
                let f = line.into_inner().next().unwrap();
                c.polarisation_direction = Some(parse_float(&f)?);
            }
            Rule::program => {
                let val = line.into_inner().next().unwrap();
                c.program = Some(val.as_str().to_string());
            }
            Rule::sim_electrons => {
                let f = line.into_inner().next().unwrap();
                c.sim_electrons = Some(parse_long(&f)?);
            }
            Rule::runs => {
                let f = line.into_inner().next().unwrap();
                c.runs = Some(parse_int(&f)?);
            }
            _ => {}
        }
    }

    // Default coefcalc to Default if not specified
    if c.coefcalc.is_none() {
        c.coefcalc = Some(CoefCalcType::Default);
    }

    Ok(c)
}

fn parse_beam(pair: Pairs<'_>) -> Result<BeamConfig, ParseError> {
    let mut b = BeamConfig::default();

    for line in pair.into_inner() {
        match line.as_rule() {
            Rule::beam_type => {
                let val = line.into_inner().next().unwrap();
                b.beam_type = Some(match val.as_str().to_lowercase().as_str() {
                    "gaussian" => BeamType::Gaussian,
                    "tophat" => BeamType::Tophat,
                    "experimental" => BeamType::Experimental,
                    "experimentalpgm" => BeamType::ExperimentalPgm,
                    other => {
                        return Err(ParseError::Grammar(format!("Unknown beam type: '{other}'")))
                    }
                });
            }
            Rule::beam_flux => {
                let f = line.into_inner().next().unwrap();
                b.flux = Some(parse_float(&f)?);
            }
            Rule::beam_fwhm => {
                let floats = collect_floats(&line)?;
                b.fwhm_x = floats.first().copied();
                b.fwhm_y = floats.get(1).copied();
            }
            Rule::beam_energy => {
                let f = line.into_inner().next().unwrap();
                b.energy = Some(parse_float(&f)?);
            }
            Rule::beam_exposure => {
                let f = line.into_inner().next().unwrap();
                b.exposure = Some(parse_float(&f)?);
            }
            Rule::beam_file => {
                let val = line.into_inner().next().unwrap();
                b.file = Some(val.as_str().to_string());
            }
            Rule::beam_pixel_size => {
                let floats = collect_floats(&line)?;
                b.pixel_size_x = floats.first().copied();
                b.pixel_size_y = floats.get(1).copied();
            }
            Rule::beam_semi_angle => {
                let f = line.into_inner().next().unwrap();
                b.semi_angle = Some(parse_float(&f)?);
            }
            Rule::beam_aperture_radius => {
                let f = line.into_inner().next().unwrap();
                b.aperture_radius = Some(parse_float(&f)?);
            }
            Rule::image_dimensions => {
                let floats = collect_floats(&line)?;
                b.image_x = floats.first().copied();
                b.image_y = floats.get(1).copied();
            }
            Rule::pulse_energy => {
                let f = line.into_inner().next().unwrap();
                b.pulse_energy = Some(parse_float(&f)?);
            }
            Rule::energy_fwhm => {
                let f = line.into_inner().next().unwrap();
                b.energy_fwhm = Some(parse_float(&f)?);
            }
            Rule::beam_collimation => {
                let text = line.as_str().to_lowercase();
                let floats = collect_floats(&line)?;
                if text.starts_with("rectangular") {
                    b.collimation = Some(Collimation::Rectangular {
                        h: floats[0],
                        v: floats[1],
                    });
                } else if text.starts_with("circular") {
                    b.collimation = Some(Collimation::Circular {
                        h: floats[0],
                        v: floats[1],
                    });
                } else if text.starts_with("horizontal") {
                    b.collimation = Some(Collimation::Horizontal { h: floats[0] });
                } else if text.starts_with("vertical") {
                    b.collimation = Some(Collimation::Vertical { v: floats[0] });
                } else {
                    b.collimation = Some(Collimation::Unspecified);
                }
            }
            _ => {}
        }
    }

    Ok(b)
}

fn parse_wedge(pair: Pairs<'_>) -> Result<WedgeConfig, ParseError> {
    let mut w = WedgeConfig::default();
    let mut inner = pair.into_inner();

    // First two floats are start and end angles
    if let Some(start) = inner.next() {
        w.start_ang = parse_float(&start)?;
    }
    if let Some(end) = inner.next() {
        w.end_ang = parse_float(&end)?;
    }

    for line in inner {
        match line.as_rule() {
            Rule::wedge_exposure => {
                let f = line.into_inner().next().unwrap();
                w.exposure_time = Some(parse_float(&f)?);
            }
            Rule::wedge_ang_res => {
                let f = line.into_inner().next().unwrap();
                w.angular_resolution = Some(parse_float(&f)?);
            }
            Rule::wedge_start_offset => {
                let floats = collect_floats(&line)?;
                w.start_offset_x = floats.first().copied();
                w.start_offset_y = floats.get(1).copied();
                w.start_offset_z = floats.get(2).copied();
            }
            Rule::wedge_translate => {
                let floats = collect_floats(&line)?;
                w.translate_x = floats.first().copied();
                w.translate_y = floats.get(1).copied();
                w.translate_z = floats.get(2).copied();
            }
            Rule::wedge_rot_ax_beam_offset => {
                let f = line.into_inner().next().unwrap();
                w.rot_ax_beam_offset = Some(parse_float(&f)?);
            }
            Rule::wedge_max_res => {
                let f = line.into_inner().next().unwrap();
                w.max_resolution = Some(parse_float(&f)?);
            }
            _ => {}
        }
    }

    Ok(w)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_input() {
        let input = r#"
Crystal
Type Cuboid
Dimensions 0.2 0.2 0.2
PixelsPerMicron 100
AbsCoefCalc SMALLMOLE
UnitCell 10.012 18.808 6.3340
SmallMoleAtoms Mg 1 O 3
NumMonomers 8

Beam
Type Gaussian
Flux 3.8e12
FWHM 10 10
Energy 12.4
Collimation Circular 30 30

Wedge 0 360
ExposureTime 60
"#;
        let config = parse(input).unwrap();
        assert_eq!(config.crystals.len(), 1);
        assert_eq!(config.beams.len(), 1);
        assert_eq!(config.wedges.len(), 1);

        let c = &config.crystals[0];
        assert_eq!(c.crystal_type, Some(CrystalType::Cuboid));
        assert_eq!(c.coefcalc, Some(CoefCalcType::SmallMole));
        assert_eq!(c.dim_x, Some(0.2));
        assert_eq!(c.dim_y, Some(0.2));
        assert_eq!(c.dim_z, Some(0.2));
        assert_eq!(c.pixels_per_micron, Some(100.0));
        assert_eq!(c.cell_a, Some(10.012));
        assert_eq!(c.num_monomers, Some(8));
        assert_eq!(c.small_mole_atoms.len(), 2);
        assert_eq!(c.small_mole_atoms[0].symbol, "Mg");
        assert_eq!(c.small_mole_atoms[1].symbol, "O");
        assert_eq!(c.small_mole_atoms[1].count, 3.0);

        let b = &config.beams[0];
        assert_eq!(b.beam_type, Some(BeamType::Gaussian));
        assert_eq!(b.flux, Some(3.8e12));
        assert_eq!(b.fwhm_x, Some(10.0));
        assert_eq!(b.energy, Some(12.4));
        assert!(matches!(
            b.collimation,
            Some(Collimation::Circular { h: 30.0, v: 30.0 })
        ));

        let w = &config.wedges[0];
        assert_eq!(w.start_ang, 0.0);
        assert_eq!(w.end_ang, 360.0);
        assert_eq!(w.exposure_time, Some(60.0));
    }

    #[test]
    fn test_parse_with_comments() {
        let input = r#"
# This is a comment
Crystal
Type Cuboid
Dimensions 100 100 100  # crystal size in microns
// Another comment
! Yet another
PixelsPerMicron 0.5

Beam
Type TopHat
Energy 12.40
Flux 6.00e12

Wedge 0 90
ExposureTime 100
"#;
        let config = parse(input).unwrap();
        assert_eq!(config.crystals.len(), 1);
        assert_eq!(config.crystals[0].dim_x, Some(100.0));
    }

    #[test]
    fn test_parse_pdb_input() {
        let input = r#"
Crystal
Type Cuboid
Dimensions 100 100 100
AbsCoefCalc EXP
Pdb 1dwa
SolventHeavyConc S 2700
PixelsPerMicron 0.2

Beam
Type TopHat
Energy 12.40
Collimation Rectangular 160 160
Flux 6.00e12

Wedge 0 90
exposureTime 100
"#;
        let config = parse(input).unwrap();
        let c = &config.crystals[0];
        assert_eq!(c.coefcalc, Some(CoefCalcType::Pdb));
        assert_eq!(c.pdb.as_deref(), Some("1dwa"));
        assert_eq!(c.heavy_solution_conc.len(), 1);
        assert_eq!(c.heavy_solution_conc[0].symbol, "S");
        assert_eq!(c.heavy_solution_conc[0].count, 2700.0);
    }

    #[test]
    fn test_parse_cif_input() {
        let input = r#"
Crystal
Type Cuboid
Dimensions 20 20 20
PixelsPerMicron 0.5
AbsCoefCalc EXPSM
CIF Fe3O4

Beam
Type Tophat
Flux 2e12
Energy 30
Collimation Rectangular 5 5

Wedge 0 0
ExposureTime 2
"#;
        let config = parse(input).unwrap();
        let c = &config.crystals[0];
        assert_eq!(c.coefcalc, Some(CoefCalcType::Cif));
        assert_eq!(c.cif.as_deref(), Some("Fe3O4"));
    }

    #[test]
    fn test_abscoefcalc_cif_keyword_rejected() {
        // Java's grammar token named CIF matches the text "EXPSM", not "CIF".
        // "AbsCoefCalc CIF" is therefore a Java parse error and must be rejected here too.
        let input = r#"
Crystal
Type Cuboid
Dimensions 20 20 20
AbsCoefCalc CIF
CIF Fe3O4

Beam
Type Tophat
Flux 2e12
Energy 12.4
Collimation Rectangular 5 5

Wedge 0 0
ExposureTime 1
"#;
        assert!(
            parse(input).is_err(),
            "AbsCoefCalc CIF must be rejected; Java keyword is EXPSM"
        );
    }

    #[test]
    fn test_parse_saxs_cylinder() {
        let input = r#"
Crystal
Type Cylinder
Dimensions 1700 1000
PixelsPerMicron 0.01
CoefCalc SAXSseq
SeqFile 4OR0.fasta
ProteinConc 2
ContainerMaterialType elemental
MaterialElements Si 1 O 2
ContainerThickness 50
ContainerDensity 2.648

Beam
Type Gaussian
Flux 2e12
FWHM 700 700
Energy 12.1
Collimation Rectangular 1000 1000

Wedge 0 0
ExposureTime 50
"#;
        let config = parse(input).unwrap();
        let c = &config.crystals[0];
        assert_eq!(c.crystal_type, Some(CrystalType::Cylinder));
        assert_eq!(c.dim_x, Some(1700.0));
        assert_eq!(c.dim_y, Some(1000.0));
        assert_eq!(c.dim_z, None);
        assert_eq!(c.coefcalc, Some(CoefCalcType::SaxsSeq));
        assert_eq!(c.seq_file.as_deref(), Some("4OR0.fasta"));
        assert_eq!(c.protein_conc, Some(2.0));
        assert_eq!(c.container_material, Some(ContainerMaterialType::Elemental));
        assert_eq!(c.container_elements.len(), 2);
        assert_eq!(c.container_elements[0].symbol, "Si");
        assert_eq!(c.container_thickness, Some(50.0));
        assert_eq!(c.container_density, Some(2.648));
    }

    #[test]
    fn test_parse_unitcell_with_angles() {
        let input = r#"
Crystal
Type Cuboid
Dimensions 100 100 100
UnitCell 78.0 78.0 78.0 90 90 90
NumResidues 200
NumMonomers 4
PixelsPerMicron 0.5

Beam
Type Gaussian
Flux 2e12
FWHM 100 100
Energy 12.4

Wedge 0 90
ExposureTime 100
"#;
        let config = parse(input).unwrap();
        let c = &config.crystals[0];
        assert_eq!(c.cell_a, Some(78.0));
        assert_eq!(c.cell_alpha, Some(90.0));
        assert_eq!(c.cell_beta, Some(90.0));
        assert_eq!(c.cell_gamma, Some(90.0));
    }
}
