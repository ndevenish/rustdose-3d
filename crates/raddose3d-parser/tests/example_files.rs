use raddose3d_parser::*;

// Java project examples are at /workspace/examples/
const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures");

#[test]
fn test_smxray_example() {
    let input =
        std::fs::read_to_string(format!("{}/SMXray_example_input.txt", FIXTURES_DIR)).unwrap();
    let config = parse(&input).unwrap();
    assert_eq!(config.crystals.len(), 1);
    assert_eq!(config.beams.len(), 1);
    assert_eq!(config.wedges.len(), 1);

    let c = &config.crystals[0];
    assert_eq!(c.crystal_type, Some(CrystalType::Cuboid));
    assert_eq!(c.coefcalc, Some(CoefCalcType::SmallMole));
    assert_eq!(c.dim_x, Some(0.2));
    assert_eq!(c.pixels_per_micron, Some(100.0));
    assert_eq!(c.small_mole_atoms.len(), 2);
    assert_eq!(c.small_mole_atoms[0].symbol, "Mg");
    assert_eq!(c.small_mole_atoms[0].count, 1.0);
    assert_eq!(c.small_mole_atoms[1].symbol, "O");
    assert_eq!(c.small_mole_atoms[1].count, 3.0);
}

#[test]
fn test_small_mole_atoms_missing_count_is_rejected() {
    // Java's ANTLR parser silently deletes the second element symbol and assigns
    // its count to the first element (e.g. "Fe O 4" → Fe×4, O dropped entirely).
    // Rust must reject this outright — element counts are mandatory.
    let input = "
Crystal
Type Cuboid
Dimensions 10.0
AbsCoefCalc SmallMole
UnitCell 25.354 40.668 20.605
SmallMoleAtoms Fe O 4
NumMonomers 5

Beam
Type Tophat
Energy 12.4
Flux 1e12
FWHM 10 10

Wedge 0 90
ExposureTime 1
";
    assert!(
        parse(input).is_err(),
        "SmallMoleAtoms with missing count must be a parse error"
    );
}

#[test]
fn test_smxray2_example() {
    let input =
        std::fs::read_to_string(format!("{}/SMXray2_example_input.txt", FIXTURES_DIR)).unwrap();
    let config = parse(&input).unwrap();
    assert_eq!(config.crystals.len(), 1);
    assert_eq!(config.beams.len(), 1);
    assert_eq!(config.wedges.len(), 1);

    let c = &config.crystals[0];
    assert_eq!(c.coefcalc, Some(CoefCalcType::Cif));
    assert_eq!(c.cif.as_deref(), Some("Fe3O4"));
}
