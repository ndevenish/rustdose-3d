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
    // "SmallMoleAtoms Mg  O 3" — Mg has no explicit count, defaults to 1
    assert_eq!(c.small_mole_atoms.len(), 2);
    assert_eq!(c.small_mole_atoms[0].symbol, "Mg");
    assert_eq!(c.small_mole_atoms[0].count, 1.0);
    assert_eq!(c.small_mole_atoms[1].symbol, "O");
    assert_eq!(c.small_mole_atoms[1].count, 3.0);
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
