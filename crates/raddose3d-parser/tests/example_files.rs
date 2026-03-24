use raddose3d_parser::*;

// Java project examples are at /workspace/examples/
const EXAMPLES_DIR: &str = "/workspace/examples";

#[test]
fn test_smxray_example() {
    // Note: The original SMXray_example_input.txt has "SmallMoleAtoms Mg  O 3"
    // which is malformed (ANTLR grammar also requires element-count pairs).
    // We fix it to "SmallMoleAtoms Mg 1 O 3" for the test.
    let input = std::fs::read_to_string(
        format!("{}/SMXray_example_input.txt", EXAMPLES_DIR)
    ).unwrap();
    let input = input.replace("SmallMoleAtoms Mg  O 3", "SmallMoleAtoms Mg 1 O 3");
    let config = parse(&input).unwrap();
    assert_eq!(config.crystals.len(), 1);
    assert_eq!(config.beams.len(), 1);
    assert_eq!(config.wedges.len(), 1);

    let c = &config.crystals[0];
    assert_eq!(c.crystal_type.as_deref(), Some("Cuboid"));
    assert_eq!(c.coefcalc, Some(CoefCalcType::SmallMole));
    assert_eq!(c.dim_x, Some(0.2));
    assert_eq!(c.pixels_per_micron, Some(100.0));
}

#[test]
fn test_smxray2_example() {
    let input = std::fs::read_to_string(
        format!("{}/SMXray2_example_input.txt", EXAMPLES_DIR)
    ).unwrap();
    let config = parse(&input).unwrap();
    assert_eq!(config.crystals.len(), 1);
    assert_eq!(config.beams.len(), 1);
    assert_eq!(config.wedges.len(), 1);

    let c = &config.crystals[0];
    assert_eq!(c.coefcalc, Some(CoefCalcType::Cif));
    assert_eq!(c.cif.as_deref(), Some("Fe3O4"));
}
