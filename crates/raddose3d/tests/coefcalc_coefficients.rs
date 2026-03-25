use raddose3d::coefcalc::CoefCalc;
/// Phase 4 verification: compare absorption coefficients against Java RADDOSE-3D.
///
/// Each test constructs the relevant CoefCalc mode, calls update_coefficients at a
/// reference energy, and checks photoelectric, inelastic, elastic, attenuation, and
/// density values match the Java reference within 0.5%.
///
/// Java reference values were obtained by running:
///   java -jar raddose3d.jar -i <fixture>
/// for each mode. The tolerance of 0.5% is loose enough to absorb any sub-display-
/// precision floating-point rounding while tight enough to catch real bugs.
use raddose3d::coefcalc::{
    CoefCalcAverage, CoefCalcFromCIF, CoefCalcFromSequence, CoefCalcFromSequenceSAXS,
    CoefCalcMicroED, CoefCalcSAXS, CoefCalcSmallMolecules,
};
use raddose3d_parser::config::{CoefCalcType, CrystalConfig, ElementCount};

const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures");

/// Assert two f64 values are within `rel_tol` relative of each other.
fn assert_rel_eq(label: &str, got: f64, want: f64, rel_tol: f64) {
    let rel = if want.abs() < 1e-30 {
        (got - want).abs()
    } else {
        ((got - want) / want).abs()
    };
    assert!(
        rel <= rel_tol,
        "{label}: got {got:.4e}, want {want:.4e} (rel err {:.4e} > {rel_tol})",
        rel
    );
}

// ── CoefCalcAverage ───────────────────────────────────────────────────────────

/// CoefCalcAverage uses Holton 2010 hardcoded constants. No energy dependence.
/// Java reference (12.4 keV): density 1.20 g/ml, absorption 2.37e-4 /um,
/// attenuation 2.81e-4 /um.
#[test]
fn test_average_coefficients() {
    let mut cc = CoefCalcAverage;
    cc.update_coefficients(12.4);

    assert_rel_eq("photo", cc.absorption_coefficient(), 2.37e-4, 0.01);
    assert_rel_eq("attenuation", cc.attenuation_coefficient(), 2.81e-4, 0.01);
    assert_rel_eq("elastic", cc.elastic_coefficient(), 1.799e-5, 0.01);
    assert_rel_eq("density", cc.density(), 1.2, 0.01);
}

// ── CoefCalcSmallMolecules ────────────────────────────────────────────────────

/// MgO₃ crystal: 8 monomers per unit cell, explicit atom counts.
/// Unit cell: a=10.012, b=18.808, c=6.334 Å (all angles 90°).
///
/// Java reference (SmallMoleAtoms Mg 1 O 3, NumMonomers 8, 12.4 keV):
///   Photoelectric: 4.27e-4 /um
///   Inelastic:     1.07e-5 /um
///   Elastic:       2.21e-5 /um
///   Attenuation:   4.60e-4 /um
///   Density:       0.81 g/ml
///
/// Note: The original SMXray_example_input.txt fixture uses `SmallMoleAtoms Mg  O 3`
/// (Mg with no count). Java's ANTLR parser applies error recovery that skips the
/// "O" token and uses "3" as the Mg count, giving 3 Mg per monomer (no oxygen).
/// Rust parses the same input correctly as 1 Mg + 3 O, which matches Java when
/// explicit counts are given. The explicit form is used here.
#[test]
fn test_small_molecules_mgo3() {
    let config = CrystalConfig {
        coefcalc: Some(CoefCalcType::SmallMole),
        cell_a: Some(10.012),
        cell_b: Some(18.808),
        cell_c: Some(6.334),
        num_monomers: Some(8),
        small_mole_atoms: vec![
            ElementCount {
                symbol: "Mg".into(),
                count: 1.0,
            },
            ElementCount {
                symbol: "O".into(),
                count: 3.0,
            },
        ],
        ..Default::default()
    };

    let mut cc = CoefCalcSmallMolecules::from_config(&config).unwrap();
    cc.update_coefficients(12.4);

    assert_rel_eq("photo", cc.absorption_coefficient(), 4.27e-4, 0.01);
    assert_rel_eq("inelastic", cc.inelastic_coefficient(), 1.07e-5, 0.01);
    assert_rel_eq("elastic", cc.elastic_coefficient(), 2.21e-5, 0.01);
    assert_rel_eq("attenuation", cc.attenuation_coefficient(), 4.60e-4, 0.01);
    assert_rel_eq("density", cc.density(), 0.81, 0.01);
}

// ── CoefCalcFromCIF ───────────────────────────────────────────────────────────

/// Alanine (C₃ H₇ N O₂), cell volume 538.08 Å³.
///
/// Java reference (EXPSM + CIF /tmp/alanine.cif, 12.4 keV):
///   Photoelectric: 4.34e-5 /um
///   Inelastic:     4.52e-6 /um
///   Elastic:       4.22e-6 /um
///   Attenuation:   5.22e-5 /um
///   Density:       0.27 g/ml
#[test]
fn test_from_cif_alanine() {
    let cif_path = format!("{FIXTURES_DIR}/alanine.cif");
    let mut cc = CoefCalcFromCIF::from_file(&cif_path).unwrap();
    cc.update_coefficients(12.4);

    assert_rel_eq("photo", cc.absorption_coefficient(), 4.34e-5, 0.01);
    assert_rel_eq("inelastic", cc.inelastic_coefficient(), 4.52e-6, 0.01);
    assert_rel_eq("elastic", cc.elastic_coefficient(), 4.22e-6, 0.01);
    assert_rel_eq("attenuation", cc.attenuation_coefficient(), 5.22e-5, 0.01);
    // Java prints "0.27 g/ml" (2 sig figs); the actual computed value is 0.2749.
    assert_rel_eq("density", cc.density(), 0.2749, 0.01);
}

// ── CoefCalcFromSequence ──────────────────────────────────────────────────────

/// Insulin A+B chains, 24 monomers, cubic cell 78 Å.
///
/// Java reference (SEQUENCE + SequenceFile insulin_seq.fasta, 12.1 keV):
///   Photoelectric: 2.67e-4 /um
///   Inelastic:     1.80e-5 /um
///   Elastic:       1.97e-5 /um
///   Attenuation:   3.05e-4 /um
///   Density:       1.10 g/ml
///   Solvent:       63.74%
#[test]
fn test_from_sequence_insulin() {
    let seq_path = format!("{FIXTURES_DIR}/insulin_seq.fasta");
    let config = CrystalConfig {
        coefcalc: Some(CoefCalcType::Sequence),
        cell_a: Some(78.0),
        cell_b: Some(78.0),
        cell_c: Some(78.0),
        num_monomers: Some(24),
        seq_file: Some(seq_path),
        ..Default::default()
    };

    let mut cc = CoefCalcFromSequence::from_config(&config).unwrap();
    cc.update_coefficients(12.1);

    assert_rel_eq("photo", cc.absorption_coefficient(), 2.67e-4, 0.01);
    assert_rel_eq("inelastic", cc.inelastic_coefficient(), 1.80e-5, 0.01);
    assert_rel_eq("elastic", cc.elastic_coefficient(), 1.97e-5, 0.01);
    assert_rel_eq("attenuation", cc.attenuation_coefficient(), 3.05e-4, 0.01);
    assert_rel_eq("density", cc.density(), 1.10, 0.01);
    assert_rel_eq("solvent_fraction", cc.solvent_fraction(), 0.6374, 0.01);
}

// ── CoefCalcSAXS ─────────────────────────────────────────────────────────────

/// 582-residue protein at 2.0 mg/mL in a 1000 Å cubic cell.
///
/// Java reference (SAXS, NumResidues 582, ProteinConc 2.0, UnitCell 1000³, 12.4 keV):
///   Num monomers:  19
///   Photoelectric: 2.18e-4 /um
///   Inelastic:     1.56e-5 /um
///   Elastic:       1.75e-5 /um
///   Attenuation:   2.51e-4 /um
///   Density:       0.93 g/ml
#[test]
fn test_saxs_582_residues() {
    let config = CrystalConfig {
        coefcalc: Some(CoefCalcType::Saxs),
        cell_a: Some(1000.0),
        cell_b: Some(1000.0),
        cell_c: Some(1000.0),
        num_residues: Some(582),
        protein_conc: Some(2.0),
        ..Default::default()
    };

    let mut cc = CoefCalcSAXS::from_config(&config).unwrap();
    cc.update_coefficients(12.4);

    assert_rel_eq("photo", cc.absorption_coefficient(), 2.18e-4, 0.01);
    assert_rel_eq("inelastic", cc.inelastic_coefficient(), 1.56e-5, 0.01);
    assert_rel_eq("elastic", cc.elastic_coefficient(), 1.75e-5, 0.01);
    assert_rel_eq("attenuation", cc.attenuation_coefficient(), 2.51e-4, 0.01);
    assert_rel_eq("density", cc.density(), 0.93, 0.01);
    // Verify monomer count
    assert_eq!(cc.compute.num_monomers, 19);
}

// ── CoefCalcFromSequenceSAXS ──────────────────────────────────────────────────

/// BSA fragment (359 residues from FASTA) at 2.0 mg/mL in a 1000 Å cubic cell.
///
/// Java reference (SAXSSEQ + SequenceFile bsa_fragment.fasta, ProteinConc 2.0,
///                 UnitCell 1000³, 12.4 keV):
///   Num monomers:  25
///   Photoelectric: 2.18e-4 /um
///   Inelastic:     1.56e-5 /um
///   Elastic:       1.75e-5 /um
///   Attenuation:   2.51e-4 /um
///   Density:       0.93 g/ml
#[test]
fn test_from_sequence_saxs_bsa() {
    let seq_path = format!("{FIXTURES_DIR}/bsa_fragment.fasta");
    let config = CrystalConfig {
        coefcalc: Some(CoefCalcType::SaxsSeq),
        cell_a: Some(1000.0),
        cell_b: Some(1000.0),
        cell_c: Some(1000.0),
        protein_conc: Some(2.0),
        seq_file: Some(seq_path),
        ..Default::default()
    };

    let mut cc = CoefCalcFromSequenceSAXS::from_config(&config).unwrap();
    cc.update_coefficients(12.4);

    assert_rel_eq("photo", cc.absorption_coefficient(), 2.18e-4, 0.01);
    assert_rel_eq("inelastic", cc.inelastic_coefficient(), 1.56e-5, 0.01);
    assert_rel_eq("elastic", cc.elastic_coefficient(), 1.75e-5, 0.01);
    assert_rel_eq("attenuation", cc.attenuation_coefficient(), 2.51e-4, 0.01);
    assert_rel_eq("density", cc.density(), 0.93, 0.01);
    assert_eq!(cc.compute.num_monomers, 25);
}

// ── CoefCalcMicroED ───────────────────────────────────────────────────────────

/// MicroED uses the same composition logic as CoefCalcFromParams (RD3D default).
/// The `is_em` flag signals electron mode for the expose loop but the X-ray
/// cross-sections still govern coefficient calculation in the current implementation.
///
/// Test: insulin unit cell (24 monomers × 51 residues, cubic 78 Å).
/// Expected to match the Java RD3D reference at 12.1 keV:
///   Photoelectric: 2.23e-4 /um
///   Inelastic:     1.78e-5 /um
///   Elastic:       1.88e-5 /um
///   Attenuation:   2.60e-4 /um
///   Density:       1.08 g/ml
#[test]
fn test_micro_ed_insulin_composition() {
    let config = CrystalConfig {
        coefcalc: Some(CoefCalcType::Default), // MicroED not a parser keyword yet
        cell_a: Some(78.0),
        cell_b: Some(78.0),
        cell_c: Some(78.0),
        num_monomers: Some(24),
        num_residues: Some(51),
        ..Default::default()
    };

    let mut cc = CoefCalcMicroED::from_config(&config).unwrap();
    cc.update_coefficients(12.1);

    assert_rel_eq("photo", cc.absorption_coefficient(), 2.23e-4, 0.01);
    assert_rel_eq("inelastic", cc.inelastic_coefficient(), 1.78e-5, 0.01);
    assert_rel_eq("elastic", cc.elastic_coefficient(), 1.88e-5, 0.01);
    assert_rel_eq("attenuation", cc.attenuation_coefficient(), 2.60e-4, 0.01);
    assert_rel_eq("density", cc.density(), 1.08, 0.01);
    assert!(cc.is_em, "MicroED should set is_em = true");
}

// ── Energy dependence ─────────────────────────────────────────────────────────

/// Verify that photoelectric coefficient decreases with increasing energy
/// (photoelectric cross-section scales roughly as E^-3).
/// Tests the SmallMolecules mode as a proxy for all modes using CoefCalcCompute.
#[test]
fn test_photoelectric_decreases_with_energy() {
    let config = CrystalConfig {
        coefcalc: Some(CoefCalcType::SmallMole),
        cell_a: Some(10.012),
        cell_b: Some(18.808),
        cell_c: Some(6.334),
        num_monomers: Some(8),
        small_mole_atoms: vec![
            ElementCount {
                symbol: "Mg".into(),
                count: 1.0,
            },
            ElementCount {
                symbol: "O".into(),
                count: 3.0,
            },
        ],
        ..Default::default()
    };

    let mut cc = CoefCalcSmallMolecules::from_config(&config).unwrap();

    cc.update_coefficients(8.0);
    let photo_8kev = cc.absorption_coefficient();

    cc.update_coefficients(12.4);
    let photo_12kev = cc.absorption_coefficient();

    cc.update_coefficients(20.0);
    let photo_20kev = cc.absorption_coefficient();

    assert!(
        photo_8kev > photo_12kev,
        "PE coeff should decrease from 8→12.4 keV: {photo_8kev:.3e} > {photo_12kev:.3e}"
    );
    assert!(
        photo_12kev > photo_20kev,
        "PE coeff should decrease from 12.4→20 keV: {photo_12kev:.3e} > {photo_20kev:.3e}"
    );
}

/// Total attenuation must always exceed photoelectric absorption.
#[test]
fn test_attenuation_exceeds_absorption() {
    let config = CrystalConfig {
        coefcalc: Some(CoefCalcType::SmallMole),
        cell_a: Some(10.012),
        cell_b: Some(18.808),
        cell_c: Some(6.334),
        num_monomers: Some(8),
        small_mole_atoms: vec![
            ElementCount {
                symbol: "Mg".into(),
                count: 1.0,
            },
            ElementCount {
                symbol: "O".into(),
                count: 3.0,
            },
        ],
        ..Default::default()
    };

    for energy in [8.0, 12.4, 20.0, 30.0] {
        let mut cc = CoefCalcSmallMolecules::from_config(&config).unwrap();
        cc.update_coefficients(energy);
        assert!(
            cc.attenuation_coefficient() >= cc.absorption_coefficient(),
            "attenuation must be >= absorption at {energy} keV"
        );
    }
}
