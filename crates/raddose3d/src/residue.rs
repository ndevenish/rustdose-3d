/// Residue types.
pub const TYPE_PROTEIN: u8 = 1;
pub const TYPE_RNA: u8 = 2;
pub const TYPE_DNA: u8 = 3;

/// Atomic composition of a residue.
#[derive(Debug, Clone, Copy)]
pub struct Residue {
    pub residue_type: u8,
    pub carbons: f64,
    pub hydrogens: f64,
    pub oxygens: f64,
    pub nitrogens: f64,
    pub sulphurs: f64,
    pub phosphoruses: f64,
    pub seleniums: f64,
    pub molecular_weight: f64,
}

impl Residue {
    const fn new(
        residue_type: u8,
        c: i32, h: i32, o: i32, n: i32, s: i32, p: i32, se: i32,
        mw: f64,
    ) -> Self {
        Residue {
            residue_type,
            carbons: c as f64,
            hydrogens: h as f64,
            oxygens: o as f64,
            nitrogens: n as f64,
            sulphurs: s as f64,
            phosphoruses: p as f64,
            seleniums: se as f64,
            molecular_weight: mw,
        }
    }
}

// ─── Protein residues (3-letter codes) ────────────────────────────────────────

pub fn get_residue_by_three_letter(code: &str) -> Option<Residue> {
    match code.trim() {
        "ALA" => Some(Residue::new(TYPE_PROTEIN, 3, 7, 2, 1, 0, 0, 0, 89.09)),
        "ARG" => Some(Residue::new(TYPE_PROTEIN, 6, 15, 1, 4, 0, 0, 0, 174.20)),
        "ASN" => Some(Residue::new(TYPE_PROTEIN, 4, 8, 2, 2, 0, 0, 0, 132.12)),
        "ASP" => Some(Residue::new(TYPE_PROTEIN, 4, 5, 3, 1, 0, 0, 0, 133.10)),
        "CYS" => Some(Residue::new(TYPE_PROTEIN, 3, 5, 1, 1, 1, 0, 0, 121.16)),
        "GLN" => Some(Residue::new(TYPE_PROTEIN, 5, 8, 2, 2, 0, 0, 0, 146.15)),
        "GLU" => Some(Residue::new(TYPE_PROTEIN, 5, 7, 3, 1, 0, 0, 0, 147.13)),
        "GLY" => Some(Residue::new(TYPE_PROTEIN, 2, 3, 1, 1, 0, 0, 0, 75.07)),
        "HIS" => Some(Residue::new(TYPE_PROTEIN, 6, 7, 1, 3, 0, 0, 0, 155.16)),
        "ILE" => Some(Residue::new(TYPE_PROTEIN, 6, 11, 1, 1, 0, 0, 0, 131.17)),
        "LEU" => Some(Residue::new(TYPE_PROTEIN, 6, 11, 1, 1, 0, 0, 0, 131.17)),
        "LYS" => Some(Residue::new(TYPE_PROTEIN, 6, 12, 1, 2, 0, 0, 0, 146.19)),
        "MET" => Some(Residue::new(TYPE_PROTEIN, 5, 9, 1, 1, 1, 0, 0, 149.21)),
        "MSE" => Some(Residue::new(TYPE_PROTEIN, 5, 9, 1, 1, 0, 0, 1, 196.11)),
        "PHE" => Some(Residue::new(TYPE_PROTEIN, 9, 9, 1, 1, 0, 0, 0, 165.19)),
        "PRO" => Some(Residue::new(TYPE_PROTEIN, 5, 7, 1, 1, 0, 0, 0, 115.13)),
        "SER" => Some(Residue::new(TYPE_PROTEIN, 3, 5, 2, 1, 0, 0, 0, 105.09)),
        "THR" => Some(Residue::new(TYPE_PROTEIN, 4, 7, 2, 1, 0, 0, 0, 119.12)),
        "TRP" => Some(Residue::new(TYPE_PROTEIN, 11, 10, 1, 2, 0, 0, 0, 204.23)),
        "TYR" => Some(Residue::new(TYPE_PROTEIN, 9, 9, 2, 1, 0, 0, 0, 181.19)),
        "VAL" => Some(Residue::new(TYPE_PROTEIN, 5, 9, 1, 1, 0, 0, 0, 117.15)),
        // RNA (3-letter codes from PDB, padded with leading spaces in Java)
        "A"   => Some(Residue::new(TYPE_RNA,     10, 12, 6, 5, 0, 1, 0, 347.2)),
        "U"   => Some(Residue::new(TYPE_RNA,     9,  11, 8, 2, 0, 1, 0, 324.2)),
        "G"   => Some(Residue::new(TYPE_RNA,     10, 12, 7, 5, 0, 1, 0, 363.2)),
        "C"   => Some(Residue::new(TYPE_RNA,     9,  12, 7, 3, 0, 1, 0, 323.2)),
        // DNA
        "DA"  => Some(Residue::new(TYPE_DNA,     10, 12, 5, 5, 0, 1, 0, 331.2)),
        "DT"  => Some(Residue::new(TYPE_DNA,     10, 11, 7, 2, 0, 1, 0, 322.2)),
        "DG"  => Some(Residue::new(TYPE_DNA,     10, 12, 6, 5, 0, 1, 0, 347.2)),
        "DC"  => Some(Residue::new(TYPE_DNA,     9,  12, 6, 3, 0, 1, 0, 307.2)),
        _ => None,
    }
}

/// Look up a residue by 1-letter code and type.
pub fn get_residue_by_one_letter(one_letter: &str, residue_type: u8) -> Option<Residue> {
    let id = one_letter.to_uppercase();
    let three = match residue_type {
        TYPE_PROTEIN => match id.as_str() {
            "A" => "ALA",
            "B" | "N" => "ASN",
            "C" => "CYS",
            "D" => "ASP",
            "E" | "Z" => "GLU",
            "F" => "PHE",
            "G" => "GLY",
            "H" => "HIS",
            "I" => "ILE",
            "J" | "L" => "LEU",
            "K" => "LYS",
            "M" => "MET",
            "P" => "PRO",
            "Q" => "GLN",
            "R" => "ARG",
            "S" => "SER",
            "T" => "THR",
            "V" => "VAL",
            "W" => "TRP",
            "X" => "ALA", // unknown → alanine
            "Y" => "TYR",
            _ => return None,
        },
        TYPE_RNA => match id.as_str() {
            "A" => "A",
            "U" => "U",
            "G" => "G",
            "C" => "C",
            _ => return None,
        },
        TYPE_DNA => match id.as_str() {
            "A" | "M" | "V" => "DA",
            "B" | "H" | "T" | "W" => "DT",
            "D" | "G" | "K" | "R" => "DG",
            "C" | "N" | "S" => "DC",
            _ => return None,
        },
        _ => return None,
    };
    get_residue_by_three_letter(three)
}
