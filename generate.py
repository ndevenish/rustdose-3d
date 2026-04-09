"""
Grammar-based and mutation-based input generators for RADDOSE-3D fuzzing.

Cost model: normalized so insulin_test.txt ≈ 1.0.
GrammarGenerator retries until estimated cost <= budget.
MutationGenerator perturbs numeric fields in existing seed text.
"""

import math
import os
import random
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

FIXTURES_DIR = Path(__file__).parent.parent / "raddose3d" / "tests" / "fixtures"

# ---------------------------------------------------------------------------
# Element pools — used by both GrammarGenerator and Hypothesis strategies
# ---------------------------------------------------------------------------

# Atoms that make up a small-molecule crystal formula unit.
# Integer stoichiometric counts (1–12).
SMALL_MOLE_ELEMENT_POOL = [
    "C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I",
    "Na", "Mg", "Al", "Si", "K", "Ca", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Se", "Mo", "Ag",
]

# Anomalously scattering atoms added to a protein.
# Counts are per-monomer occupancies (fractional allowed, 0.1–10).
HEAVY_PROTEIN_ELEMENT_POOL = [
    "Se", "S", "Zn", "Fe", "Cu", "Mn", "Ca", "Co", "Ni",
    "Br", "I", "Hg", "Pt", "Au", "Mo",
]

# Heavy atoms dissolved in the solvent, given as concentration in mM (1–2000).
SOLVENT_HEAVY_ELEMENT_POOL = [
    "Na", "K", "Mg", "Ca", "Mn", "Fe", "Zn", "Cu", "Se",
    "Br", "I", "P", "Cl", "Rb", "Sr", "Cs", "Ba",
]


def _element_counts_str(
    rng: random.Random,
    pool: list[str],
    n_lo: int,
    n_hi: int,
    count_lo: float,
    count_hi: float,
    integer_counts: bool = True,
) -> str:
    """Return a 'El count El count …' string with 1–n_hi unique elements."""
    n = rng.randint(n_lo, min(n_hi, len(pool)))
    elements = rng.sample(pool, n)
    parts = []
    for el in elements:
        if integer_counts:
            count: float = rng.randint(int(count_lo), int(count_hi))
        else:
            count = round(rng.uniform(count_lo, count_hi), 3)
        parts.append(f"{el} {count}")
    return " ".join(parts)

# ---------------------------------------------------------------------------
# Cost model — normalized so insulin_test.txt ≈ 1.0 (Java ~12s)
#
# Calibration measurements (Java, aarch64):
#   insulin_test:  ~12s   → 125k voxels × 180 angle-steps = 22.5M units
#   mc_test:      ~583s   → 125k voxels × 9 steps + 1e7 electrons MC
#                           base≈0.6s, MC overhead≈582s
#                           → ~31 units per simulated electron (not 0.7)
#   xfel_test:    ~268s   → 125 voxels × 1 step, ExposureTime=1s
#                           XFEL cost ≈ voxels × ExposureTime × 4M units/s
#                           → XFEL_PER_VOXEL_PER_SECOND ≈ 0.178 (normalized)
#   microed_test:  ~1s    → tiny; MICROED_MULTIPLIER kept small
# ---------------------------------------------------------------------------
INSULIN_BASE_COST = 22_500_000.0
MC_COST_PER_ELECTRON = 31.0      # work units per simulated electron (measured)
XFEL_PER_VOXEL_PER_SECOND = 0.178  # normalized cost per voxel per second of exposure
MICROED_MULTIPLIER = 2.0
DEFAULT_BUDGET = 2.0             # ~2x insulin ≈ 24s Java max


# ---------------------------------------------------------------------------
# Per-segment dataclasses
# ---------------------------------------------------------------------------

@dataclass
class BeamConfig:
    """Parameters for a single Beam block."""
    beam_type: str = "Gaussian"     # Gaussian | Tophat
    flux: float = 2e12
    fwhm_x: float = 100.0
    fwhm_y: float = 100.0
    energy: float = 12.1
    collimation_x: float = 100.0
    collimation_y: float = 100.0
    collimation_type: str = "Rectangular"  # Rectangular | Circular
    pulse_energy: float = 2e-6             # XFEL only


@dataclass
class WedgeConfig:
    """Parameters for a single Wedge block."""
    start: float = 0.0
    end: float = 360.0
    exposure_time: float = 100.0
    angular_resolution: float = 2.0


# ---------------------------------------------------------------------------
# Top-level Config
# ---------------------------------------------------------------------------

@dataclass
class Config:
    # Crystal geometry
    crystal_type: str = "Cuboid"   # Cuboid | Cylinder | Spherical | Polyhedron
    dim_x: float = 100.0
    dim_y: float = 100.0
    dim_z: float = 100.0
    pixels_per_micron: float = 0.5

    # Absorption coefficient mode
    coefcalc: str = "RD3D"   # RD3D | SMALLMOLE | CIF | SAXSseq | MicroED

    # Standard protein composition (RD3D / MicroED)
    unit_cell_a: float = 78.0
    unit_cell_b: float = 78.0
    unit_cell_c: float = 78.0
    num_monomers: int = 24
    num_residues: int = 51
    num_rna: int = 0
    num_dna: int = 0
    heavy_protein_atoms: str = ""      # e.g. "Zn 0.333 S 6"
    solvent_heavy_conc: str = ""       # e.g. "P 425"
    solvent_fraction: float = 0.64

    # Small-molecule composition (SMALLMOLE)
    small_mole_atoms: str = "Mg O 3"

    # CIF mode
    cif: str = "Fe3O4"

    # SAXS mode
    seq_file: str = ""
    protein_conc: float = 2.0
    saxs_container: bool = False
    container_elements: str = "Si 1 O 2"
    container_thickness: float = 50.0
    container_density: float = 2.648

    # Dose decay model
    ddm: str = "Simple"        # Simple | Linear | Leal | Bfactor
    gamma_param: float = 0.5
    b0_param: float = 1.0
    beta_param: float = 0.1    # must be positive (Gumbel-sign issue)

    # Subprogram
    subprogram: str = ""        # "" | MONTECARLO | XFEL | EMSP
    runs: int = 1
    sim_electrons: int = 1_000_000
    calculate_pe_escape: bool = False
    calculate_fl_escape: bool = False

    # Beam+Wedge segments: list of (BeamConfig, [WedgeConfig, ...]) pairs.
    # One entry per Beam block; each entry may have multiple consecutive Wedge blocks.
    segments: list = field(default_factory=lambda: [
        (BeamConfig(), [WedgeConfig()])
    ])


# ---------------------------------------------------------------------------
# Cost estimation helpers
# ---------------------------------------------------------------------------

def _pe_cost_factor(escape_active: bool) -> float:
    """
    Multiplicative cost factor for PE/FL escape in MONTECARLO mode.
    When PE_ANGLE_RESOLUTION=N is set (default 1), escape calculations
    cost O(N²) because the track count scales as N² over the sphere.
    """
    if not escape_active:
        return 1.0
    pe_res = int(os.environ.get("PE_ANGLE_RESOLUTION", "1"))
    return float(pe_res * pe_res)


# ---------------------------------------------------------------------------
# Cost estimation
# ---------------------------------------------------------------------------

def estimate_cost(cfg: Config) -> float:
    """Return cost relative to insulin_test.txt (≈ 1.0)."""
    ppm3 = cfg.pixels_per_micron ** 3
    if cfg.crystal_type == "Cuboid":
        voxels = cfg.dim_x * cfg.dim_y * cfg.dim_z * ppm3
    elif cfg.crystal_type == "Cylinder":
        r = cfg.dim_y / 2.0
        voxels = math.pi * r * r * cfg.dim_x * ppm3
    elif cfg.crystal_type == "Spherical":
        r = cfg.dim_x / 2.0
        voxels = (4.0 / 3.0) * math.pi * r * r * r * ppm3
    else:  # Polyhedron
        voxels = cfg.dim_x * cfg.dim_y * cfg.dim_z * ppm3
    voxels = min(voxels, 1_000_000)  # Java auto-caps at 1M

    all_wedges = [w for _, wedges in cfg.segments for w in wedges]

    def _steps(w: WedgeConfig) -> float:
        span = abs(w.end - w.start)
        return max(1.0, span / w.angular_resolution) if span > 0 else 1.0

    if cfg.subprogram == "XFEL":
        # XFEL_PER_VOXEL_PER_SECOND is already normalized, no INSULIN_BASE_COST division.
        total = sum(voxels * w.exposure_time * XFEL_PER_VOXEL_PER_SECOND for w in all_wedges)
        return total * max(1, cfg.runs)

    total_steps = sum(_steps(w) for w in all_wedges)
    base = voxels * total_steps

    if cfg.subprogram == "MONTECARLO":
        pe_factor = _pe_cost_factor(cfg.calculate_pe_escape or cfg.calculate_fl_escape)
        mc = cfg.runs * cfg.sim_electrons * MC_COST_PER_ELECTRON * pe_factor
        return (base + mc) / INSULIN_BASE_COST
    elif cfg.subprogram == "EMSP":
        return base * MICROED_MULTIPLIER / INSULIN_BASE_COST
    else:
        return base / INSULIN_BASE_COST


# ---------------------------------------------------------------------------
# Renderer
# ---------------------------------------------------------------------------

def render(cfg: Config) -> str:
    """Render a Config to RADDOSE-3D input file text."""
    lines = ["Crystal"]
    lines.append(f"Type {cfg.crystal_type}")

    if cfg.crystal_type == "Cuboid":
        lines.append(f"Dimensions {cfg.dim_x} {cfg.dim_y} {cfg.dim_z}")
    elif cfg.crystal_type == "Cylinder":
        lines.append(f"Dimensions {cfg.dim_x} {cfg.dim_y}")
    elif cfg.crystal_type == "Spherical":
        lines.append(f"Dimensions {cfg.dim_x}")
    elif cfg.crystal_type == "Polyhedron":
        lines.append(f"Dimensions {cfg.dim_x} {cfg.dim_y} {cfg.dim_z}")
        lines.append(f"Wireframetype obj")
        lines.append(f"ModelFile {FIXTURES_DIR / 'cube.obj'}")

    lines.append(f"PixelsPerMicron {cfg.pixels_per_micron}")
    # Java's grammar token named CIF matches the text "EXPSM", not "CIF".
    coefcalc_keyword = "EXPSM" if cfg.coefcalc == "CIF" else cfg.coefcalc
    lines.append(f"AbsCoefCalc {coefcalc_keyword}")

    if cfg.coefcalc in ("RD3D", "MicroED"):
        lines.append(f"UnitCell {cfg.unit_cell_a} {cfg.unit_cell_b} {cfg.unit_cell_c}")
        lines.append(f"NumMonomers {cfg.num_monomers}")
        lines.append(f"NumResidues {cfg.num_residues}")
        if cfg.num_rna > 0:
            lines.append(f"NumRNA {cfg.num_rna}")
        if cfg.num_dna > 0:
            lines.append(f"NumDNA {cfg.num_dna}")
        if cfg.heavy_protein_atoms:
            lines.append(f"ProteinHeavyAtoms {cfg.heavy_protein_atoms}")
        if cfg.solvent_heavy_conc:
            lines.append(f"SolventHeavyConc {cfg.solvent_heavy_conc}")
        if cfg.coefcalc == "RD3D":
            lines.append(f"SolventFraction {cfg.solvent_fraction}")
    elif cfg.coefcalc == "SMALLMOLE":
        lines.append(f"UnitCell {cfg.unit_cell_a} {cfg.unit_cell_b} {cfg.unit_cell_c}")
        lines.append(f"SmallMoleAtoms {cfg.small_mole_atoms}")
        lines.append(f"NumMonomers {cfg.num_monomers}")
    elif cfg.coefcalc == "CIF":
        lines.append(f"CIF {cfg.cif}")
    elif cfg.coefcalc == "SAXSseq":
        lines.append(f"SeqFile {cfg.seq_file}")
        lines.append(f"ProteinConc {cfg.protein_conc}")
        if cfg.saxs_container:
            lines.append(f"ContainerMaterialType elemental")
            lines.append(f"MaterialElements {cfg.container_elements}")
            lines.append(f"ContainerThickness {cfg.container_thickness}")
            lines.append(f"ContainerDensity {cfg.container_density}")

    if cfg.ddm != "Simple":
        lines.append(f"DDM {cfg.ddm}")
        if cfg.ddm in ("Leal", "Bfactor"):
            lines.append(f"DecayParam {cfg.gamma_param} {cfg.b0_param} {cfg.beta_param}")

    if cfg.subprogram:
        lines.append(f"Subprogram {cfg.subprogram}")
        if cfg.subprogram == "MONTECARLO":
            lines.append(f"Runs {cfg.runs}")
            lines.append(f"SimElectrons {cfg.sim_electrons}")
            if cfg.calculate_pe_escape:
                lines.append("CalculatePEEscape TRUE")
            if cfg.calculate_fl_escape:
                lines.append("CalculateFLEscape TRUE")
        elif cfg.subprogram == "XFEL":
            lines.append(f"Runs {cfg.runs}")

    # ---- One Beam block + N Wedge blocks per segment ----
    for beam, wedges in cfg.segments:
        lines.append("")
        lines.append("Beam")
        lines.append(f"Type {beam.beam_type}")
        lines.append(f"Energy {beam.energy}")

        if cfg.subprogram != "EMSP":
            lines.append(f"Flux {beam.flux:.3e}")

        if beam.beam_type == "Gaussian":
            lines.append(f"FWHM {beam.fwhm_x} {beam.fwhm_y}")

        lines.append(f"Collimation {beam.collimation_type} {beam.collimation_x} {beam.collimation_y}")

        if cfg.subprogram == "XFEL":
            lines.append(f"PulseEnergy {beam.pulse_energy:.3e}")

        for w in wedges:
            lines.append("")
            lines.append(f"Wedge {w.start} {w.end}")
            lines.append(f"ExposureTime {w.exposure_time}")
            lines.append(f"AngularResolution {w.angular_resolution}")

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Grammar-based generator
# ---------------------------------------------------------------------------

class GrammarGenerator:
    """
    Samples random valid configs from the grammar.
    Retries until estimated_cost(config) <= budget.
    Feature flags (subprogram, PE/FL escape, container, DDM, etc.) are
    sampled categorically so all combinations are reachable.
    """

    def __init__(self, budget: float = DEFAULT_BUDGET, rng: Optional[random.Random] = None):
        self.budget = budget
        self.rng = rng or random.Random()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _u(self, lo: float, hi: float) -> float:
        return self.rng.uniform(lo, hi)

    def _i(self, lo: int, hi: int) -> int:
        return self.rng.randint(lo, hi)

    def _c(self, seq):
        return self.rng.choice(seq)

    def _flip(self, p: float = 0.5) -> bool:
        return self.rng.random() < p

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate(self, max_retries: int = 200) -> Config:
        for _ in range(max_retries):
            cfg = self._sample()
            if estimate_cost(cfg) <= self.budget:
                return cfg
        return self._minimal_fallback()

    # ------------------------------------------------------------------
    # Sampling
    # ------------------------------------------------------------------

    def _sample(self) -> Config:
        cfg = Config()

        # ---- Subprogram (controls major simulation mode) ----
        cfg.subprogram = self._c(["", "", "MONTECARLO", "XFEL", "EMSP"])

        # ---- Crystal type ----
        if cfg.subprogram == "EMSP":
            cfg.crystal_type = "Cuboid"
        else:
            cfg.crystal_type = self._c(["Cuboid", "Cuboid",
                                         "Cylinder", "Spherical", "Polyhedron"])

        # ---- CoefCalc (must match crystal type and subprogram) ----
        if cfg.subprogram == "EMSP":
            cfg.coefcalc = "MicroED"
        elif cfg.crystal_type == "Cylinder":
            cfg.coefcalc = self._c(["SAXSseq", "SAXSseq", "RD3D", "CIF"])
        elif cfg.crystal_type == "Polyhedron":
            cfg.coefcalc = self._c(["RD3D", "RD3D", "SMALLMOLE", "CIF"])
        else:
            cfg.coefcalc = self._c(["RD3D", "RD3D", "SMALLMOLE", "CIF"])

        # ---- Crystal dimensions ----
        if cfg.crystal_type == "Cuboid":
            cfg.dim_x = round(self._u(0.5, 2000), 2)
            cfg.dim_y = round(self._u(0.5, 2000), 2)
            cfg.dim_z = round(self._u(0.5, 2000), 2)
        elif cfg.crystal_type == "Cylinder":
            cfg.dim_x = round(self._u(10, 10000), 1)   # height
            cfg.dim_y = round(self._u(5, 5000), 1)     # diameter
            cfg.dim_z = cfg.dim_y
        elif cfg.crystal_type == "Spherical":
            d = round(self._u(0.5, 2000), 2)
            cfg.dim_x = cfg.dim_y = cfg.dim_z = d
        else:  # Polyhedron — cube.obj defines geometry; dims set bounding box for cost model
            d = round(self._u(5, 500), 2)
            cfg.dim_x = cfg.dim_y = cfg.dim_z = d

        # ---- PixelsPerMicron (major cost driver — budget enforced by retry loop) ----
        cfg.pixels_per_micron = round(self._u(0.001, 10.0), 4)

        # ---- Composition ----
        if cfg.coefcalc in ("RD3D", "MicroED"):
            cfg.unit_cell_a = round(self._u(3, 800), 2)
            cfg.unit_cell_b = round(self._u(3, 800), 2)
            cfg.unit_cell_c = round(self._u(3, 800), 2)
            cfg.num_monomers = self._i(1, 200)
            cfg.num_residues = self._i(1, 5000)
            cfg.num_rna = self._i(0, 5) if self._flip(0.1) else 0
            cfg.num_dna = self._i(0, 5) if self._flip(0.1) else 0
            cfg.solvent_fraction = round(self._u(0.01, 0.99), 4)
            cfg.heavy_protein_atoms = (
                _element_counts_str(self.rng, HEAVY_PROTEIN_ELEMENT_POOL, 1, 3, 0.1, 10.0, integer_counts=False)
                if self._flip(0.5) else ""
            )
            cfg.solvent_heavy_conc = (
                _element_counts_str(self.rng, SOLVENT_HEAVY_ELEMENT_POOL, 1, 3, 1, 2000, integer_counts=True)
                if self._flip(0.5) else ""
            )
        elif cfg.coefcalc == "SMALLMOLE":
            cfg.unit_cell_a = round(self._u(5, 50), 3)
            cfg.unit_cell_b = round(self._u(5, 50), 3)
            cfg.unit_cell_c = round(self._u(5, 50), 3)
            cfg.small_mole_atoms = _element_counts_str(
                self.rng, SMALL_MOLE_ELEMENT_POOL, 1, 4, 1, 12, integer_counts=True
            )
            cfg.num_monomers = self._i(1, 16)
        elif cfg.coefcalc == "CIF":
            cfg.cif = str(FIXTURES_DIR / "alanine.cif")
        elif cfg.coefcalc == "SAXSseq":
            cfg.seq_file = str(FIXTURES_DIR / self._c(["rcsb_pdb_4OR0.fasta",
                                                         "insulin_seq.fasta",
                                                         "bsa_fragment.fasta"]))
            cfg.protein_conc = round(self._u(0.01, 500.0), 3)
            cfg.saxs_container = self._flip(0.5)

        # ---- Dose decay model ----
        cfg.ddm = self._c(["Simple", "Linear", "Leal", "Bfactor"])
        if cfg.ddm in ("Leal", "Bfactor"):
            cfg.gamma_param = round(self._u(0.001, 100.0), 4)
            cfg.b0_param = round(self._u(0.001, 500.0), 4)
            cfg.beta_param = round(max(0.001, self._u(0.001, 50.0)), 4)  # must stay positive

        # ---- Subprogram parameters ----
        if cfg.subprogram == "MONTECARLO":
            cfg.runs = self._i(1, 3)
            # Decide escape flags first — they affect the cost factor
            cfg.calculate_pe_escape = self._flip(0.4)
            cfg.calculate_fl_escape = self._flip(0.3)
            pe_factor = _pe_cost_factor(cfg.calculate_pe_escape or cfg.calculate_fl_escape)
            # Cap electrons so MC overhead fits within budget
            max_electrons = int(
                self.budget * INSULIN_BASE_COST / MC_COST_PER_ELECTRON
                / max(1, cfg.runs) / pe_factor
            )
            cfg.sim_electrons = self._i(10_000, min(500_000, max(10_000, max_electrons)))
        elif cfg.subprogram == "XFEL":
            cfg.runs = self._i(1, 3)

        # ---- Segments (Beam+Wedge) ----
        cfg.segments = self._sample_segments(cfg)

        return cfg

    def _sample_beam(self, cfg: Config) -> BeamConfig:
        """Sample a BeamConfig appropriate for cfg's subprogram and crystal type."""
        beam = BeamConfig()

        if cfg.subprogram == "EMSP":
            beam.energy = round(self._u(20, 1000), 1)
            beam.beam_type = "Gaussian"
            beam.flux = 10 ** self._u(3, 10)
        else:
            beam.energy = round(self._u(1.0, 100.0), 3)
            beam.flux = 10 ** self._u(6, 16)
            beam.beam_type = self._c(["Gaussian", "Tophat"])

        if cfg.crystal_type == "Cylinder":
            max_fwhm = max(1.0, cfg.dim_y * 0.5)
            hw = round(self._u(max(0.5, cfg.dim_y * 0.05), max_fwhm), 1)
            beam.fwhm_x = beam.fwhm_y = hw
            beam.collimation_x = round(self._u(hw, max(hw + 1.0, cfg.dim_y * 0.8)), 1)
            beam.collimation_y = beam.collimation_x
        else:
            beam.fwhm_x = round(self._u(5, max(10, cfg.dim_x * 1.5)), 2)
            beam.fwhm_y = round(self._u(5, max(10, cfg.dim_y * 1.5)), 2)
            beam.collimation_x = round(self._u(5, max(10, cfg.dim_x * 2)), 2)
            beam.collimation_y = round(self._u(5, max(10, cfg.dim_y * 2)), 2)

        beam.collimation_type = self._c(["Rectangular", "Circular"])

        if cfg.subprogram == "XFEL":
            beam.pulse_energy = 10 ** self._u(-9, -3)

        return beam

    def _sample_segments(self, cfg: Config) -> list:
        """
        Sample Beam+Wedge segments.

        Special modes (XFEL, EMSP, Cylinder, MONTECARLO): single segment, single wedge.
        Standard mode: 1–3 segments with consecutive angular ranges.

        Returns list of (BeamConfig, [WedgeConfig]) tuples.
        """
        ppm3 = cfg.pixels_per_micron ** 3
        vox = min(cfg.dim_x * cfg.dim_y * cfg.dim_z * ppm3, 1_000_000)

        if cfg.subprogram == "XFEL" or cfg.subprogram == "EMSP" or cfg.crystal_type == "Cylinder":
            beam = self._sample_beam(cfg)
            if cfg.subprogram == "XFEL":
                max_exp = self.budget / max(vox * XFEL_PER_VOXEL_PER_SECOND * cfg.runs, 1e-9)
                max_exp = max(0.0001, min(max_exp, 0.01))
                exposure_time = round(self._u(0.0001, max_exp), 6)
            elif cfg.subprogram == "EMSP":
                exposure_time = 1.0  # EMSP doesn't use exposure time, placeholder
            else:
                exposure_time = round(self._u(1.0, 200.0), 2)
            wedge = WedgeConfig(start=0.0, end=0.0, angular_resolution=1.0,
                                exposure_time=exposure_time)
            return [(beam, [wedge])]

        if cfg.subprogram == "MONTECARLO":
            beam = self._sample_beam(cfg)
            wedge_end = float(self._c([0, 45, 90, 180, 360]))
            res = round(self._u(0.1, 90.0), 2) if wedge_end > 0 else 1.0
            wedge = WedgeConfig(
                start=0.0, end=wedge_end, angular_resolution=res,
                exposure_time=round(self._u(1.0, 200.0), 2),
            )
            return [(beam, [wedge])]

        # Standard mode: 1–3 beam segments with consecutive angle ranges.
        n_segs = self._c([1, 1, 2, 2, 3])
        total_end = float(self._c([0, 45, 90, 180, 360]))
        res = round(self._u(0.1, 90.0), 2)

        if n_segs == 1 or total_end == 0:
            boundaries = [0.0, total_end]
            n_segs = 1
        else:
            pts = sorted(round(self._u(0.0, total_end), 1) for _ in range(n_segs - 1))
            boundaries = [0.0] + pts + [total_end]

        segments = []
        for j in range(n_segs):
            seg_start = boundaries[j]
            seg_end = boundaries[j + 1]
            beam = self._sample_beam(cfg)
            wedge = WedgeConfig(
                start=seg_start, end=seg_end,
                angular_resolution=res if seg_end > seg_start else 1.0,
                exposure_time=round(self._u(1.0, 200.0), 2),
            )
            segments.append((beam, [wedge]))
        return segments

    def _minimal_fallback(self) -> Config:
        """Guaranteed-cheap config used when retries are exhausted."""
        cfg = Config(
            crystal_type="Cuboid", dim_x=50, dim_y=50, dim_z=50,
            pixels_per_micron=0.5, coefcalc="RD3D",
            unit_cell_a=78, unit_cell_b=78, unit_cell_c=78,
            num_monomers=24, num_residues=51,
        )
        cfg.segments = [(
            BeamConfig(
                beam_type="Gaussian", flux=2e12,
                fwhm_x=100, fwhm_y=100, energy=12.1,
                collimation_x=100, collimation_y=100,
            ),
            [WedgeConfig(start=0, end=90, exposure_time=10, angular_resolution=5)],
        )]
        return cfg


# ---------------------------------------------------------------------------
# Mutation-based generator
# ---------------------------------------------------------------------------

# Matches integers and floats including scientific notation.
# Excludes numbers that are part of element symbols like "Fe3O4".
_NUMBER_RE = re.compile(
    r'(?<![A-Za-z])(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(?![A-Za-z])'
)


def mutate_text(text: str, mutation_rate: float = 0.25,
                rng: Optional[random.Random] = None) -> str:
    """
    Perturb numeric values in a seed input file.
    Uses log-normal noise so values stay in the same order of magnitude.
    Returns the mutated text; does not check cost budgets.
    """
    if rng is None:
        rng = random.Random()

    def perturb(m: re.Match) -> str:
        if rng.random() > mutation_rate:
            return m.group(0)
        raw = m.group(1)
        val = float(raw)
        if val == 0.0:
            # Zero is special: flip to a small positive value occasionally
            return "0" if rng.random() < 0.7 else f"{rng.uniform(0.01, 1.0):.3g}"
        factor = math.exp(rng.gauss(0.0, 0.4))
        new_val = val * factor
        # Preserve integer appearance for small whole numbers
        if '.' not in raw and 'e' not in raw.lower() and val < 10000:
            return str(max(1, int(round(new_val))))
        # Use scientific notation when the original did, or for extreme values
        if 'e' in raw.lower() or abs(new_val) > 1e5 or abs(new_val) < 0.001:
            return f"{new_val:.3e}"
        return f"{new_val:.4g}"

    return _NUMBER_RE.sub(perturb, text)
