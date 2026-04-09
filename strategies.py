"""
Hypothesis strategies for generating RADDOSE-3D Config objects.

Mirrors GrammarGenerator logic but expressed as Hypothesis composites so
that Hypothesis can shrink failing examples automatically.
"""

from hypothesis import assume
from hypothesis import strategies as st

from generate import (
    BeamConfig, WedgeConfig, Config, estimate_cost, DEFAULT_BUDGET,
    FIXTURES_DIR, INSULIN_BASE_COST, MC_COST_PER_ELECTRON, XFEL_PER_VOXEL_PER_SECOND,
    _pe_cost_factor,
    SMALL_MOLE_ELEMENT_POOL, HEAVY_PROTEIN_ELEMENT_POOL, SOLVENT_HEAVY_ELEMENT_POOL,
)

# ---------------------------------------------------------------------------
# Leaf-level strategies
# ---------------------------------------------------------------------------

@st.composite
def _element_counts_st(
    draw,
    pool: list[str],
    n_lo: int,
    n_hi: int,
    count_lo: float,
    count_hi: float,
    integer_counts: bool = True,
) -> str:
    """Draw a 'El count El count …' string with unique elements from pool."""
    n = draw(st.integers(n_lo, min(n_hi, len(pool))))
    elements = draw(st.lists(st.sampled_from(pool), min_size=n, max_size=n, unique=True))
    parts = []
    for el in elements:
        if integer_counts:
            count: float = draw(st.integers(int(count_lo), int(count_hi)))
        else:
            count = draw(st.floats(count_lo, count_hi, allow_nan=False, allow_infinity=False))
            count = round(count, 3)
        parts.append(f"{el} {count}")
    return " ".join(parts)


_ELEMENTS_WITH_COUNTS = _element_counts_st(
    HEAVY_PROTEIN_ELEMENT_POOL, n_lo=1, n_hi=3, count_lo=0.1, count_hi=10.0, integer_counts=False
)

_SOLVENT_HEAVY = st.one_of(
    st.just(""),
    _element_counts_st(SOLVENT_HEAVY_ELEMENT_POOL, n_lo=1, n_hi=3, count_lo=1, count_hi=2000, integer_counts=True),
)

_SMALL_MOLE_ATOMS = _element_counts_st(
    SMALL_MOLE_ELEMENT_POOL, n_lo=1, n_hi=4, count_lo=1, count_hi=12, integer_counts=True
)

_CIF_CHOICES = st.sampled_from([
    "Fe3O4", "CaCO3", "NaCl",
    str(FIXTURES_DIR / "alanine.cif"),
])

_FASTA_CHOICES = st.sampled_from([
    str(FIXTURES_DIR / "rcsb_pdb_4OR0.fasta"),
    str(FIXTURES_DIR / "insulin_seq.fasta"),
    str(FIXTURES_DIR / "bsa_fragment.fasta"),
])


# ---------------------------------------------------------------------------
# Main composite strategy
# ---------------------------------------------------------------------------

@st.composite
def raddose_config(draw, budget: float = DEFAULT_BUDGET) -> Config:
    """
    Draw a valid Config with estimated_cost <= budget.
    Uses assume() to discard over-budget examples so Hypothesis learns
    to avoid them via shrinking.
    """
    cfg = Config()

    # ---- Subprogram (controls major simulation mode) ----
    cfg.subprogram = draw(st.sampled_from(
        ["", "", "", "", "MONTECARLO", "XFEL", "EMSP"]
    ))

    # ---- Crystal type ----
    if cfg.subprogram == "EMSP":
        cfg.crystal_type = "Cuboid"
    else:
        cfg.crystal_type = draw(st.sampled_from(
            ["Cuboid", "Cuboid", "Cuboid", "Cylinder", "Spherical", "Polyhedron"]
        ))

    # ---- CoefCalc ----
    if cfg.subprogram == "EMSP":
        cfg.coefcalc = "MicroED"
    elif cfg.crystal_type == "Cylinder":
        cfg.coefcalc = "SAXSseq"
    elif cfg.crystal_type == "Polyhedron":
        cfg.coefcalc = "RD3D"
    else:
        cfg.coefcalc = draw(st.sampled_from(["RD3D", "RD3D", "RD3D", "SMALLMOLE", "CIF"]))

    # ---- Crystal dimensions ----
    if cfg.crystal_type == "Cuboid":
        cfg.dim_x = draw(st.floats(5, 300, allow_nan=False, allow_infinity=False))
        cfg.dim_y = draw(st.floats(5, 300, allow_nan=False, allow_infinity=False))
        cfg.dim_z = draw(st.floats(5, 300, allow_nan=False, allow_infinity=False))
    elif cfg.crystal_type == "Cylinder":
        cfg.dim_x = draw(st.floats(200, 3000, allow_nan=False, allow_infinity=False))
        cfg.dim_y = draw(st.floats(100, 2000, allow_nan=False, allow_infinity=False))
        cfg.dim_z = cfg.dim_y
    elif cfg.crystal_type == "Spherical":
        d = draw(st.floats(5, 300, allow_nan=False, allow_infinity=False))
        cfg.dim_x = cfg.dim_y = cfg.dim_z = d
    else:  # Polyhedron — geometry from cube.obj, dims ignored
        cfg.dim_x = cfg.dim_y = cfg.dim_z = 30.0

    # ---- PixelsPerMicron ----
    if cfg.crystal_type == "Cylinder":
        cfg.pixels_per_micron = draw(st.sampled_from([0.005, 0.01, 0.02, 0.05]))
    elif cfg.coefcalc == "SMALLMOLE":
        cfg.pixels_per_micron = draw(st.sampled_from([1.0, 2.0, 5.0]))
    else:
        cfg.pixels_per_micron = draw(st.sampled_from([0.1, 0.2, 0.5, 1.0, 2.0]))

    # ---- Composition ----
    if cfg.coefcalc in ("RD3D", "MicroED"):
        cfg.unit_cell_a = draw(st.floats(20, 200, allow_nan=False, allow_infinity=False))
        cfg.unit_cell_b = draw(st.floats(20, 200, allow_nan=False, allow_infinity=False))
        cfg.unit_cell_c = draw(st.floats(20, 200, allow_nan=False, allow_infinity=False))
        cfg.num_monomers = draw(st.integers(1, 48))
        cfg.num_residues = draw(st.integers(20, 500))
        cfg.num_rna = draw(st.integers(0, 5)) if draw(st.booleans()) else 0
        cfg.num_dna = draw(st.integers(0, 5)) if draw(st.booleans()) else 0
        cfg.solvent_fraction = draw(st.floats(0.3, 0.85, allow_nan=False, allow_infinity=False))
        cfg.heavy_protein_atoms = draw(st.one_of(st.just(""), _ELEMENTS_WITH_COUNTS))
        cfg.solvent_heavy_conc = draw(_SOLVENT_HEAVY)
    elif cfg.coefcalc == "SMALLMOLE":
        cfg.unit_cell_a = draw(st.floats(5, 50, allow_nan=False, allow_infinity=False))
        cfg.unit_cell_b = draw(st.floats(5, 50, allow_nan=False, allow_infinity=False))
        cfg.unit_cell_c = draw(st.floats(5, 50, allow_nan=False, allow_infinity=False))
        cfg.small_mole_atoms = draw(_SMALL_MOLE_ATOMS)
        cfg.num_monomers = draw(st.integers(1, 16))
    elif cfg.coefcalc == "CIF":
        cfg.cif = draw(_CIF_CHOICES)
    elif cfg.coefcalc == "SAXSseq":
        cfg.seq_file = draw(_FASTA_CHOICES)
        cfg.protein_conc = draw(st.floats(0.5, 50.0, allow_nan=False, allow_infinity=False))
        cfg.saxs_container = draw(st.booleans())

    # ---- Dose decay model ----
    cfg.ddm = draw(st.sampled_from(["Simple", "Simple", "Linear", "Leal", "Bfactor"]))
    if cfg.ddm in ("Leal", "Bfactor"):
        cfg.gamma_param = draw(st.floats(0.1, 2.0, allow_nan=False, allow_infinity=False))
        cfg.b0_param = draw(st.floats(0.1, 10.0, allow_nan=False, allow_infinity=False))
        cfg.beta_param = draw(st.floats(0.01, 2.0, allow_nan=False, allow_infinity=False))

    # ---- Subprogram parameters ----
    if cfg.subprogram == "MONTECARLO":
        cfg.runs = draw(st.integers(1, 3))
        # Decide escape flags first — they affect the cost factor
        cfg.calculate_pe_escape = draw(st.booleans())
        cfg.calculate_fl_escape = draw(st.booleans())
        pe_factor = _pe_cost_factor(cfg.calculate_pe_escape or cfg.calculate_fl_escape)
        max_electrons = int(budget * INSULIN_BASE_COST / MC_COST_PER_ELECTRON / cfg.runs / pe_factor)
        cfg.sim_electrons = draw(st.integers(10_000, min(500_000, max(10_000, max_electrons))))
    elif cfg.subprogram == "XFEL":
        cfg.runs = draw(st.integers(1, 3))

    # ---- Helper: draw a BeamConfig appropriate for cfg ----
    def _draw_beam() -> BeamConfig:
        beam = BeamConfig()
        if cfg.subprogram == "EMSP":
            beam.energy = float(draw(st.sampled_from([100, 200, 300, 400])))
            beam.beam_type = "Gaussian"
            beam.flux = draw(st.floats(1e5, 1e8, allow_nan=False, allow_infinity=False))
        else:
            beam.energy = draw(st.floats(5.0, 25.0, allow_nan=False, allow_infinity=False))
            beam.flux = draw(st.floats(1e10, 1e14, allow_nan=False, allow_infinity=False))
            beam.beam_type = draw(st.sampled_from(["Gaussian", "Tophat"]))

        max_dim = max(cfg.dim_x, cfg.dim_y)
        if cfg.crystal_type == "Cylinder":
            hw = draw(st.floats(100, min(cfg.dim_y, 1000), allow_nan=False, allow_infinity=False))
            beam.fwhm_x = beam.fwhm_y = hw
            beam.collimation_x = draw(st.floats(hw, cfg.dim_y * 0.9, allow_nan=False, allow_infinity=False))
            beam.collimation_y = beam.collimation_x
        else:
            beam.fwhm_x = draw(st.floats(5, max(10, max_dim * 1.5), allow_nan=False, allow_infinity=False))
            beam.fwhm_y = draw(st.floats(5, max(10, max_dim * 1.5), allow_nan=False, allow_infinity=False))
            beam.collimation_x = draw(st.floats(5, max(10, max_dim * 2), allow_nan=False, allow_infinity=False))
            beam.collimation_y = draw(st.floats(5, max(10, max_dim * 2), allow_nan=False, allow_infinity=False))

        beam.collimation_type = draw(st.sampled_from(["Rectangular", "Circular"]))

        if cfg.subprogram == "XFEL":
            beam.pulse_energy = draw(st.floats(1e-9, 1e-3, allow_nan=False, allow_infinity=False))

        return beam

    # ---- Segments (Beam+Wedge pairs) ----
    if cfg.subprogram in ("XFEL", "EMSP") or cfg.crystal_type == "Cylinder":
        # Single segment, single wedge with end=0 (single-image or SAXS)
        beam = _draw_beam()
        if cfg.subprogram == "XFEL":
            ppm3 = cfg.pixels_per_micron ** 3
            vox = min(cfg.dim_x * cfg.dim_y * cfg.dim_z * ppm3, 1_000_000)
            max_exp = budget / max(vox * XFEL_PER_VOXEL_PER_SECOND * cfg.runs, 1e-9)
            max_exp = max(0.0001, min(max_exp, 0.01))
            exposure_time = draw(st.floats(0.0001, max_exp, allow_nan=False, allow_infinity=False))
        elif cfg.subprogram == "EMSP":
            exposure_time = 1.0
        else:
            exposure_time = draw(st.floats(1.0, 200.0, allow_nan=False, allow_infinity=False))
        wedge = WedgeConfig(start=0.0, end=0.0, angular_resolution=1.0, exposure_time=exposure_time)
        cfg.segments = [(beam, [wedge])]

    elif cfg.subprogram == "MONTECARLO":
        # Single segment for MC
        beam = _draw_beam()
        wedge_end = float(draw(st.sampled_from([0, 45, 90, 180, 360])))
        res = float(draw(st.sampled_from([0.5, 1.0, 2.0, 5.0, 10.0]))) if wedge_end > 0 else 1.0
        wedge = WedgeConfig(
            start=0.0, end=wedge_end, angular_resolution=res,
            exposure_time=draw(st.floats(1.0, 200.0, allow_nan=False, allow_infinity=False)),
        )
        cfg.segments = [(beam, [wedge])]

    else:
        # Standard mode: 1–3 beam segments with consecutive angle ranges.
        n_segs = draw(st.integers(min_value=1, max_value=3))
        total_end = float(draw(st.sampled_from([0, 45, 90, 180, 360])))
        res = float(draw(st.sampled_from([0.5, 1.0, 2.0, 5.0, 10.0])))

        if n_segs == 1 or total_end == 0:
            n_segs = 1
            boundaries = [0.0, total_end]
        else:
            pts = sorted([
                draw(st.floats(0.0, total_end, allow_nan=False, allow_infinity=False))
                for _ in range(n_segs - 1)
            ])
            boundaries = [0.0] + pts + [total_end]

        cfg.segments = []
        for j in range(n_segs):
            seg_start = boundaries[j]
            seg_end = boundaries[j + 1]
            beam = _draw_beam()
            wedge = WedgeConfig(
                start=seg_start, end=seg_end,
                angular_resolution=res if seg_end > seg_start else 1.0,
                exposure_time=draw(st.floats(1.0, 200.0, allow_nan=False, allow_infinity=False)),
            )
            cfg.segments.append((beam, [wedge]))

    # ---- Enforce budget ----
    assume(estimate_cost(cfg) <= budget)

    return cfg
