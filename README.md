# RADDOSE-3D Differential Fuzzer

Finds cases where the Java and Rust implementations of RADDOSE-3D produce different outputs or crash differently. Two modes:

- **`fuzz.py`** — parallel batch fuzzer, good for overnight runs
- **`test_differential.py`** — Hypothesis-based tests with automatic shrinking

## Setup

Requires Python 3.14+ and [uv](https://docs.astral.sh/uv/).

```bash
cd fuzz/
uv sync          # installs hypothesis + pytest into .venv
```

Both binaries must be built before running:

```bash
# From repo root
ant jar                             # builds raddose3d.jar
cd raddose3d && cargo build --release   # builds raddose3d/target/release/raddose3d
```

## Quick start

```bash
# 50 iterations, 4 parallel workers
uv run python fuzz.py -n 50 --workers 4

# Hypothesis tests (finds + shrinks failures automatically)
uv run pytest test_differential.py -v
```

## Batch fuzzer (`fuzz.py`)

```
uv run python fuzz.py [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-n`, `--iterations` | 100 | Number of inputs to test |
| `--strategy` | `both` | `grammar`, `mutate`, or `both` (50/50) |
| `--workers`, `-j` | 4 | Parallel workers. Each runs one Java+Rust pair. Can exceed CPU count since workers spend most time waiting on subprocesses. |
| `--budget` | 2.0 | Max cost relative to the insulin fixture (~24s Java). Controls how large/complex generated inputs can be. |
| `--timeout` | 60.0 | Per-run wall-clock timeout in seconds. Applies to each Java or Rust process individually. |
| `--seed` | random | RNG seed. Inputs are generated deterministically before dispatch, so the same seed always produces the same inputs regardless of worker count. |
| `--save-all` | off | Save every input, not just interesting ones. |
| `--minor-diff-threshold` | 0.01 | Relative difference threshold for `MAJOR_DIFF` (default 1%). |
| `--java-jar` | `../raddose3d.jar` | Path to Java jar. |
| `--rust-bin` | `../raddose3d/target/release/raddose3d` | Path to Rust binary. |

### Input strategies

**`grammar`** — Generates structurally valid inputs by sampling from a model of the full input grammar. Controls which features are on/off: subprogram (`MONTECARLO`, `XFEL`, `EMSP`), coefcalc mode (`RD3D`, `SMALLMOLE`, `CIF`, `SAXSseq`, `MicroED`), crystal type, beam type, DDM model, container, escape flags.

**`mutate`** — Perturbs numeric values in existing seed files from `corpus/seeds/`. Uses log-normal noise so values stay in the same order of magnitude. Does not mutate keywords or structure.

**`both`** — Interleaves the two strategies 50/50.

### Cost model and budget

The budget prevents generating inputs that would time out. Costs are normalized so the insulin fixture ≈ 1.0 (~12s Java). Key calibration points:

| Mode | Cost formula |
|------|-------------|
| Standard | `voxels × angle_steps / 22.5M` |
| MONTECARLO | `(base + runs × sim_electrons × 31) / 22.5M` |
| XFEL | `voxels × exposure_time × 0.178` |
| EMSP/MicroED | `base × 2 / 22.5M` |

The timeout still applies as a safety net for miscalibrated inputs, hangs, or performance divergences.

### Output

Progress is printed per completed iteration (out-of-order due to parallelism):

```
[   1/100]! MAJOR_DIFF  max_diff=7.99e-02 [Diffraction Efficiency]  java=3.8s rust=3.8s  -> ...txt
[   2/100]  MATCH  max_diff=9.83e-05 [Average DWD]  java=5.2s rust=7.1s
```

`!` marks interesting cases. All runs are logged to `results/run_<timestamp>.jsonl`. Interesting inputs are saved under `corpus/`:

```
corpus/
├── seeds/          # Seed inputs for mutation (copied from fixtures)
├── diffs/          # MAJOR_DIFF and NAN_INF
├── crashes/
│   ├── java/       # Java crashed, Rust succeeded
│   ├── rust/       # Rust crashed, Java succeeded
│   └── both/       # Both crashed
├── perf_diverge/   # One timed out, the other succeeded
└── timeouts/       # Both timed out (or one timed out, other also failed)
```

Each saved input has a `.txt` (the input file) and `.json` (metadata: category, diff values, exit codes, stderr snippets).

### Result categories

| Category | Meaning |
|----------|---------|
| `MATCH` | All metrics agree within 0.01% |
| `MINOR_DIFF` | Max relative diff between 0.01% and 1% |
| `MAJOR_DIFF` | Max relative diff ≥ 1% |
| `NAN_INF` | NaN or Inf in one output but not the other |
| `JAVA_CRASH` | Java non-zero exit, Rust succeeded |
| `RUST_CRASH` | Rust non-zero exit, Java succeeded |
| `BOTH_CRASH` | Both non-zero exit |
| `JAVA_TIMEOUT` | Java timed out, Rust also timed out or errored |
| `RUST_TIMEOUT` | Rust timed out, Java also timed out or errored |
| `BOTH_TIMEOUT` | Both timed out — usually means bad budget calibration |
| `PERF_DIVERGE` | One timed out, the other succeeded |
| `PARSE_ERROR` | Output files missing or unparseable |
| `HARNESS_ERROR` | Binary not found or failed to launch |

The comparison is done on `Summary.csv` metrics: Average DWD, Last DWD, Elastic Yield, Diffraction Efficiency, AD-WC, Max Dose, Wedge Absorbed Energy, Dose Inefficiency, Dose Inefficiency PE.

## Hypothesis tests (`test_differential.py`)

Three property-based tests, each targeting a different failure mode. When a test fails, Hypothesis automatically shrinks the input to the smallest `Config` that still triggers the failure, then saves it to `.hypothesis/` for automatic replay on future runs.

```bash
# Run all three tests (100 examples each by default)
uv run pytest test_differential.py -v

# Run more examples
RADDOSE_FUZZ_EXAMPLES=500 uv run pytest test_differential.py -v

# Fix the RNG seed (deterministic, but Hypothesis still shrinks)
uv run pytest test_differential.py -v --hypothesis-seed=42

# Run a single test
uv run pytest test_differential.py::test_no_crashes -v

# Replay previously failing examples from .hypothesis/ database
uv run pytest test_differential.py -v    # failing examples always re-run first
```

| Environment variable | Default | Description |
|----------------------|---------|-------------|
| `RADDOSE_FUZZ_EXAMPLES` | 100 | Examples per test |
| `RADDOSE_FUZZ_TIMEOUT` | 60.0 | Per-run timeout (seconds) |
| `RADDOSE_FUZZ_BUDGET` | 2.0 | Cost budget |

### Tests

**`test_outputs_match`** — Fails on `MAJOR_DIFF`, `NAN_INF`, and all crash categories. The broadest test; use this for general differential coverage.

**`test_no_crashes`** — Fails only on crash categories (`JAVA_CRASH`, `RUST_CRASH`, `BOTH_CRASH`). Use when investigating crashes in isolation, so Hypothesis shrinks toward the crash rather than the diff.

**`test_no_nan_inf`** — Fails only on `NAN_INF`. Use when investigating divergent floating-point edge cases.

### Replaying a Hypothesis failure

When a test fails, Hypothesis prints the exact `Config` dataclass values. You can reproduce the run manually:

```python
# In a Python shell from fuzz/
from generate import Config, render
from harness import run_both
from pathlib import Path
import tempfile

cfg = Config(
    crystal_type='Cuboid', dim_x=5.0, dim_y=5.0, dim_z=5.0,
    pixels_per_micron=0.1, coefcalc='MicroED',
    # ... paste full Config from Hypothesis output
)

with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp) / 'input.txt'
    p.write_text(render(cfg))
    print(render(cfg))
    java_r, rust_r = run_both(p, Path(tmp))
    print('Java exit:', java_r.exit_code, 'stderr:', java_r.stderr[:200])
    print('Rust exit:', rust_r.exit_code, 'stderr:', rust_r.stderr[:200])
```

## Coverage-driven corpus distillation (`coverage.py`)

Finds the smallest subset of inputs that together maximise Rust source-level
code coverage. Run this as a post-processing step after collecting a pool of
inputs with the fuzzer.

### Setup

Build the Rust binary with coverage instrumentation (once, separate from the
normal release binary):

```bash
cd raddose3d
RUSTFLAGS="-C instrument-coverage" cargo build --release
```

Install the LLVM tools if not already present:

```bash
rustup component add llvm-tools-preview
```

### Usage

```bash
# Analyse an entire corpus directory (recursively finds all .txt files)
uv run python coverage.py corpus/

# Analyse specific subdirectories or files
uv run python coverage.py corpus/seeds/ corpus/diffs/

# Copy the selected subset to a new directory and write a JSON report
uv run python coverage.py corpus/ --out-dir corpus/coverage-selected/ --report coverage.json
```

| Option | Default | Description |
|--------|---------|-------------|
| `--rust-bin` | `../raddose3d/target/release/raddose3d` | Path to the **instrumented** binary |
| `--timeout` | `0` (none) | Per-input timeout in seconds. **Leave at 0**: LLVM only writes coverage data on a clean process exit, so killed processes produce no data. Only set a timeout if you have genuinely unbounded inputs. |
| `--out-dir` | — | Copy selected inputs here (numbered by selection order) |
| `--report` | — | Write a JSON report with per-file line/region coverage |
| `--keep-profraws` | off | Retain raw `.profraw` files for debugging |
| `--workers`, `-j` | 1 | Parallel workers. Each uses an isolated profraw file. |
| `--patch-sim-electrons` | `1000` | Override `SimElectrons`/`SimPhotons` to N before running. Reduces MONTECARLO/XFEL runtime from hours to seconds without changing which code paths are exercised. Non-MC/XFEL inputs are unaffected. Set to `0` to disable. |

### Algorithm

For each input, the script runs the instrumented binary with `LLVM_PROFILE_FILE`
set to capture a `.profraw` file, then calls `llvm-cov export` to extract the
set of covered source regions (skipping registry crates and compiler internals).
A greedy set-cover pass then selects inputs in order of marginal coverage gain
until no new regions can be added.

### Output

Progress is printed per input, followed by a set-cover summary:

```
Phase 1: collecting coverage (1 worker(s), timeout=0.0s each)
------------------------------------------------------------
[   1/15] OK  10.0s  standard_cuboid.txt
[   2/15] OK   0.4s  ddm_bfactor.txt
...

Phase 2: greedy set-cover
------------------------------------------------------------
Total regions reachable across all inputs: 6,761
Inputs with usable coverage:               12
Selected by greedy set-cover:              9 input(s)
Regions covered by selection:              6,761 (100.0%)

Selected inputs (in selection order):
    1.  + 4,512 regions  corpus/seeds/standard_cylinder.txt
    2.  + 1,300 regions  corpus/seeds/microed_basic.txt
    ...
```

The JSON report (if requested) includes per-file line and region coverage
percentages from the merged profile of all selected inputs.

### Typical workflow

```bash
# 1. Collect a large pool of inputs (saves every input, not just failures)
uv run python fuzz.py -n 500 --save-all --workers 8

# 2. Run coverage distillation on the collected pool
uv run python coverage.py corpus/ --out-dir corpus/coverage-selected/ --report coverage.json
```

## Files

| File | Purpose |
|------|---------|
| `generate.py` | `Config` dataclass, `GrammarGenerator`, `mutate_text`, `estimate_cost`, `render` |
| `strategies.py` | Hypothesis composite strategy wrapping `GrammarGenerator` logic |
| `harness.py` | Launches both binaries in parallel, collects `RunResult` |
| `compare.py` | Parses `Summary.csv`, returns classified `Comparison` |
| `fuzz.py` | CLI entry point, `ThreadPoolExecutor` loop, corpus saving |
| `coverage.py` | Offline coverage-driven corpus distillation |
| `test_differential.py` | Hypothesis-based pytest tests |
| `corpus/seeds/` | Seed inputs (copied from `raddose3d/tests/fixtures/*.txt`) |
| `results/` | Per-run `.jsonl` logs and `_stats.json` summaries |
