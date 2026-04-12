# CLAUDE.md — RADDOSE-3D Differential Fuzzer

Finds cases where the Java and Rust implementations of RADDOSE-3D produce different numerical outputs or crash differently. Lives at `fuzz/` inside the Rust repo root. The Rust source is at `../` and the Java source is at `../java/`.

## Quick Start

```bash
cd fuzz
uv run python fuzz.py              # Batch fuzzer (grammar-based, runs until Ctrl-C)
uv run python fuzz.py --strategy mutate   # Mutate seeds from corpus/seeds/
uv run python test_differential.py        # Single-shot differential test run
```

Python dependencies are managed with `uv` (see `pyproject.toml`).

## Scripts

| Script | Purpose |
|--------|---------|
| `fuzz.py` | Parallel batch fuzzer — generates/mutates inputs, compares Java vs Rust, saves diffs |
| `test_differential.py` | Pytest-based differential test suite (good for CI) |
| `harness.py` | Low-level runner: invokes Java and Rust on a single input, returns structured results |
| `compare.py` | Numerical comparison logic — parses Summary.txt, computes relative differences |
| `generate.py` | Grammar-based input generator — produces random but valid RADDOSE-3D inputs |
| `strategies.py` | Fuzzing strategies: grammar generation, seed mutation, hybrid |
| `coverage.py` | Coverage-guided seed selection helpers |
| `refuzz.py` | Re-runs saved corpus inputs to recheck known diffs after code changes |

## Corpus Layout

```
corpus/
├── seeds/          # Hand-crafted seed inputs (committed)
├── diffs/          # Inputs with significant numerical differences (Java ≠ Rust)
├── crashes/
│   ├── java/       # Java crashes, Rust succeeds
│   ├── rust/       # Rust crashes, Java succeeds
│   └── both/       # Both crash
└── results/        # Per-run JSON logs and aggregate stats
```

## How It Works

1. `generate.py` / `strategies.py` produces a RADDOSE-3D input string
2. `harness.py` runs both `java -jar ../java/raddose3d.jar` and `../target/release/raddose3d` on the input
3. `compare.py` parses both `Summary.txt` outputs and computes relative differences for key metrics (Max Dose, DWD, Absorbed Energy, etc.)
4. Inputs exceeding the tolerance threshold (default ~1%) are saved to `corpus/diffs/`
5. Crashes are saved by which side crashed

## Known Expected Differences

Some variance between Java and Rust is intentional (see root `CLAUDE.md` for full details):

- **findDepth deduplication bug** (~0.3% Max Dose on coarse cuboid grids): Java's `==` reference equality on boxed Doubles is a no-op, causing boundary voxels to get depth=0. Rust matches this behavior but FP path differences mean the exact affected voxel set varies.
- **Gaussian beam circular collimation**: Both use the same polar integration algorithm; differences here indicate a regression.
- **PE/FL escape**: ~2% random noise from single-track sampling; Max Dose can vary ±10%.
- **`Type SphericalAnalytic`**: Rust-only crystal type — exclude from differential fuzzing.

## Running Against a Specific Input

```bash
uv run python harness.py path/to/input.txt
```

Prints a side-by-side comparison of Java and Rust outputs.
