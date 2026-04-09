#!/usr/bin/env python3
"""
Differential fuzzer for RADDOSE-3D Java vs Rust.

Usage:
    python fuzz.py [options]

Strategies:
    grammar    -- generate inputs from the full grammar (default)
    mutate     -- mutate seed files from corpus/seeds/
    both       -- interleave grammar and mutation (50/50)

Output is written under corpus/:
    corpus/diffs/         -- inputs with significant numerical differences
    corpus/crashes/java/  -- inputs where Java crashes, Rust succeeds
    corpus/crashes/rust/  -- inputs where Rust crashes, Java succeeds
    corpus/crashes/both/  -- inputs where both crash
    results/              -- per-run JSON logs and aggregate stats
"""

import argparse
import json
import os
import random
import re
import shutil
import sys
import tempfile
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, FIRST_COMPLETED, wait
from pathlib import Path

# Add fuzz/ dir to sys.path so sibling imports work
sys.path.insert(0, str(Path(__file__).parent))

from compare import Category, compare
from generate import (
    DEFAULT_BUDGET, Config, BeamConfig, WedgeConfig, GrammarGenerator,
    estimate_cost, mutate_text, render, generate_chaos, CoverageTracker,
)
from harness import (
    DEFAULT_JAVA_JAR, DEFAULT_RUST_BIN, DEFAULT_TIMEOUT, run_both,
)

FUZZ_DIR = Path(__file__).parent
CORPUS_DIR = FUZZ_DIR / "corpus"
RESULTS_DIR = FUZZ_DIR / "results"
SEEDS_DIR = CORPUS_DIR / "seeds"


# ---------------------------------------------------------------------------
# Skip counter
# ---------------------------------------------------------------------------

class SkipCounter:
    """
    Tracks BOTH_TIMEOUT occurrences per (subprogram, coefcalc, crystal_type)
    triple and progressively tightens the effective budget for that triple.

    Each BOTH_TIMEOUT halves the budget multiplier for the triple (floor 0.0625,
    i.e. at most 4 halvings before the triple is effectively frozen out).
    The multiplier recovers by 25% per non-timeout result so the generator
    can explore again if the budget was just temporarily miscalibrated.

    Thread-safe: the fuzzer's result-handling thread and the generation thread
    are the same (generation happens in the main thread before dispatch), so
    the lock is only there for future safety.
    """

    MIN_MULTIPLIER = 0.0625   # 1/16 — four halvings

    def __init__(self):
        self._lock = threading.Lock()
        self._multipliers: dict[tuple, float] = defaultdict(lambda: 1.0)
        self._timeout_counts: Counter = Counter()

    def record_timeout(self, triple: tuple) -> None:
        with self._lock:
            cur = self._multipliers[triple]
            self._multipliers[triple] = max(self.MIN_MULTIPLIER, cur * 0.5)
            self._timeout_counts[triple] += 1

    def record_ok(self, triple: tuple) -> None:
        """Partial recovery: nudge multiplier back toward 1.0 on non-timeout."""
        with self._lock:
            cur = self._multipliers[triple]
            if cur < 1.0:
                self._multipliers[triple] = min(1.0, cur * 1.25)

    def effective_budget(self, triple: tuple, base_budget: float) -> float:
        with self._lock:
            return base_budget * self._multipliers[triple]

    def summary(self) -> list[tuple]:
        """Return list of (triple, multiplier, timeout_count) sorted by most-penalised."""
        with self._lock:
            return sorted(
                [(t, self._multipliers[t], self._timeout_counts[t])
                 for t in self._timeout_counts],
                key=lambda x: x[1],
            )


def _triple_from_text(text: str) -> tuple[str, str, str]:
    """Extract (subprogram, coefcalc, crystal_type) from raw input text."""
    sub = re.search(r'Subprogram\s+(\w+)', text, re.IGNORECASE)
    coef = re.search(r'(?:AbsCoefCalc|CoefCalc)\s+(\w+)', text, re.IGNORECASE)
    ctype = re.search(r'Type\s+(\w+)', text, re.IGNORECASE)
    return (
        sub.group(1).upper() if sub else "",
        coef.group(1).upper() if coef else "RD3D",
        ctype.group(1).capitalize() if ctype else "Cuboid",
    )


def main():
    ap = argparse.ArgumentParser(description="RADDOSE-3D differential fuzzer")
    ap.add_argument("--iterations", "-n", type=int, default=100,
                    help="Number of fuzz iterations (default: 100)")
    ap.add_argument("--strategy",
                    choices=["grammar", "mutate", "both", "structural", "chaos"],
                    default="both", help="Input generation strategy (default: both)")
    ap.add_argument("--budget", type=float, default=DEFAULT_BUDGET,
                    help=f"Max cost relative to insulin (default: {DEFAULT_BUDGET})")
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT,
                    help=f"Per-run wall-clock timeout in seconds (default: {DEFAULT_TIMEOUT})")
    ap.add_argument("--java-jar", type=Path, default=DEFAULT_JAVA_JAR)
    ap.add_argument("--rust-bin", type=Path, default=DEFAULT_RUST_BIN)
    ap.add_argument("--seed", type=int, default=None,
                    help="RNG seed for reproducibility")
    ap.add_argument("--save-all", action="store_true",
                    help="Save every input, not just interesting ones")
    ap.add_argument("--minor-diff-threshold", type=float, default=1e-2,
                    help="Relative diff threshold for MAJOR_DIFF (default: 0.01 = 1%%)")
    ap.add_argument("--workers", "-j", type=int, default=4,
                    help="Parallel workers (default: 4). Workers wait on subprocesses so "
                         "this can exceed CPU count. Each worker runs one Java+Rust pair.")
    args = ap.parse_args()

    # Override compare module threshold if requested
    import compare as compare_mod
    compare_mod.TOL_MINOR = args.minor_diff_threshold

    rng = random.Random(args.seed)
    grammar_gen = GrammarGenerator(budget=args.budget, rng=rng)
    seeds = _load_seeds()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    run_id = time.strftime("%Y%m%d_%H%M%S")
    log_path = RESULTS_DIR / f"run_{run_id}.jsonl"
    stats_path = RESULTS_DIR / f"run_{run_id}_stats.json"

    category_counts: Counter = Counter()
    interesting_saved = 0
    total_java_time = 0.0
    total_rust_time = 0.0
    skip_counter = SkipCounter()
    coverage_tracker = CoverageTracker(bias_after=50)

    print(f"RADDOSE-3D differential fuzzer  |  {args.iterations} iterations  |  "
          f"strategy={args.strategy}  |  budget={args.budget}x  |  "
          f"timeout={args.timeout}s  |  workers={args.workers}")
    print(f"Java:  {args.java_jar}")
    print(f"Rust:  {args.rust_bin}")
    print(f"Log:   {log_path}")
    print()

    def _run_one(item):
        i, source, input_text, cfg_cost, coverage_key = item
        with tempfile.TemporaryDirectory(prefix="raddose_fuzz_") as tmp:
            tmp_path = Path(tmp)
            input_path = tmp_path / "input.txt"
            input_path.write_text(input_text)
            java_r, rust_r = run_both(
                input_path, tmp_path,
                java_jar=args.java_jar,
                rust_bin=args.rust_bin,
                timeout=args.timeout,
            )
            result = compare(java_r, rust_r)
        return i, result, java_r, rust_r, source, cfg_cost, input_text, coverage_key

    def _next_item(i: int) -> tuple:
        """Generate the i-th work item, using skip-counter-adjusted budget."""
        source, input_text, cfg_cost = _generate(
            args.strategy, grammar_gen, seeds, rng, args.budget,
            skip_counter=skip_counter,
            coverage_tracker=coverage_tracker,
        )
        coverage_key = _coverage_key_from_text(input_text)
        return (i, source, input_text, cfg_cost, coverage_key)

    with open(log_path, "w") as log_f, \
         ThreadPoolExecutor(max_workers=args.workers) as pool:

        # Sliding window: keep at most `workers` futures in flight.
        # Generate inputs lazily so skip_counter adjustments affect next submissions.
        next_i = 0
        active: dict = {}   # future → work_item

        # Prime the pool
        while next_i < min(args.workers, args.iterations):
            item = _next_item(next_i)
            active[pool.submit(_run_one, item)] = item
            next_i += 1

        done_count = 0
        while active:
            finished, _ = wait(active, return_when=FIRST_COMPLETED)

            for fut in finished:
                item = active.pop(fut)
                i, result, java_r, rust_r, source, cfg_cost, input_text, coverage_key = fut.result()
                done_count += 1

                # ---- Update coverage tracker ----
                coverage_tracker.record_key(coverage_key)

                # ---- Update skip counter ----
                triple = _triple_from_text(input_text)
                if result.category == Category.BOTH_TIMEOUT:
                    skip_counter.record_timeout(triple)
                else:
                    skip_counter.record_ok(triple)

                # ---- Submit next item while there's work left ----
                if next_i < args.iterations:
                    new_item = _next_item(next_i)
                    active[pool.submit(_run_one, new_item)] = new_item
                    next_i += 1

                # ---- Accumulate stats ----
                category_counts[result.category.name] += 1
                total_java_time += result.java_time
                total_rust_time += result.rust_time

                # ---- Log ----
                log_entry = {
                    "iter": i,
                    "source": source,
                    "estimated_cost": round(cfg_cost, 3),
                    "category": result.category.name,
                    "max_rel_diff": result.max_rel_diff if not (
                        result.max_rel_diff == float("inf")
                    ) else "inf",
                    "java_time": round(result.java_time, 2),
                    "rust_time": round(result.rust_time, 2),
                    "note": result.note,
                    "triple": list(triple),
                }
                if result.diffs:
                    log_entry["diffs"] = [
                        {"metric": d.name,
                         "java": d.java_val,
                         "rust": d.rust_val,
                         "rel_diff": d.rel_diff if d.rel_diff != float("inf") else "inf"}
                        for d in result.diffs
                    ]
                log_f.write(json.dumps(log_entry) + "\n")
                log_f.flush()

                # ---- Save interesting inputs ----
                save_path = None
                if result.category.interesting or args.save_all:
                    save_path = _save_input(
                        input_text, result, i, run_id, java_r, rust_r
                    )
                    interesting_saved += 1

                # ---- Progress line ----
                marker = "!" if result.category.interesting else " "
                print(f"[{done_count:4d}/{args.iterations}]{marker} {result.summary_line()}"
                      + (f"  -> {save_path.name}" if save_path else ""))

    # ---- Final stats ----
    skip_summary = skip_counter.summary()
    stats = {
        "run_id": run_id,
        "iterations": args.iterations,
        "strategy": args.strategy,
        "budget": args.budget,
        "timeout": args.timeout,
        "interesting_saved": interesting_saved,
        "avg_java_time": total_java_time / args.iterations,
        "avg_rust_time": total_rust_time / args.iterations,
        "categories": dict(category_counts),
        "skip_counter": [
            {"triple": list(t), "multiplier": m, "timeouts": n}
            for t, m, n in skip_summary
        ],
        "coverage": coverage_tracker.summary(),
    }
    stats_path.write_text(json.dumps(stats, indent=2))

    print()
    print("=" * 60)
    print(f"Done.  {interesting_saved} interesting cases saved.")
    print("Category breakdown:")
    for cat, count in sorted(category_counts.items(), key=lambda x: -x[1]):
        bar = "#" * min(count, 40)
        print(f"  {cat:<22} {count:4d}  {bar}")
    if skip_summary:
        print("Skip counter (penalised triples):")
        for triple, mult, n in skip_summary:
            sub, coef, ctype = triple
            print(f"  ({sub or 'none':12} {coef:8} {ctype:10})  "
                  f"multiplier={mult:.4f}  timeouts={n}")
    print(f"Stats: {stats_path}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _coverage_key_from_text(text: str) -> tuple:
    """Extract (subprogram, crystal_type, coefcalc, ddm, has_container, n_segments)."""
    sub = re.search(r'Subprogram\s+(\w+)', text, re.IGNORECASE)
    coef = re.search(r'(?:AbsCoefCalc|CoefCalc)\s+(\w+)', text, re.IGNORECASE)
    ctype = re.search(r'Type\s+(\w+)', text, re.IGNORECASE)
    ddm_m = re.search(r'\bDDM\s+(\w+)', text, re.IGNORECASE)
    has_container = bool(re.search(r'ContainerMaterialType', text, re.IGNORECASE))
    n_segs = max(1, len(re.findall(r'^Beam\s*$', text, re.IGNORECASE | re.MULTILINE)))
    return (
        sub.group(1).upper() if sub else "",
        ctype.group(1).capitalize() if ctype else "Cuboid",
        coef.group(1).upper() if coef else "RD3D",
        ddm_m.group(1).capitalize() if ddm_m else "Simple",
        has_container,
        n_segs,
    )


def _generate(
    strategy: str,
    grammar_gen: GrammarGenerator,
    seeds: list[str],
    rng: random.Random,
    budget: float,
    skip_counter: "SkipCounter | None" = None,
    coverage_tracker: "CoverageTracker | None" = None,
) -> tuple[str, str, float]:
    """Return (source_description, input_text, estimated_cost)."""

    # ---- Chaos mode (dedicated or ~10% of "both") ----
    if strategy == "chaos" or (strategy == "both" and rng.random() < 0.10):
        text = generate_chaos(rng)
        return "chaos", text, _estimate_cost_from_text(text)

    # ---- Structural mode (dedicated) ----
    if strategy == "structural":
        forced = None
        if coverage_tracker and coverage_tracker.should_bias() and rng.random() < 0.3:
            forced = coverage_tracker.least_covered_subprogram()
        cfg = grammar_gen.generate(forced_subprogram=forced)
        cfg = grammar_gen._structural_mutate(cfg)
        return "structural", render(cfg), estimate_cost(cfg)

    use_grammar = (
        strategy == "grammar"
        or (strategy == "both" and rng.random() < 0.5)
        or not seeds
    )

    # ---- Coverage-biased generation ----
    forced_subprogram = None
    if (use_grammar and coverage_tracker and coverage_tracker.should_bias()
            and rng.random() < 0.3):
        forced_subprogram = coverage_tracker.least_covered_subprogram()

    if use_grammar:
        if skip_counter is not None:
            cfg = grammar_gen.generate(forced_subprogram=forced_subprogram)
            triple = (cfg.subprogram, cfg.coefcalc, cfg.crystal_type)
            eff = skip_counter.effective_budget(triple, budget)
            if eff < budget and estimate_cost(cfg) > eff:
                old = grammar_gen.budget
                grammar_gen.budget = eff
                cfg = grammar_gen.generate(forced_subprogram=forced_subprogram)
                grammar_gen.budget = old
        else:
            cfg = grammar_gen.generate(forced_subprogram=forced_subprogram)
        text = render(cfg)
        cost = estimate_cost(cfg)
        return "grammar", text, cost
    else:
        seed_text = rng.choice(seeds)
        mutated = mutate_text(seed_text, mutation_rate=0.2, rng=rng)
        cost = _estimate_cost_from_text(mutated)
        # If cost exceeds budget (or skip-adjusted budget), scale down mutation
        triple = _triple_from_text(mutated)
        eff = skip_counter.effective_budget(triple, budget) if skip_counter else budget
        if cost > eff * 2:
            mutated = mutate_text(seed_text, mutation_rate=0.05, rng=rng)
            cost = _estimate_cost_from_text(mutated)
        return "mutate", mutated, cost


def _estimate_cost_from_text(text: str) -> float:
    """
    Rough cost estimate from raw input text.
    Extracts key cost-driving parameters with regex.
    Handles multiple Wedge blocks by summing their costs.
    """
    import re
    from generate import MC_COST_PER_ELECTRON, XFEL_PER_VOXEL_PER_SECOND, MICROED_MULTIPLIER, INSULIN_BASE_COST

    def _find(pattern, default):
        m = re.search(pattern, text, re.IGNORECASE)
        return float(m.group(1)) if m else default

    # Dimensions: "Dimensions <x> <y> <z>" or "Dimensions <x> <y>" or "Dimensions <x>"
    dim_m = re.search(r'Dimensions\s+([\d.eE+\-]+)\s*([\d.eE+\-]*)\s*([\d.eE+\-]*)',
                       text, re.IGNORECASE)
    if dim_m:
        dx = float(dim_m.group(1)) if dim_m.group(1) else 100.0
        dy = float(dim_m.group(2)) if dim_m.group(2) else dx
        dz = float(dim_m.group(3)) if dim_m.group(3) else dy
    else:
        dx = dy = dz = 100.0

    ppm = _find(r'PixelsPerMicron\s+([\d.eE+\-]+)', 0.5)
    voxels = min(dx * dy * dz * ppm ** 3, 1_000_000)

    # Subprogram
    sub_m = re.search(r'Subprogram\s+(\w+)', text, re.IGNORECASE)
    subprogram = sub_m.group(1).upper() if sub_m else ""

    # Find all Wedge blocks and their AngularResolution values (in order)
    wedge_matches = list(re.finditer(
        r'Wedge\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)', text, re.IGNORECASE
    ))
    res_values = [
        float(m.group(1))
        for m in re.finditer(r'AngularResolution\s+([\d.eE+\-]+)', text, re.IGNORECASE)
    ]

    if subprogram == "XFEL":
        runs = int(_find(r'Runs\s+([\d]+)', 1))
        # Sum all ExposureTime values (one per Wedge block)
        exposure_values = [
            float(m.group(1))
            for m in re.finditer(r'ExposureTime\s+([\d.eE+\-]+)', text, re.IGNORECASE)
        ]
        total_exposure = sum(exposure_values) if exposure_values else 1.0
        return voxels * total_exposure * XFEL_PER_VOXEL_PER_SECOND * runs

    # Sum angle-steps across all Wedge blocks
    if not wedge_matches:
        total_steps = 180.0  # fallback: 360° / 2° resolution
    else:
        total_steps = 0.0
        for i, wm in enumerate(wedge_matches):
            span = abs(float(wm.group(2)) - float(wm.group(1)))
            res = res_values[i] if i < len(res_values) else 2.0
            total_steps += max(1.0, span / res) if span > 0 else 1.0

    base = voxels * total_steps

    if subprogram == "MONTECARLO":
        runs = int(_find(r'Runs\s+([\d]+)', 1))
        sim_e = _find(r'SimElectrons\s+([\d.eE+\-]+)', 1_000_000)
        has_escape = bool(re.search(
            r'Calculate(?:PE|FL)Escape\s+TRUE', text, re.IGNORECASE
        ))
        from generate import _pe_cost_factor
        pe_factor = _pe_cost_factor(has_escape)
        return (base + runs * sim_e * MC_COST_PER_ELECTRON * pe_factor) / INSULIN_BASE_COST
    elif subprogram in ("EMSP", "MICROED"):
        return base * MICROED_MULTIPLIER / INSULIN_BASE_COST
    else:
        return base / INSULIN_BASE_COST


def _save_input(
    input_text: str,
    result,
    iteration: int,
    run_id: str,
    java_r,
    rust_r,
) -> Path:
    """Save an interesting input and its outputs to corpus/."""
    cat = result.category

    if cat in (Category.JAVA_CRASH, Category.RUST_CRASH, Category.BOTH_CRASH):
        subdir = (
            "crashes/java" if cat == Category.JAVA_CRASH else
            "crashes/rust" if cat == Category.RUST_CRASH else
            "crashes/both"
        )
    elif cat in (Category.JAVA_TIMEOUT, Category.RUST_TIMEOUT, Category.BOTH_TIMEOUT):
        subdir = "timeouts"
    elif cat == Category.PERF_DIVERGE:
        subdir = "perf_diverge"
    else:
        subdir = "diffs"

    save_dir = CORPUS_DIR / subdir
    save_dir.mkdir(parents=True, exist_ok=True)

    base = f"{run_id}_{iteration:04d}_{cat.name}"
    input_save = save_dir / f"{base}.txt"
    input_save.write_text(input_text)

    # Save a short metadata sidecar
    meta = {
        "category": cat.name,
        "max_rel_diff": result.max_rel_diff if result.max_rel_diff != float("inf") else "inf",
        "note": result.note,
        "java_time": result.java_time,
        "rust_time": result.rust_time,
        "java_exit": java_r.exit_code,
        "rust_exit": rust_r.exit_code,
        "java_stderr": java_r.stderr[-500:] if java_r.stderr else "",
        "rust_stderr": rust_r.stderr[-500:] if rust_r.stderr else "",
    }
    if result.diffs:
        meta["diffs"] = [
            {"metric": d.name, "java": d.java_val, "rust": d.rust_val,
             "rel_diff": d.rel_diff if d.rel_diff != float("inf") else "inf"}
            for d in result.diffs
        ]
    (save_dir / f"{base}.json").write_text(json.dumps(meta, indent=2))

    return input_save


def _load_seeds() -> list[str]:
    """Load all .txt files from corpus/seeds/."""
    seeds = []
    if SEEDS_DIR.exists():
        for p in sorted(SEEDS_DIR.glob("*.txt")):
            try:
                seeds.append(p.read_text())
            except Exception:
                pass
    return seeds


if __name__ == "__main__":
    main()
