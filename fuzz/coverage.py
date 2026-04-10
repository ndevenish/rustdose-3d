#!/usr/bin/env python3
"""
Offline coverage-driven corpus distillation for RADDOSE-3D Rust.

Workflow
--------
1. Run the existing fuzzer with --save-all to collect a pool of .txt inputs.
2. Run this script on that pool to find which inputs together maximise Rust
   source-level code coverage (line + region).

Algorithm
---------
For each input, run the instrumented Rust binary to collect a .profraw file,
then call llvm-cov export to get a coverage bitmap (set of covered regions).
After all inputs are processed a greedy set-cover pass selects the smallest
subset that covers the most regions.

Usage
-----
    # Build instrumented binary (once)
    cd /path/to/raddose3d
    RUSTFLAGS="-C instrument-coverage" cargo build --release

    # Run coverage analysis on a corpus directory
    cd /path/to/fuzz
    uv run python coverage.py corpus/  [options]

    # Or point at any directory tree of .txt inputs
    uv run python coverage.py corpus/seeds/ corpus/diffs/

Requirements
------------
- Rust binary built with RUSTFLAGS="-C instrument-coverage"
- llvm-tools-preview installed: rustup component add llvm-tools-preview
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).parent.parent
DEFAULT_RUST_BIN = REPO_ROOT / "raddose3d" / "target" / "release" / "raddose3d"
DEFAULT_TIMEOUT = 0.0  # seconds per input; 0 = no timeout


def _find_llvm_tool(name: str) -> Path:
    """Locate an llvm tool from the active rustup toolchain."""
    # Try PATH first (cargo-llvm-cov installs wrappers there)
    import shutil
    if p := shutil.which(name):
        return Path(p)

    # Fall back to rustup toolchain bin directory
    try:
        result = subprocess.run(
            ["rustup", "run", "stable", "--", "rustup", "which", "rustc"],
            capture_output=True, text=True, check=True,
        )
        rustc = Path(result.stdout.strip())
        bin_dir = rustc.parent.parent.parent / "lib" / "rustlib" / \
                  subprocess.run(
                      ["rustc", "-vV"], capture_output=True, text=True
                  ).stdout.split("host: ")[1].split()[0] / "bin"
        candidate = bin_dir / name
        if candidate.exists():
            return candidate
    except Exception:
        pass

    # Direct search under ~/.rustup
    import glob as _glob
    pattern = str(Path.home() / ".rustup" / "toolchains" / "*" / "lib" /
                  "rustlib" / "*" / "bin" / name)
    matches = sorted(_glob.glob(pattern))
    if matches:
        return Path(matches[-1])  # take lexicographically last (newest toolchain)

    raise FileNotFoundError(
        f"Cannot find {name}. Run: rustup component add llvm-tools-preview"
    )


# ---------------------------------------------------------------------------
# Input patching
# ---------------------------------------------------------------------------

import re as _re

def patch_sim_electrons(text: str, n: int) -> str:
    """
    Replace the value on any SimElectrons/SimPhotons line with n.
    If no such line exists and the input contains a MONTECARLO or XFEL
    subprogram block, insert 'SimElectrons <n>' after the Subprogram line.
    """
    # Replace existing SimElectrons/SimPhotons lines
    patched, count = _re.subn(
        r'(?i)^(\s*(?:SimElectrons|SimPhotons)\s+)[^\s#!]+',
        lambda m: m.group(1) + str(n),
        text,
        flags=_re.MULTILINE,
    )
    if count:
        return patched

    # Insert after Subprogram line if MC/XFEL
    def _insert(m):
        return m.group(0) + f"\nSimElectrons {n}"

    patched, count = _re.subn(
        r'(?i)^(\s*Subprogram\s+(?:MONTECARLO|GOS|XFEL)\s*)$',
        _insert,
        text,
        flags=_re.MULTILINE,
    )
    return patched if count else text


# ---------------------------------------------------------------------------
# Per-input coverage collection
# ---------------------------------------------------------------------------

def run_instrumented(
    rust_bin: Path,
    input_path: Path,
    profraw_path: Path,
    timeout: float,
    patched_text: "str | None" = None,
) -> tuple[bool, float, str]:
    """
    Run the instrumented binary on input_path, writing coverage to profraw_path.

    If patched_text is provided it is written to a temp file and used instead
    of input_path directly (used by --patch-sim-electrons).

    Returns (succeeded, wall_time, stderr_snippet).
    """
    env = os.environ.copy()
    env["LLVM_PROFILE_FILE"] = str(profraw_path)

    with tempfile.TemporaryDirectory(prefix="raddose_cov_") as tmp:
        if patched_text is not None:
            actual_input = Path(tmp) / input_path.name
            actual_input.write_text(patched_text)
        else:
            actual_input = input_path.resolve()
        cmd = [
            str(rust_bin),
            "-i", str(actual_input),
            "-p", str(Path(tmp) / "out-"),
        ]
        t0 = time.monotonic()
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                timeout=timeout if timeout > 0 else None,
                env=env,
                cwd=tmp,
            )
            wall = time.monotonic() - t0
            succeeded = proc.returncode == 0
            stderr = proc.stderr.decode(errors="replace")[-400:]
            return succeeded, wall, stderr
        except subprocess.TimeoutExpired:
            wall = time.monotonic() - t0
            return False, wall, "timed out"
        except Exception as e:
            return False, 0.0, str(e)


def merge_profraws(
    profraw_paths: list[Path],
    out_profdata: Path,
    llvm_profdata: Path,
) -> bool:
    """Merge a list of .profraw files into a single .profdata."""
    if not profraw_paths:
        return False
    cmd = [str(llvm_profdata), "merge", "--sparse",
           *[str(p) for p in profraw_paths],
           "-o", str(out_profdata)]
    result = subprocess.run(cmd, capture_output=True)
    return result.returncode == 0


def get_coverage_bitmap(
    rust_bin: Path,
    profdata: Path,
    llvm_cov: Path,
) -> frozenset[tuple]:
    """
    Run llvm-cov export and return a frozenset of covered region entry points.

    Each ID is (filename, line, col) for each segment that is a region-entry
    with execution count > 0. Only source files from the workspace are included
    (registry crates and compiler internals are skipped).
    """
    cmd = [
        str(llvm_cov), "export",
        "--instr-profile", str(profdata),
        "--format", "text",
        str(rust_bin),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return frozenset()

    covered = set()
    try:
        data = json.loads(result.stdout)
        for file_entry in data.get("data", [{}])[0].get("files", []):
            fname = file_entry.get("filename", "")
            # Skip compiler internals, registry crates, and test files
            if ("/rustc/" in fname
                    or "/.cargo/registry/" in fname
                    or fname.endswith("_test.rs")):
                continue
            for segment in file_entry.get("segments", []):
                # Segment format: [line, col, count, has_count, is_region_entry, is_gap]
                # is_region_entry (index 4) marks the start of a new coverage region.
                if len(segment) >= 5 and segment[4] and segment[2] > 0:
                    covered.add((fname, segment[0], segment[1]))
    except (json.JSONDecodeError, KeyError, IndexError):
        pass

    return frozenset(covered)


# ---------------------------------------------------------------------------
# Greedy set-cover
# ---------------------------------------------------------------------------

def greedy_set_cover(
    inputs: list[Path],
    bitmaps: dict[Path, frozenset],
    wall_times: "dict[Path, float] | None" = None,
) -> list[Path]:
    """
    Greedy maximum-coverage selection.

    Iteratively picks the input that covers the most uncovered regions,
    until no new regions can be added. Among inputs with equal marginal gain,
    the faster one (by wall_time) is preferred. Returns paths in selection order.
    """
    universe: set[tuple] = set()
    for bm in bitmaps.values():
        universe |= bm

    if not universe:
        return []

    times = wall_times or {}
    covered: set[tuple] = set()
    selected: list[Path] = []
    remaining = list(inputs)

    while remaining and covered != universe:
        # Primary key: most new regions. Tiebreaker: shortest wall time.
        best = max(remaining, key=lambda p: (len(bitmaps[p] - covered), -times.get(p, 0.0)))
        gain = len(bitmaps[best] - covered)
        if gain == 0:
            break
        covered |= bitmaps[best]
        selected.append(best)
        remaining.remove(best)

    return selected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def collect_inputs(sources: list[Path]) -> list[Path]:
    """Recursively find all .txt files under the given paths."""
    inputs = []
    for src in sources:
        if src.is_file() and src.suffix == ".txt":
            inputs.append(src)
        elif src.is_dir():
            inputs.extend(sorted(src.rglob("*.txt")))
    return inputs


def main():
    ap = argparse.ArgumentParser(
        description="Coverage-driven corpus distillation for RADDOSE-3D Rust",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "sources", nargs="+", type=Path,
        help="Directory or .txt file(s) to analyse. Directories are searched recursively.",
    )
    ap.add_argument(
        "--rust-bin", type=Path, default=DEFAULT_RUST_BIN,
        help=f"Path to instrumented Rust binary (default: {DEFAULT_RUST_BIN})",
    )
    ap.add_argument(
        "--timeout", type=float, default=DEFAULT_TIMEOUT,
        help="Per-input timeout in seconds. 0 = no timeout (default). "
             "Coverage requires a clean process exit to flush the profraw, so "
             "killed processes produce no coverage data. Only set a timeout if "
             "you have genuinely unbounded inputs in your corpus.",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=None,
        help="Write selected inputs here (copies, not moves). Default: print list only.",
    )
    ap.add_argument(
        "--report", type=Path, default=None,
        help="Write a JSON coverage report to this path.",
    )
    ap.add_argument(
        "--keep-profraws", action="store_true",
        help="Keep per-input .profraw files in a temp dir (printed at end).",
    )
    ap.add_argument(
        "--workers", "-j", type=int, default=1,
        help="Parallel workers for running the instrumented binary (default: 1). "
             "Each worker needs its own LLVM_PROFILE_FILE so files don't collide.",
    )
    ap.add_argument(
        "--patch-sim-electrons", type=int, default=1000, metavar="N",
        help="Override SimElectrons/SimPhotons to N before running (default: 1000). "
             "Reduces MONTECARLO/XFEL runtime from hours to seconds without changing "
             "which code paths are exercised. Non-MC/XFEL inputs are unaffected. "
             "Set to 0 to disable.",
    )
    args = ap.parse_args()

    # Validate binary exists
    if not args.rust_bin.exists():
        print(f"ERROR: Rust binary not found: {args.rust_bin}", file=sys.stderr)
        print("Build it with:", file=sys.stderr)
        print("  cd raddose3d", file=sys.stderr)
        print('  RUSTFLAGS="-C instrument-coverage" cargo build --release', file=sys.stderr)
        sys.exit(1)

    # Find llvm tools
    try:
        llvm_profdata = _find_llvm_tool("llvm-profdata")
        llvm_cov = _find_llvm_tool("llvm-cov")
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"llvm-profdata: {llvm_profdata}")
    print(f"llvm-cov:      {llvm_cov}")
    print()

    # Collect inputs
    inputs = collect_inputs(args.sources)
    if not inputs:
        print("No .txt inputs found.", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(inputs)} input(s).")

    # Work dir for profraw files
    work_dir_obj = tempfile.TemporaryDirectory(prefix="raddose_cov_work_")
    work_dir = Path(work_dir_obj.name)

    # ---------------------------------------------------------------------------
    # Phase 1: run each input and collect individual coverage bitmaps
    # ---------------------------------------------------------------------------
    print(f"\nPhase 1: collecting coverage ({args.workers} worker(s), "
          f"timeout={args.timeout}s each)")
    print("-" * 60)

    bitmaps: dict[Path, frozenset] = {}
    wall_times: dict[Path, float] = {}
    skipped: list[Path] = []
    per_input_profraw: dict[Path, Path] = {}

    # Sequential implementation (workers=1) is simplest and avoids profraw
    # collision issues. Parallel path uses one profraw file per worker slot.
    if args.workers <= 1:
        for idx, inp in enumerate(inputs, 1):
            profraw = work_dir / f"{idx:05d}.profraw"
            profdata = work_dir / f"{idx:05d}.profdata"

            patched = (
                patch_sim_electrons(inp.read_text(), args.patch_sim_electrons)
                if args.patch_sim_electrons else None
            )
            succeeded, wall, stderr = run_instrumented(
                args.rust_bin, inp, profraw, args.timeout, patched_text=patched
            )

            wall_times[inp] = wall
            status = "OK" if succeeded else "FAIL"
            print(f"[{idx:4d}/{len(inputs)}] {status}  {wall:.1f}s  {inp.name}")
            if not succeeded and stderr:
                print(f"           stderr: {stderr[:120]}")

            if profraw.exists():
                if merge_profraws([profraw], profdata, llvm_profdata):
                    bm = get_coverage_bitmap(args.rust_bin, profdata, llvm_cov)
                    bitmaps[inp] = bm
                    per_input_profraw[inp] = profraw
                else:
                    skipped.append(inp)
                    bitmaps[inp] = frozenset()
            else:
                skipped.append(inp)
                bitmaps[inp] = frozenset()
    else:
        # Parallel: use ThreadPoolExecutor; each thread gets a unique slot index
        # for its profraw filename so files can't collide.
        import threading
        from concurrent.futures import ThreadPoolExecutor, as_completed

        slot_lock = threading.Lock()
        slot_counter = [0]

        def _run_one(inp: Path):
            with slot_lock:
                slot = slot_counter[0]
                slot_counter[0] += 1
            profraw = work_dir / f"{slot:05d}.profraw"
            profdata = work_dir / f"{slot:05d}.profdata"
            patched = (
                patch_sim_electrons(inp.read_text(), args.patch_sim_electrons)
                if args.patch_sim_electrons else None
            )
            succeeded, wall, stderr = run_instrumented(
                args.rust_bin, inp, profraw, args.timeout, patched_text=patched
            )
            bm = frozenset()
            if profraw.exists():
                if merge_profraws([profraw], profdata, llvm_profdata):
                    bm = get_coverage_bitmap(args.rust_bin, profdata, llvm_cov)
            return inp, succeeded, wall, stderr, bm, profraw

        done = 0
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(_run_one, inp): inp for inp in inputs}
            for fut in as_completed(futures):
                inp, succeeded, wall, stderr, bm, profraw = fut.result()
                done += 1
                status = "OK" if succeeded else "FAIL"
                print(f"[{done:4d}/{len(inputs)}] {status}  {wall:.1f}s  {inp.name}")
                if not succeeded and stderr:
                    print(f"           stderr: {stderr[:120]}")
                bitmaps[inp] = bm
                wall_times[inp] = wall
                per_input_profraw[inp] = profraw
                if not bm:
                    skipped.append(inp)

    # ---------------------------------------------------------------------------
    # Phase 2: greedy set-cover
    # ---------------------------------------------------------------------------
    print()
    print("Phase 2: greedy set-cover")
    print("-" * 60)

    all_regions = set()
    for bm in bitmaps.values():
        all_regions |= bm
    total_regions = len(all_regions)

    usable = [p for p in inputs if bitmaps[p]]
    selected = greedy_set_cover(usable, bitmaps, wall_times)

    selected_regions: set[tuple] = set()
    for p in selected:
        selected_regions |= bitmaps[p]

    print(f"Total regions reachable across all inputs: {total_regions:,}")
    print(f"Inputs with usable coverage:               {len(usable):,}")
    print(f"Inputs skipped (no profraw / merge fail):  {len(skipped):,}")
    print(f"Selected by greedy set-cover:              {len(selected):,} input(s)")
    print(f"Regions covered by selection:              {len(selected_regions):,} "
          f"({100*len(selected_regions)/max(total_regions,1):.1f}%)")

    print()
    print("Selected inputs (in selection order):")
    cumulative = set()
    for i, p in enumerate(selected, 1):
        gain = len(bitmaps[p] - cumulative)
        cumulative |= bitmaps[p]
        t = wall_times.get(p, 0.0)
        print(f"  {i:3d}.  +{gain:6,} regions  {t:5.1f}s  {p}")

    # ---------------------------------------------------------------------------
    # Phase 3: merge selected profraws → total program coverage (always)
    # ---------------------------------------------------------------------------
    merged_all_profdata = work_dir / "all_selected.profdata"
    selected_profraws = [
        per_input_profraw[p] for p in selected
        if p in per_input_profraw and per_input_profraw[p].exists()
    ]
    cov_data: dict = {}
    if selected_profraws and merge_profraws(selected_profraws, merged_all_profdata, llvm_profdata):
        cmd = [
            str(llvm_cov), "export",
            "--instr-profile", str(merged_all_profdata),
            "--format", "text",
            str(args.rust_bin),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            try:
                cov_data = json.loads(result.stdout)
            except json.JSONDecodeError:
                pass

    totals = cov_data.get("data", [{}])[0].get("totals", {})
    if totals:
        def _pct(d: dict) -> str:
            return f"{d.get('covered', 0):,}/{d.get('count', 0):,} ({d.get('percent', 0.0):.1f}%)"
        print()
        print("Total program coverage (selected inputs merged):")
        print(f"  Lines:    {_pct(totals.get('lines', {}))}")
        print(f"  Regions:  {_pct(totals.get('regions', {}))}")
        print(f"  Branches: {_pct(totals.get('branches', {}))}")

    # ---------------------------------------------------------------------------
    # Phase 4: copy selected inputs to --out-dir (optional)
    # ---------------------------------------------------------------------------
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        import shutil
        for i, p in enumerate(selected, 1):
            dest = args.out_dir / f"{i:03d}_{p.name}"
            shutil.copy2(p, dest)
        print(f"\nCopied {len(selected)} input(s) to {args.out_dir}")

    # ---------------------------------------------------------------------------
    # Phase 5: JSON report (optional)
    # ---------------------------------------------------------------------------
    if args.report:
        file_coverage: dict[str, dict] = {}
        for file_entry in cov_data.get("data", [{}])[0].get("files", []):
            fname = file_entry.get("filename", "")
            if "/rustc/" in fname or "/.cargo/registry/" in fname or fname.endswith("_test.rs"):
                continue
            summary = file_entry.get("summary", {})
            file_coverage[fname] = {
                "lines": summary.get("lines", {}),
                "regions": summary.get("regions", {}),
                "branches": summary.get("branches", {}),
            }

        report = {
            "inputs_analysed": len(inputs),
            "inputs_skipped": len(skipped),
            "inputs_selected": len(selected),
            "total_regions": total_regions,
            "covered_regions": len(selected_regions),
            "coverage_pct": round(100 * len(selected_regions) / max(total_regions, 1), 2),
            "totals": totals,
            "selected": [
                {
                    "path": str(p),
                    "wall_time": wall_times.get(p, 0.0),
                    "regions_covered": len(bitmaps[p]),
                    "marginal_gain": len(bitmaps[p] - (
                        set().union(*[bitmaps[q] for q in selected[:i]])
                    )),
                }
                for i, p in enumerate(selected)
            ],
            "per_file_coverage": file_coverage,
        }
        args.report.write_text(json.dumps(report, indent=2))
        print(f"Report written to {args.report}")

    # Cleanup (unless --keep-profraws)
    if args.keep_profraws:
        work_dir_obj._finalizer.detach()  # prevent cleanup
        print(f"\nProfraw files retained in: {work_dir}")
    else:
        work_dir_obj.cleanup()


if __name__ == "__main__":
    main()
