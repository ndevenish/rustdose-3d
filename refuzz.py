#!/usr/bin/env python3
"""
Re-run corpus inputs through both implementations and re-categorize.

For each .txt input found, reruns Java and Rust, then:
  - Removes inputs that now produce MATCH (fixed)
  - Moves inputs to the correct subdirectory if the category changed
  - Updates the .json sidecar with fresh results
  - Reports a diff of old → new categories

Usage:
    uv run python refuzz.py corpus/crashes/java/
    uv run python refuzz.py corpus/
    uv run python refuzz.py 'corpus/crashes/**/*.txt'
    uv run python refuzz.py corpus/crashes/java/ --dry-run
    uv run python refuzz.py corpus/ --workers 8 --timeout 120
"""

import argparse
import glob as glob_module
import json
import sys
import tempfile
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from compare import Category, Comparison, compare
from harness import DEFAULT_JAVA_JAR, DEFAULT_RUST_BIN, DEFAULT_TIMEOUT, run_both

FUZZ_DIR = Path(__file__).parent
CORPUS_DIR = FUZZ_DIR / "corpus"

# Re-run these categories — skip PARSE_ERROR / HARNESS_ERROR inputs from previous runs
# (they might have been saved manually or be malformed)
_SKIP_CATEGORIES = {Category.PARSE_ERROR, Category.HARNESS_ERROR}


def _subdir_for(cat: Category) -> str | None:
    """Return corpus subdirectory for category, or None to delete (fixed)."""
    if cat == Category.MATCH:
        return None  # fixed — remove
    if cat == Category.JAVA_CRASH:
        return "crashes/java"
    if cat == Category.RUST_CRASH:
        return "crashes/rust"
    if cat == Category.BOTH_CRASH:
        return "crashes/both"
    if cat in (Category.JAVA_TIMEOUT, Category.RUST_TIMEOUT, Category.BOTH_TIMEOUT):
        return "timeouts"
    if cat == Category.PERF_DIVERGE:
        return "perf_diverge"
    return "diffs"  # MAJOR_DIFF, NAN_INF, MINOR_DIFF


def _old_category_name(txt_path: Path) -> str | None:
    """Read old category from .json sidecar, else infer from filename stem."""
    json_path = txt_path.with_suffix(".json")
    if json_path.exists():
        try:
            return json.loads(json_path.read_text()).get("category")
        except Exception:
            pass
    stem = txt_path.stem
    for cat in Category:
        if stem.endswith(f"_{cat.name}"):
            return cat.name
    return None


def _rename_stem(stem: str, new_cat_name: str) -> str:
    """Replace the category token in a filename stem with new_cat_name.

    Handles collision suffixes like _1 or _1_1 that may follow the category.
    """
    for cat in Category:
        marker = f"_{cat.name}"
        idx = stem.find(marker)
        if idx != -1:
            suffix = stem[idx + len(marker):]  # e.g. "", "_1", "_1_1"
            prefix = stem[:idx]
            return f"{prefix}_{new_cat_name}{suffix}"
    # Filename doesn't follow convention — append new category
    return f"{stem}_{new_cat_name}"


def _collect_inputs(paths: list[str]) -> list[Path]:
    seen: set[Path] = set()
    result: list[Path] = []
    for p in paths:
        path = Path(p)
        if path.is_dir():
            candidates = sorted(path.rglob("*.txt"))
        else:
            candidates = sorted(Path(m) for m in glob_module.glob(p, recursive=True)
                                if Path(m).suffix == ".txt")
        for c in candidates:
            if c not in seen:
                seen.add(c)
                result.append(c)
    return result


def _outcome_label(old: str | None, new: Category) -> str:
    if new == Category.MATCH:
        return "FIXED"
    if new in (Category.PARSE_ERROR, Category.HARNESS_ERROR):
        return "ERROR"
    if old is None:
        return "new"
    try:
        old_cat = Category[old]
    except KeyError:
        return "?"
    if old_cat == new:
        return "-"
    # Rank categories by severity (higher = worse)
    severity = {
        Category.MATCH: 0,
        Category.MINOR_DIFF: 1,
        Category.PERF_DIVERGE: 2,
        Category.MAJOR_DIFF: 3,
        Category.NAN_INF: 3,
        Category.BOTH_CRASH: 4,
        Category.JAVA_CRASH: 4,
        Category.RUST_CRASH: 4,
        Category.BOTH_TIMEOUT: 4,
        Category.JAVA_TIMEOUT: 4,
        Category.RUST_TIMEOUT: 4,
        Category.PARSE_ERROR: 5,
        Category.HARNESS_ERROR: 5,
    }
    old_sev = severity.get(old_cat, 3)
    new_sev = severity.get(new, 3)
    if new_sev < old_sev:
        return "improved"
    if new_sev > old_sev:
        return "REGRESSED"
    return "changed"  # same severity, different category


def _update_corpus(
    txt_path: Path,
    new_result: Comparison,
    java_r,
    rust_r,
    dry_run: bool,
) -> tuple[str, Path | None]:
    """
    Move/delete txt_path (and its .json sidecar) to match new category.
    Returns (action, new_path_or_None).
    """
    new_subdir = _subdir_for(new_result.category)
    json_path = txt_path.with_suffix(".json")

    # Build updated sidecar metadata
    new_meta = {
        "category": new_result.category.name,
        "max_rel_diff": (new_result.max_rel_diff
                         if new_result.max_rel_diff != float("inf") else "inf"),
        "note": new_result.note,
        "java_time": new_result.java_time,
        "rust_time": new_result.rust_time,
        "java_exit": java_r.exit_code,
        "rust_exit": rust_r.exit_code,
        "java_stderr": java_r.stderr[-500:] if java_r.stderr else "",
        "rust_stderr": rust_r.stderr[-500:] if rust_r.stderr else "",
    }
    if new_result.diffs:
        new_meta["diffs"] = [
            {"metric": d.name, "java": d.java_val, "rust": d.rust_val,
             "rel_diff": d.rel_diff if d.rel_diff != float("inf") else "inf"}
            for d in new_result.diffs
        ]

    if new_subdir is None:
        # Fixed — delete
        if not dry_run:
            txt_path.unlink(missing_ok=True)
            json_path.unlink(missing_ok=True)
        return "deleted", None

    dest_dir = CORPUS_DIR / new_subdir
    new_stem = _rename_stem(txt_path.stem, new_result.category.name)
    dest_txt = dest_dir / f"{new_stem}.txt"

    if dest_txt.resolve() == txt_path.resolve():
        # Same location — just update the sidecar
        if not dry_run:
            json_path.write_text(json.dumps(new_meta, indent=2))
        return "updated", txt_path

    # Different location — move files
    if not dry_run:
        dest_dir.mkdir(parents=True, exist_ok=True)
        # Avoid clobbering an existing file with the same name
        if dest_txt.exists():
            for i in range(1, 1000):
                candidate = dest_dir / f"{new_stem}_{i}.txt"
                if not candidate.exists():
                    dest_txt = candidate
                    break
        txt_path.rename(dest_txt)
        if json_path.exists():
            json_path.rename(dest_txt.with_suffix(".json"))
        dest_txt.with_suffix(".json").write_text(json.dumps(new_meta, indent=2))
    return "moved", dest_txt


def main():
    ap = argparse.ArgumentParser(description="Re-run corpus inputs and re-categorize")
    ap.add_argument("paths", nargs="+",
                    help="Directories or glob patterns of .txt inputs")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would happen without making changes")
    ap.add_argument("--workers", "-j", type=int, default=4)
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT)
    ap.add_argument("--java-jar", type=Path, default=DEFAULT_JAVA_JAR)
    ap.add_argument("--rust-bin", type=Path, default=DEFAULT_RUST_BIN)
    args = ap.parse_args()

    inputs = _collect_inputs(args.paths)
    if not inputs:
        print("No .txt inputs found.")
        sys.exit(0)

    print(f"Re-fuzzing {len(inputs)} input(s)  |  workers={args.workers}"
          f"  |  timeout={args.timeout}s"
          + ("  |  DRY RUN" if args.dry_run else ""))
    print()

    outcome_counts: Counter = Counter()
    old_to_new: Counter = Counter()  # (old_name, new_name) counts

    def _run_one(txt_path: Path):
        input_text = txt_path.read_text()
        with tempfile.TemporaryDirectory(prefix="raddose_refuzz_") as tmp:
            tmp_path = Path(tmp)
            input_copy = tmp_path / "input.txt"
            input_copy.write_text(input_text)
            java_r, rust_r = run_both(
                input_copy, tmp_path,
                java_jar=args.java_jar,
                rust_bin=args.rust_bin,
                timeout=args.timeout,
            )
            result = compare(java_r, rust_r)
        return txt_path, result, java_r, rust_r

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        active = {pool.submit(_run_one, p): p for p in inputs[:args.workers]}
        next_idx = args.workers
        done_count = 0

        while active:
            finished, _ = wait(active, return_when=FIRST_COMPLETED)
            for fut in finished:
                active.pop(fut)
                done_count += 1

                if next_idx < len(inputs):
                    p = inputs[next_idx]
                    active[pool.submit(_run_one, p)] = p
                    next_idx += 1

                txt_path, new_result, java_r, rust_r = fut.result()
                old_cat_name = _old_category_name(txt_path)
                new_cat = new_result.category
                label = _outcome_label(old_cat_name, new_cat)

                # Update corpus
                if new_cat not in _SKIP_CATEGORIES:
                    action, new_path = _update_corpus(
                        txt_path, new_result, java_r, rust_r, args.dry_run
                    )
                else:
                    action, new_path = "skipped", txt_path

                outcome_counts[label] += 1
                old_to_new[(old_cat_name or "?", new_cat.name)] += 1

                # Progress line
                old_str = f"{old_cat_name or '?':15}"
                new_str = f"{new_cat.name:15}"
                marker = "!" if label in ("FIXED", "REGRESSED") else " "
                extra = ""
                if action == "moved" and new_path:
                    extra = f"  → {new_path.relative_to(CORPUS_DIR)}"
                elif action == "deleted":
                    extra = "  (removed)"
                print(f"[{done_count:4d}/{len(inputs)}]{marker}"
                      f" {old_str} → {new_str} [{label}]"
                      f"  {txt_path.name[:50]}{extra}")

    # Summary
    print()
    print("=" * 60)
    for label in ("FIXED", "improved", "changed", "REGRESSED", "-", "ERROR", "skipped", "new"):
        count = outcome_counts[label]
        if count:
            print(f"  {label:<12} {count:4d}")

    if old_to_new:
        print()
        print("Transitions (old → new):")
        for (old, new), count in sorted(old_to_new.items(), key=lambda x: -x[1]):
            if old != new:
                print(f"  {old:20} → {new:20}  ×{count}")


if __name__ == "__main__":
    main()
