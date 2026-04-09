"""
Compare Java and Rust RADDOSE-3D output files and classify the result.
"""

import csv
import math
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Optional

from harness import RunResult


class Category(Enum):
    MATCH = auto()           # All metrics agree within tolerance
    MINOR_DIFF = auto()      # Max relative diff in (tol_minor, tol_major)
    MAJOR_DIFF = auto()      # Max relative diff >= tol_major
    NAN_INF = auto()         # NaN or Inf in one but not the other
    JAVA_CRASH = auto()      # Java non-zero exit, Rust succeeded
    RUST_CRASH = auto()      # Rust non-zero exit, Java succeeded
    BOTH_CRASH = auto()      # Both non-zero exit
    JAVA_TIMEOUT = auto()
    RUST_TIMEOUT = auto()
    BOTH_TIMEOUT = auto()
    PERF_DIVERGE = auto()    # One timed out, other succeeded
    PARSE_ERROR = auto()     # Outputs missing or unparseable
    HARNESS_ERROR = auto()   # Binary not found etc.

    @property
    def interesting(self) -> bool:
        """Categories worth saving to corpus."""
        return self in {
            Category.MAJOR_DIFF, Category.NAN_INF,
            Category.JAVA_CRASH, Category.RUST_CRASH,
            Category.JAVA_TIMEOUT, Category.RUST_TIMEOUT,
            Category.BOTH_TIMEOUT, Category.PERF_DIVERGE,
        }


# Metrics from Summary.csv that we compare.
# Keys match the CSV header (stripped).
COMPARE_METRICS = [
    "Average DWD",
    "Last DWD",
    "Elastic Yield (wedge)",
    "Diffraction Efficiency",
    "AD-WC",
    "Max Dose",
    "Wedge Absorbed Energy",
    "Dose Inefficiency",
    "Dose Inefficiency PE",
]

TOL_MATCH = 1e-4      # 0.01%  — differences below this are "match"
TOL_MINOR = 1e-2      # 1%     — differences below this are "minor"
# differences >= TOL_MINOR are "major"
PERF_DIVERGE_RATIO = 3.0   # wall-time ratio to flag as PERF_DIVERGE


@dataclass
class MetricDiff:
    name: str
    java_val: float
    rust_val: float
    rel_diff: float        # |java - rust| / max(|java|, |rust|, 1e-30)
    is_nan_inf: bool


@dataclass
class Comparison:
    category: Category
    max_rel_diff: float = 0.0
    diffs: list[MetricDiff] = field(default_factory=list)
    java_time: float = 0.0
    rust_time: float = 0.0
    note: str = ""

    def summary_line(self) -> str:
        parts = [self.category.name]
        if self.diffs:
            worst = max(self.diffs, key=lambda d: d.rel_diff)
            parts.append(f"max_diff={self.max_rel_diff:.2e} [{worst.name}]")
        if self.note:
            parts.append(self.note)
        parts.append(f"java={self.java_time:.1f}s rust={self.rust_time:.1f}s")
        return "  ".join(parts)


def compare(java_result: RunResult, rust_result: RunResult) -> Comparison:
    """
    Compare two RunResults and return a classified Comparison.
    """
    jt = java_result.wall_time
    rt = rust_result.wall_time

    # ---- Timeout / crash classification first ----
    if java_result.error or rust_result.error:
        note = " | ".join(filter(None, [java_result.error, rust_result.error]))
        return Comparison(Category.HARNESS_ERROR, note=note, java_time=jt, rust_time=rt)

    if java_result.timed_out and rust_result.timed_out:
        return Comparison(Category.BOTH_TIMEOUT, java_time=jt, rust_time=rt)

    if java_result.timed_out:
        return Comparison(Category.JAVA_TIMEOUT,
                          note="java timed out, rust succeeded" if rust_result.succeeded else "",
                          java_time=jt, rust_time=rt)

    if rust_result.timed_out:
        return Comparison(Category.RUST_TIMEOUT,
                          note="rust timed out, java succeeded" if java_result.succeeded else "",
                          java_time=jt, rust_time=rt)

    if java_result.crashed and rust_result.crashed:
        return Comparison(Category.BOTH_CRASH,
                          note=f"java exit={java_result.exit_code} rust exit={rust_result.exit_code}",
                          java_time=jt, rust_time=rt)

    if java_result.crashed:
        return Comparison(Category.JAVA_CRASH,
                          note=f"exit={java_result.exit_code}",
                          java_time=jt, rust_time=rt)

    if rust_result.crashed:
        return Comparison(Category.RUST_CRASH,
                          note=f"exit={rust_result.exit_code}",
                          java_time=jt, rust_time=rt)

    # ---- Both succeeded: compare Summary.csv ----
    java_csv = java_result.summary_csv_path()
    rust_csv = rust_result.summary_csv_path()

    if java_csv is None or rust_csv is None:
        missing = []
        if java_csv is None:
            missing.append("java Summary.csv")
        if rust_csv is None:
            missing.append("rust Summary.csv")
        return Comparison(Category.PARSE_ERROR,
                          note=f"missing: {', '.join(missing)}",
                          java_time=jt, rust_time=rt)

    try:
        java_rows = _parse_summary_csv(java_csv)
        rust_rows = _parse_summary_csv(rust_csv)
    except Exception as e:
        return Comparison(Category.PARSE_ERROR, note=str(e),
                          java_time=jt, rust_time=rt)

    if len(java_rows) != len(rust_rows):
        return Comparison(
            Category.PARSE_ERROR,
            note=f"wedge row count mismatch: java={len(java_rows)} rust={len(rust_rows)}",
            java_time=jt, rust_time=rt,
        )

    multi_wedge = len(java_rows) > 1
    diffs = []
    for row_idx, (java_metrics, rust_metrics) in enumerate(zip(java_rows, rust_rows)):
        for name in COMPARE_METRICS:
            j = java_metrics.get(name)
            r = rust_metrics.get(name)
            if j is None or r is None:
                continue
            label = f"{name}[{row_idx}]" if multi_wedge else name
            nan_inf = _is_nan_inf(j) or _is_nan_inf(r)
            if nan_inf and not (_is_nan_inf(j) and _is_nan_inf(r)):
                # One has NaN/Inf, other doesn't
                diffs.append(MetricDiff(label, j, r, rel_diff=float("inf"), is_nan_inf=True))
            elif not (math.isnan(j) or math.isnan(r)):
                rel = _rel_diff(j, r)
                diffs.append(MetricDiff(label, j, r, rel_diff=rel, is_nan_inf=False))

    if not diffs:
        return Comparison(Category.PARSE_ERROR, note="no comparable metrics found",
                          java_time=jt, rust_time=rt)

    # Check for NaN/Inf disagreements first
    nan_inf_diffs = [d for d in diffs if d.is_nan_inf]
    if nan_inf_diffs:
        return Comparison(Category.NAN_INF, diffs=diffs,
                          max_rel_diff=float("inf"),
                          java_time=jt, rust_time=rt)

    max_rel = max(d.rel_diff for d in diffs)
    if max_rel < TOL_MATCH:
        cat = Category.MATCH
    elif max_rel < TOL_MINOR:
        cat = Category.MINOR_DIFF
    else:
        cat = Category.MAJOR_DIFF

    # Check for large runtime divergence even when outputs match.
    # Only flag if at least one run exceeds 10s (avoids noise on fast inputs).
    if jt > 0.5 and rt > 0.5 and max(jt, rt) > 10.0:
        ratio = max(jt, rt) / min(jt, rt)
        if ratio >= PERF_DIVERGE_RATIO:
            slower = "java" if jt > rt else "rust"
            return Comparison(Category.PERF_DIVERGE, diffs=diffs, max_rel_diff=max_rel,
                              note=f"{slower} {ratio:.1f}x slower",
                              java_time=jt, rust_time=rt)

    return Comparison(cat, diffs=diffs, max_rel_diff=max_rel,
                      java_time=jt, rust_time=rt)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_summary_csv(path: Path) -> list[dict[str, float]]:
    """Parse a RADDOSE-3D Summary.csv into a list of {metric_name: value} dicts, one per wedge row."""
    with open(path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        raise ValueError(f"Empty CSV: {path}")

    result = []
    for row in rows:
        d = {}
        for key, val in row.items():
            if key is None:
                continue
            clean_key = key.strip()
            try:
                d[clean_key] = float(val.strip())
            except (ValueError, AttributeError):
                pass
        result.append(d)
    return result


def _rel_diff(a: float, b: float) -> float:
    """Relative difference, safe for near-zero values.

    Absolute differences below ABS_FLOOR are treated as zero — this avoids
    false positives when Java's 6 d.p. CSV output truncates a tiny value to
    0.000000 while Rust emits full precision (e.g. 3.88e-7 vs 0.0).
    """
    if abs(a - b) < 5e-7:
        return 0.0
    denom = max(abs(a), abs(b), 1e-30)
    return abs(a - b) / denom


def _is_nan_inf(v: float) -> bool:
    return math.isnan(v) or math.isinf(v)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <java-Summary.csv> <rust-Summary.csv>", file=sys.stderr)
        sys.exit(1)

    java_path = Path(sys.argv[1])
    rust_path = Path(sys.argv[2])

    # Build minimal RunResult objects that point at the given CSV files.
    # summary_csv_path() looks for "{impl}-Summary.csv" in output_dir, so
    # place a symlink-free workaround by subclassing is not needed — instead
    # we temporarily patch the paths via a thin wrapper.
    from harness import RunResult

    class _FileRunResult(RunResult):
        _csv: Optional[Path] = None

        def summary_csv_path(self) -> Optional[Path]:
            return self._csv

    def _make(impl: str, csv_path: Path) -> "_FileRunResult":
        r = _FileRunResult.__new__(_FileRunResult)
        r.impl = impl
        r.exit_code = 0
        r.wall_time = 0.0
        r.stdout = ""
        r.stderr = ""
        r.output_dir = csv_path.parent
        r.timed_out = False
        r.error = ""
        r._csv = csv_path
        return r

    java_r = _make("java", java_path)
    rust_r = _make("rust", rust_path)

    result = compare(java_r, rust_r)
    print(result.summary_line())
