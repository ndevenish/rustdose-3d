"""
Hypothesis-based differential test for RADDOSE-3D Java vs Rust.

Run with:
    uv run pytest test_differential.py -v
    uv run pytest test_differential.py -v --hypothesis-seed=12345

Hypothesis automatically shrinks failing examples to the minimal
reproducer and persists them in .hypothesis/ for replay.

Settings can be tuned via environment variables:
    RADDOSE_FUZZ_EXAMPLES=500    number of examples per test (default 100)
    RADDOSE_FUZZ_TIMEOUT=60      per-run timeout in seconds (default 60)
    RADDOSE_FUZZ_BUDGET=5.0      cost budget relative to insulin (default 5.0)
"""

import os
import sys
import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings, assume
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).parent))

from compare import Category, compare
from generate import DEFAULT_BUDGET, render
from harness import DEFAULT_JAVA_JAR, DEFAULT_RUST_BIN, DEFAULT_TIMEOUT, run_both
from strategies import raddose_config

# ---------------------------------------------------------------------------
# Configuration via environment
# ---------------------------------------------------------------------------

MAX_EXAMPLES = int(os.environ.get("RADDOSE_FUZZ_EXAMPLES", 100))
TIMEOUT = float(os.environ.get("RADDOSE_FUZZ_TIMEOUT", DEFAULT_TIMEOUT))
BUDGET = float(os.environ.get("RADDOSE_FUZZ_BUDGET", DEFAULT_BUDGET))

# Categories that constitute a test failure — these are the interesting ones
FAIL_CATEGORIES = {
    Category.MAJOR_DIFF,
    Category.NAN_INF,
    Category.JAVA_CRASH,
    Category.RUST_CRASH,
    Category.BOTH_CRASH,
    Category.JAVA_TIMEOUT,   # one timed out, other succeeded
    Category.RUST_TIMEOUT,
}

# Categories to skip (both timed out = bad budget, not a real failure)
SKIP_CATEGORIES = {
    Category.BOTH_TIMEOUT,
    Category.HARNESS_ERROR,
}

# ---------------------------------------------------------------------------
# Shared settings profile
# ---------------------------------------------------------------------------

fuzz_settings = settings(
    max_examples=MAX_EXAMPLES,
    deadline=None,          # each example runs subprocesses; no per-example time limit
    suppress_health_check=[
        HealthCheck.too_slow,
        HealthCheck.filter_too_much,
    ],
)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _run_and_compare(cfg):
    """Render config, run both implementations, return Comparison."""
    text = render(cfg)
    with tempfile.TemporaryDirectory(prefix="raddose_hyp_") as tmp:
        input_path = Path(tmp) / "input.txt"
        input_path.write_text(text)
        java_r, rust_r = run_both(
            input_path, Path(tmp),
            java_jar=DEFAULT_JAVA_JAR,
            rust_bin=DEFAULT_RUST_BIN,
            timeout=TIMEOUT,
        )
        return compare(java_r, rust_r)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@given(raddose_config(budget=BUDGET))
@fuzz_settings
def test_outputs_match(cfg):
    """
    Java and Rust must produce numerically equivalent Summary.csv outputs.
    Hypothesis shrinks to the minimal failing Config on failure.
    """
    result = _run_and_compare(cfg)

    # Skip uninformative outcomes rather than failing on them
    assume(result.category not in SKIP_CATEGORIES)

    assert result.category not in FAIL_CATEGORIES, (
        f"\n{result.summary_line()}\n\n"
        f"--- input ---\n{render(cfg)}"
    )


@given(raddose_config(budget=BUDGET))
@fuzz_settings
def test_no_crashes(cfg):
    """
    Neither implementation should crash on any valid-looking input.
    Separated from test_outputs_match so shrinking focuses only on crashes.
    """
    result = _run_and_compare(cfg)
    assume(result.category not in SKIP_CATEGORIES)

    assert result.category not in {Category.JAVA_CRASH, Category.RUST_CRASH, Category.BOTH_CRASH}, (
        f"\n{result.summary_line()}\n\n"
        f"--- input ---\n{render(cfg)}"
    )


@given(raddose_config(budget=BUDGET))
@fuzz_settings
def test_no_nan_inf(cfg):
    """
    Neither implementation should produce NaN or Inf in summary metrics
    unless the other does too.
    """
    result = _run_and_compare(cfg)
    assume(result.category not in SKIP_CATEGORIES)

    assert result.category != Category.NAN_INF, (
        f"\n{result.summary_line()}\n\n"
        f"--- input ---\n{render(cfg)}"
    )
