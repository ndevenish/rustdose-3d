"""
Run both Java and Rust RADDOSE-3D on the same input file and collect results.
Both processes run in parallel. If one finishes before the deadline, the other
gets a straggler grace period so near-misses aren't misclassified as timeouts.
"""

import re
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).parent.parent
DEFAULT_JAVA_JAR = REPO_ROOT / "raddose3d.jar"
DEFAULT_RUST_BIN = REPO_ROOT / "raddose3d" / "target" / "release" / "raddose3d"
DEFAULT_TIMEOUT = 60.0    # seconds — shared deadline for both processes
DEFAULT_GRACE   = 30.0    # extra seconds given to the straggler after the first finishes
# Absolute cap on total wall time = timeout + grace, regardless of when the first finishes.


@dataclass
class RunResult:
    impl: str                   # "java" | "rust"
    exit_code: Optional[int]    # None = timed out
    wall_time: float            # seconds from launch to finish/kill
    stdout: str
    stderr: str
    output_dir: Path
    timed_out: bool = False
    error: str = ""             # harness-level error (binary not found, etc.)

    @property
    def _java_parse_error(self) -> bool:
        """Java exits 0 even on ANTLR parse errors — detect by stdout pattern."""
        combined = self.stdout + self.stderr
        return bool(re.search(r'InputException|Parser found \d+ errors', combined))

    @property
    def succeeded(self) -> bool:
        return (self.exit_code == 0 and not self.timed_out
                and not self.error and not self._java_parse_error)

    @property
    def crashed(self) -> bool:
        if self.timed_out or self.error:
            return False
        if self.exit_code != 0:
            return True
        # Java exits 0 on parse errors but runs a broken simulation anyway
        return self._java_parse_error

    def summary_csv_path(self) -> Optional[Path]:
        p = self.output_dir / f"{self.impl}-Summary.csv"
        return p if p.exists() else None


def run_both(
    input_path: Path,
    work_dir: Path,
    java_jar: Path = DEFAULT_JAVA_JAR,
    rust_bin: Path = DEFAULT_RUST_BIN,
    timeout: float = DEFAULT_TIMEOUT,
    grace: float = DEFAULT_GRACE,
) -> tuple[RunResult, RunResult]:
    """
    Run Java and Rust on input_path simultaneously.

    Both processes share a hard deadline of `timeout` seconds. If one finishes
    early the other gets up to `grace` additional seconds (capped so total wall
    time never exceeds `timeout + grace`).

    Each implementation writes its output files into its own subdirectory of
    work_dir so the two runs don't clobber each other.
    """
    java_out_dir = work_dir / "java_out"
    rust_out_dir = work_dir / "rust_out"
    java_out_dir.mkdir(parents=True, exist_ok=True)
    rust_out_dir.mkdir(parents=True, exist_ok=True)

    java_cmd = [
        "java", "-jar", str(java_jar),
        "-i", str(input_path),
        "-p", str(java_out_dir / "java-"),
    ]
    rust_cmd = [
        str(rust_bin),
        "-i", str(input_path),
        "-p", str(rust_out_dir / "rust-"),
    ]

    t0 = time.monotonic()
    hard_deadline   = t0 + timeout          # original shared deadline
    absolute_cap    = t0 + timeout + grace  # nothing runs past this

    # Shared state: when the first process finishes, it posts its finish time
    # so the second can compute its extended deadline.
    first_done_at: list[Optional[float]] = [None]
    first_done_lock = threading.Lock()

    java_proc = _launch(java_cmd, work_dir)
    rust_proc = _launch(rust_cmd, work_dir)

    results: list[Optional[RunResult]] = [None, None]

    def collect(idx: int, proc, impl: str, output_dir: Path) -> None:
        results[idx] = _collect_with_grace(
            proc, impl, output_dir,
            t0, hard_deadline, absolute_cap,
            first_done_at, first_done_lock,
        )

    t_java = threading.Thread(target=collect, args=(0, java_proc, "java", java_out_dir))
    t_rust = threading.Thread(target=collect, args=(1, rust_proc, "rust", rust_out_dir))
    t_java.start()
    t_rust.start()
    t_java.join()
    t_rust.join()

    return results[0], results[1]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_POLL_INTERVAL = 0.05   # seconds between poll() checks


def _launch(cmd: list[str], cwd: Path):
    try:
        return subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=str(cwd),
        )
    except OSError as e:
        return e


def _collect_with_grace(
    proc,
    impl: str,
    output_dir: Path,
    t0: float,
    hard_deadline: float,
    absolute_cap: float,
    first_done_at: list,
    first_done_lock: threading.Lock,
) -> RunResult:
    """
    Poll proc until it exits or its deadline expires.

    Deadline starts as hard_deadline.  Once the other process finishes
    (signalled via first_done_at), the deadline extends to
    min(first_done_at + grace, absolute_cap).  We never extend past
    absolute_cap so the total wall time is bounded.

    Uses poll() rather than communicate() so both processes can be
    monitored concurrently from separate threads without blocking each
    other.  RADDOSE-3D writes its actual output to files, so the pipe
    buffers (stdout/stderr) only carry progress messages and are
    unlikely to fill; communicate() is called once after the process
    exits to drain them safely.
    """
    if isinstance(proc, Exception):
        return RunResult(
            impl=impl, exit_code=None, wall_time=0.0,
            stdout="", stderr="", output_dir=output_dir,
            error=str(proc),
        )

    while True:
        now = time.monotonic()

        # Recompute deadline each iteration in case the partner just finished.
        # The straggler always gets at least until hard_deadline; if the partner
        # finished early it gets up to grace extra seconds beyond that, capped
        # by absolute_cap.
        with first_done_lock:
            partner_t = first_done_at[0]
        if partner_t is not None:
            my_deadline = min(max(hard_deadline, partner_t + (absolute_cap - hard_deadline)),
                              absolute_cap)
        else:
            my_deadline = hard_deadline

        if now >= my_deadline:
            proc.kill()
            try:
                stdout_b, stderr_b = proc.communicate(timeout=5)
            except Exception:
                stdout_b, stderr_b = b"", b""
            return RunResult(
                impl=impl, exit_code=None, wall_time=now - t0,
                stdout=stdout_b.decode(errors="replace"),
                stderr=stderr_b.decode(errors="replace"),
                output_dir=output_dir, timed_out=True,
            )

        rc = proc.poll()
        if rc is not None:
            # Process exited — drain pipes (safe: process is already done).
            try:
                stdout_b, stderr_b = proc.communicate(timeout=10)
            except Exception:
                stdout_b, stderr_b = b"", b""
            finish_time = time.monotonic()

            # Signal the partner it can extend its deadline.
            with first_done_lock:
                if first_done_at[0] is None:
                    first_done_at[0] = finish_time

            return RunResult(
                impl=impl, exit_code=rc, wall_time=finish_time - t0,
                stdout=stdout_b.decode(errors="replace"),
                stderr=stderr_b.decode(errors="replace"),
                output_dir=output_dir,
            )

        # Sleep until just before the next deadline check, waking early
        # if needed so we don't overshoot my_deadline.
        sleep_for = min(_POLL_INTERVAL, max(0.0, my_deadline - time.monotonic()))
        if sleep_for > 0:
            time.sleep(sleep_for)
