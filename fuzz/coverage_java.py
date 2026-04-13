#!/usr/bin/env python3
"""
Offline coverage-driven corpus distillation for RADDOSE-3D Java.

Workflow
--------
1. Collect a pool of .txt inputs (corpus/seeds/, fuzz output, etc.).
2. Run this script on that pool to find which inputs together maximise
   Java bytecode coverage (line + branch + method) as measured by JaCoCo.

Algorithm
---------
For each input, run the Java JAR under the JaCoCo agent to collect a .exec
file, then use jacococli to export a coverage XML.  After all inputs are
processed a greedy set-cover pass selects the smallest subset that covers
the most covered lines.  Mirrors coverage.py (Rust) but uses JaCoCo instead
of LLVM instrumentation.

Usage
-----
    cd fuzz/
    uv run python coverage_java.py corpus/seeds/ [options]

    # Copy selected inputs and write a JSON report
    uv run python coverage_java.py corpus/seeds/ corpus/diffs/ \\
        --out-dir corpus/java_selected/ --report java_coverage.json

    # HTML source-annotated report (needs the Java source tree)
    uv run python coverage_java.py corpus/seeds/ \\
        --out-dir corpus/java_selected/ --html

Requirements
------------
- java on PATH (8+)
- JaCoCo jars are fetched automatically from Maven Central on first run and
  cached under fuzz/.jacoco/  (override with --jacoco-dir)
- raddose3d.jar at ../java/raddose3d.jar  (override with --java-jar)
"""

import argparse
import json
import os
import re as _re
import subprocess
import sys
import tempfile
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

# ---------------------------------------------------------------------------
# Defaults / constants
# ---------------------------------------------------------------------------

JACOCO_VERSION = "0.8.11"
REPO_ROOT = Path(__file__).parent.parent
JAVA_DIR = REPO_ROOT / "java"
FIXTURES_DIR = REPO_ROOT / "tests" / "fixtures"
DEFAULT_JAVA_JAR = JAVA_DIR / "raddose3d.jar"
DEFAULT_CLASSFILES = JAVA_DIR / "bin"
DEFAULT_JACOCO_DIR = Path(__file__).parent / ".jacoco"
DEFAULT_TIMEOUT = 0.0  # 0 = no timeout

_AGENT_URL = (
    f"https://repo1.maven.org/maven2/org/jacoco/org.jacoco.agent"
    f"/{JACOCO_VERSION}/org.jacoco.agent-{JACOCO_VERSION}-runtime.jar"
)
_CLI_URL = (
    f"https://repo1.maven.org/maven2/org/jacoco/org.jacoco.cli"
    f"/{JACOCO_VERSION}/org.jacoco.cli-{JACOCO_VERSION}-nodeps.jar"
)


# ---------------------------------------------------------------------------
# JaCoCo setup
# ---------------------------------------------------------------------------

def ensure_jacoco(jacoco_dir: Path) -> tuple[Path, Path]:
    """
    Return (agent_jar, cli_jar), downloading from Maven Central if absent.

    Files are cached under jacoco_dir so subsequent runs are instant.
    """
    jacoco_dir.mkdir(parents=True, exist_ok=True)
    agent = jacoco_dir / f"jacocoagent-{JACOCO_VERSION}.jar"
    cli   = jacoco_dir / f"jacococli-{JACOCO_VERSION}.jar"

    for jar, url in [(agent, _AGENT_URL), (cli, _CLI_URL)]:
        if not jar.exists():
            print(f"Downloading {jar.name} ...")
            try:
                urllib.request.urlretrieve(url, jar)
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to download {url}: {exc}\n"
                    "Supply --jacoco-dir pointing to a directory that already "
                    "contains the jars, or download them manually."
                ) from exc
            print(f"  -> {jar}")

    return agent, cli


# ---------------------------------------------------------------------------
# Input patching (shared with coverage.py)
# ---------------------------------------------------------------------------

def ensure_debug_classfiles(classfiles_dir: Path, java_dir: Path) -> Path:
    """
    Ensure java/bin/ contains debug-compiled .class files for line coverage.

    The pre-built raddose3d.jar is compiled with debug="false" so JaCoCo
    cannot map execution data back to source lines.  We run `ant build-debug`
    to compile with debug="true" into java/bin/, then pass that directory
    to jacococli report as --classfiles.  The JAR itself is still used at
    runtime — only the report step needs the debug classfiles.

    Skips the build if java/bin/se/ already exists (safe to re-run: ant is
    incremental, but we avoid the overhead when not needed).  Pass
    --rebuild-classfiles to force a rebuild.
    """
    sentinel = classfiles_dir / "se"
    if sentinel.is_dir():
        return classfiles_dir

    print(f"Compiling Java with debug info (ant build-debug in {java_dir}) ...")
    result = subprocess.run(
        ["ant", "build-debug"],
        cwd=java_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"ant build-debug failed:\n{result.stdout[-1000:]}\n{result.stderr[-500:]}"
        )
    if not sentinel.is_dir():
        raise RuntimeError(
            f"ant build-debug succeeded but {classfiles_dir} was not created"
        )
    print(f"  -> {classfiles_dir}")
    return classfiles_dir


def rewrite_fixture_paths(text: str) -> str:
    """
    Rewrite any path containing tests/fixtures/<file> to the absolute path
    on this machine (REPO_ROOT/tests/fixtures/<file>).

    Seeds were authored with hardcoded absolute paths like
      /Users/nickd/.../raddose3d/tests/fixtures/cube.obj
    which break when Java runs from a temp working directory.  We replace
    the whole path component with the correct local absolute path so that
    both Java and Rust can find the file regardless of where they're invoked.

    Handles all RADDOSE-3D keywords that take a file path argument:
      File, ModelFile, PDB, CIF, SeqFile
    and also bare paths on those lines (the path is the only token after the
    keyword).  Matching is case-insensitive on the keyword.
    """
    # Match: optional leading whitespace, a file-referencing keyword,
    # whitespace, then a path that contains "tests/fixtures/" somewhere.
    # Capture group 1 = everything up to the start of the path,
    # capture group 2 = the filename after the last slash.
    pattern = _re.compile(
        r'(?i)^(\s*(?:File|ModelFile|PDB|CIF|SeqFile)\s+)'  # keyword + spaces
        r'(?:\S+/)?tests/fixtures/(\S+)',                    # ...tests/fixtures/<name>
        _re.MULTILINE,
    )

    def _replace(m: _re.Match) -> str:
        keyword_part = m.group(1)
        filename     = m.group(2)
        local_path   = FIXTURES_DIR / filename
        return keyword_part + str(local_path)

    return pattern.sub(_replace, text)


def patch_sim_electrons(text: str, n: int) -> str:
    """
    Replace SimElectrons/SimPhotons value with n to keep MC/XFEL inputs
    fast without changing which Java code paths are exercised.
    """
    patched, count = _re.subn(
        r'(?i)^(\s*(?:SimElectrons|SimPhotons)\s+)[^\s#!]+',
        lambda m: m.group(1) + str(n),
        text,
        flags=_re.MULTILINE,
    )
    if count:
        return patched

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

def run_with_jacoco(
    java_jar: Path,
    classfiles: Path,
    jacoco_agent: Path,
    input_path: Path,
    exec_path: Path,
    timeout: float,
    patched_text: "str | None" = None,
) -> tuple[bool, float, str]:
    """
    Run RADDOSE-3D under the JaCoCo agent, writing coverage to exec_path.

    Uses `-cp classfiles:java_jar` rather than `-jar java_jar` so that the
    JaCoCo agent instruments the debug-compiled classes in classfiles/.  Those
    classes produce the same bytecode fingerprints that `jacococli report
    --classfiles classfiles/` expects, resolving the "class files do not match"
    warning that occurs when the pre-built JAR (no debug info) is used at
    runtime but the recompiled bin/ directory (with debug info) is used for
    reporting.

    Returns (succeeded, wall_time_seconds, stderr_snippet).
    """
    with tempfile.TemporaryDirectory(prefix="raddose_jcov_") as tmp:
        tmp_path = Path(tmp)

        if patched_text is not None:
            actual_input = tmp_path / input_path.name
            actual_input.write_text(patched_text)
        else:
            actual_input = input_path.resolve()

        jvmarg = (
            f"-javaagent:{jacoco_agent}"
            f"=destfile={exec_path}"
            f",append=false"
        )
        # classfiles/ first on classpath so the debug-compiled RADDOSE classes
        # shadow the ones bundled (without debug info) inside java_jar.
        # Library classes (antlr, commons-math3, selenium) still come from
        # java_jar since they don't exist in classfiles/.
        classpath = f"{classfiles}{os.pathsep}{java_jar}"
        out_prefix = str(tmp_path / "out-")
        # Always request output modules that are never triggered by default,
        # so their code is exercised on every run.  Use "-" as the destination
        # (WriterConsole/stdout) to avoid the prefix-prepend that
        # parseOutputDestinations applies to every filename argument.
        # FinalDoseStateRPreview is omitted: it writes megabytes of R code to
        # stdout for large crystals, perturbing JaCoCo coverage data.
        extra_outputs = [
            "-o", "ProgressEstimate:-",
            "-o", "FluencePerDoseHistCSV:-",
        ]
        cmd = [
            "java", jvmarg,
            "-cp", classpath,
            "se.raddo.raddose3D.RD3D",
            "-i", str(actual_input),
            "-p", out_prefix,
            *extra_outputs,
        ]
        t0 = time.monotonic()
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                timeout=timeout if timeout > 0 else None,
                cwd=tmp,
            )
            wall = time.monotonic() - t0
            succeeded = proc.returncode == 0
            # Java exits 0 even on parse errors — detect by stderr pattern
            combined = proc.stdout.decode(errors="replace") + proc.stderr.decode(errors="replace")
            if _re.search(r'InputException|Parser found \d+ errors', combined):
                succeeded = False
            stderr = combined[-400:]
            return succeeded, wall, stderr
        except subprocess.TimeoutExpired:
            return False, time.monotonic() - t0, "timed out"
        except Exception as exc:
            return False, 0.0, str(exc)


def exec_to_xml(
    classfiles: Path,
    jacoco_cli: Path,
    exec_path: Path,
    xml_path: Path,
) -> bool:
    """
    Convert a JaCoCo .exec file to an XML coverage report via jacococli.

    classfiles should point to the debug-compiled java/bin/ directory so that
    JaCoCo can map execution data back to source lines.
    """
    cmd = [
        "java", "-jar", str(jacoco_cli),
        "report", str(exec_path),
        "--classfiles", str(classfiles),
        "--xml", str(xml_path),
    ]
    result = subprocess.run(cmd, capture_output=True)
    return result.returncode == 0 and xml_path.exists()


def merge_execs(
    jacoco_cli: Path,
    exec_paths: list[Path],
    out_exec: Path,
) -> bool:
    """Merge multiple JaCoCo .exec files into one."""
    if not exec_paths:
        return False
    cmd = [
        "java", "-jar", str(jacoco_cli),
        "merge", *[str(p) for p in exec_paths],
        "--destfile", str(out_exec),
    ]
    result = subprocess.run(cmd, capture_output=True)
    return result.returncode == 0


def get_coverage_bitmap(xml_path: Path) -> frozenset[tuple]:
    """
    Parse a JaCoCo XML report and return a frozenset of covered RADDOSE lines.

    Each entry is (package_name, sourcefile_name, line_nr) for every line in a
    se/raddo/raddose3D package with at least one covered instruction (ci > 0).

    Requires the report to have been generated from debug-compiled classfiles
    (java/bin/ built with `ant build-debug`) so that JaCoCo can map execution
    data back to source lines.
    """
    covered = set()
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for package in root.iter("package"):
            pkg = package.get("name", "")
            if not pkg.startswith("se/raddo/raddose3D"):
                continue
            for sourcefile in package.iter("sourcefile"):
                sf = sourcefile.get("name", "")
                for line in sourcefile.iter("line"):
                    nr = int(line.get("nr", 0))
                    ci = int(line.get("ci", 0))  # covered instructions
                    if ci > 0:
                        covered.add((pkg, sf, nr))
    except Exception:
        pass
    return frozenset(covered)


def get_coverage_totals(xml_path: Path) -> dict:
    """
    Extract the aggregate <counter> elements from a JaCoCo XML report.

    Returns a dict keyed by counter type (instruction, line, branch, method,
    class) with covered/missed/count/percent sub-keys.
    """
    totals: dict[str, dict] = {}
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for counter in root.findall("counter"):
            ctype   = counter.get("type", "").lower()
            covered = int(counter.get("covered", 0))
            missed  = int(counter.get("missed",  0))
            total   = covered + missed
            totals[ctype] = {
                "covered": covered,
                "missed":  missed,
                "count":   total,
                "percent": 100.0 * covered / total if total else 0.0,
            }
    except Exception:
        pass
    return totals


# ---------------------------------------------------------------------------
# Greedy set-cover (identical algorithm to coverage.py)
# ---------------------------------------------------------------------------

def greedy_set_cover(
    inputs: list[Path],
    bitmaps: dict[Path, frozenset],
    wall_times: "dict[Path, float] | None" = None,
) -> list[Path]:
    """
    Iteratively pick the input covering the most uncovered lines.

    Tiebreak: prefer the faster input (lower wall time).
    Stops when no input adds new coverage.
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
        best = max(
            remaining,
            key=lambda p: (len(bitmaps[p] - covered), -times.get(p, 0.0)),
        )
        gain = len(bitmaps[best] - covered)
        if gain == 0:
            break
        covered |= bitmaps[best]
        selected.append(best)
        remaining.remove(best)

    return selected


# ---------------------------------------------------------------------------
# Input collection
# ---------------------------------------------------------------------------

def collect_inputs(sources: list[Path]) -> list[Path]:
    """Recursively find all .txt files under the given paths."""
    inputs: list[Path] = []
    for src in sources:
        if src.is_file() and src.suffix == ".txt":
            inputs.append(src)
        elif src.is_dir():
            inputs.extend(sorted(src.rglob("*.txt")))
    return inputs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Coverage-driven corpus distillation for RADDOSE-3D Java",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "sources", nargs="+", type=Path,
        help="Directory or .txt file(s) to analyse. Directories searched recursively.",
    )
    ap.add_argument(
        "--java-jar", type=Path, default=DEFAULT_JAVA_JAR,
        help=f"Path to raddose3d.jar (default: {DEFAULT_JAVA_JAR})",
    )
    ap.add_argument(
        "--jacoco-dir", type=Path, default=DEFAULT_JACOCO_DIR,
        help=f"Cache directory for JaCoCo jars (default: {DEFAULT_JACOCO_DIR})",
    )
    ap.add_argument(
        "--timeout", type=float, default=DEFAULT_TIMEOUT,
        help="Per-input timeout in seconds. 0 = no timeout (default). "
             "A killed process writes no .exec data, so set conservatively.",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=None,
        help="Copy selected inputs here. Default: print list only.",
    )
    ap.add_argument(
        "--report", type=Path, default=None,
        help="Write a JSON coverage report to this path.",
    )
    ap.add_argument(
        "--workers", "-j", type=int, default=1,
        help="Parallel workers (default: 1). Each gets its own .exec file slot.",
    )
    ap.add_argument(
        "--patch-sim-electrons", type=int, default=1000, metavar="N",
        help="Override SimElectrons/SimPhotons to N before running (default: 1000). "
             "Keeps MC/XFEL inputs fast without changing covered code paths. "
             "Set to 0 to disable.",
    )
    ap.add_argument(
        "--html", action="store_true",
        help="Generate an HTML coverage report (written to --out-dir/html/ or "
             "html/ in the current directory). Requires --source-dir for line "
             "annotations; without it the report shows only class-level data.",
    )
    ap.add_argument(
        "--source-dir", type=Path,
        default=REPO_ROOT / "java" / "src",
        metavar="DIR",
        help=f"Java source root for HTML report annotations "
             f"(default: {REPO_ROOT / 'java' / 'src'})",
    )
    ap.add_argument(
        "--classfiles", type=Path, default=DEFAULT_CLASSFILES,
        help=f"Directory of debug-compiled .class files for jacococli report "
             f"(default: {DEFAULT_CLASSFILES}). Built automatically with "
             f"`ant build-debug` if not present.",
    )
    ap.add_argument(
        "--rebuild-classfiles", action="store_true",
        help="Force `ant build-debug` even if java/bin/ already exists.",
    )
    ap.add_argument(
        "--save-exec", type=Path, default=None, metavar="PATH",
        help="Save the merged .exec file here for onward jacococli processing.",
    )
    args = ap.parse_args()

    import shutil as _shutil

    # Validate JAR
    if not args.java_jar.exists():
        print(f"ERROR: JAR not found: {args.java_jar}", file=sys.stderr)
        print("Supply --java-jar or build it: cd java && ant jar", file=sys.stderr)
        sys.exit(1)

    # Obtain JaCoCo
    try:
        jacoco_agent, jacoco_cli = ensure_jacoco(args.jacoco_dir)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    # Ensure debug-compiled classfiles for line coverage
    if args.rebuild_classfiles and (args.classfiles / "se").is_dir():
        import shutil as _rmshutil
        _rmshutil.rmtree(args.classfiles / "se")
    try:
        classfiles = ensure_debug_classfiles(args.classfiles, JAVA_DIR)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"JaCoCo agent: {jacoco_agent}")
    print(f"JaCoCo CLI:   {jacoco_cli}")
    print(f"Classfiles:   {classfiles}")
    print()

    # Collect inputs
    inputs = collect_inputs(args.sources)
    if not inputs:
        print("No .txt inputs found.", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(inputs)} input(s).")

    # Temp work dir for .exec and .xml files
    work_dir_obj = tempfile.TemporaryDirectory(prefix="raddose_jcov_work_")
    work_dir = Path(work_dir_obj.name)

    # -------------------------------------------------------------------------
    # Phase 1: run each input under JaCoCo, collect individual coverage bitmaps
    # -------------------------------------------------------------------------
    print(f"\nPhase 1: collecting coverage ({args.workers} worker(s), "
          f"timeout={args.timeout}s each)")
    print("-" * 60)

    bitmaps:        dict[Path, frozenset] = {}
    wall_times:     dict[Path, float]    = {}
    skipped:        list[Path]           = []
    per_input_exec: dict[Path, Path]     = {}

    def _process_one(slot: int, inp: Path) -> tuple:
        exec_path = work_dir / f"{slot:05d}.exec"
        xml_path  = work_dir / f"{slot:05d}.xml"

        text = rewrite_fixture_paths(inp.read_text())
        patched = patch_sim_electrons(text, args.patch_sim_electrons) if args.patch_sim_electrons else text
        succeeded, wall, stderr = run_with_jacoco(
            args.java_jar, classfiles, jacoco_agent, inp, exec_path, args.timeout,
            patched_text=patched,
        )

        bm = frozenset()
        if exec_path.exists():
            if exec_to_xml(classfiles, jacoco_cli, exec_path, xml_path):
                bm = get_coverage_bitmap(xml_path)

        return succeeded, wall, stderr, bm, exec_path

    if args.workers <= 1:
        for idx, inp in enumerate(inputs, 1):
            succeeded, wall, stderr, bm, exec_path = _process_one(idx, inp)
            wall_times[inp] = wall
            status = "OK" if succeeded else "FAIL"
            print(f"[{idx:4d}/{len(inputs)}] {status}  {wall:.1f}s  {len(bm):5,} lines  {inp.name}")
            if not succeeded and stderr:
                print(f"           stderr: {stderr[:120]}")
            bitmaps[inp] = bm
            if bm:
                per_input_exec[inp] = exec_path
            else:
                skipped.append(inp)
    else:
        import threading
        from concurrent.futures import ThreadPoolExecutor, as_completed

        slot_lock    = threading.Lock()
        slot_counter = [0]

        def _run_parallel(inp: Path):
            with slot_lock:
                slot = slot_counter[0]
                slot_counter[0] += 1
            return inp, *_process_one(slot + 1, inp)

        done = 0
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(_run_parallel, inp): inp for inp in inputs}
            for fut in as_completed(futures):
                inp, succeeded, wall, stderr, bm, exec_path = fut.result()
                done += 1
                status = "OK" if succeeded else "FAIL"
                print(f"[{done:4d}/{len(inputs)}] {status}  {wall:.1f}s  {len(bm):4,} methods  {inp.name}")
                if not succeeded and stderr:
                    print(f"           stderr: {stderr[:120]}")
                bitmaps[inp]    = bm
                wall_times[inp] = wall
                if bm:
                    per_input_exec[inp] = exec_path
                else:
                    skipped.append(inp)

    # -------------------------------------------------------------------------
    # Phase 2: greedy set-cover
    # -------------------------------------------------------------------------
    print()
    print("Phase 2: greedy set-cover")
    print("-" * 60)

    all_lines   = set().union(*bitmaps.values())
    total_lines = len(all_lines)

    usable   = [p for p in inputs if bitmaps[p]]
    selected = greedy_set_cover(usable, bitmaps, wall_times)

    selected_lines: set[tuple] = set()
    for p in selected:
        selected_lines |= bitmaps[p]

    print(f"Total RADDOSE lines across all inputs:   {total_lines:,}")
    print(f"Inputs with usable coverage:             {len(usable):,}")
    print(f"Inputs skipped (no exec / XML fail):     {len(skipped):,}")
    print(f"Selected by greedy set-cover:            {len(selected):,} input(s)")
    print(f"Lines covered by selection:              {len(selected_lines):,} "
          f"({100 * len(selected_lines) / max(total_lines, 1):.1f}%)")

    print()
    print("Selected inputs (in selection order):")
    cumulative: set[tuple] = set()
    for i, p in enumerate(selected, 1):
        gain = len(bitmaps[p] - cumulative)
        cumulative |= bitmaps[p]
        t = wall_times.get(p, 0.0)
        print(f"  {i:3d}.  +{gain:4,} lines  {t:5.1f}s  {p}")

    # -------------------------------------------------------------------------
    # Phase 3: merge selected execs → aggregate coverage report
    # -------------------------------------------------------------------------
    merged_exec = work_dir / "all_selected.exec"
    merged_xml  = work_dir / "all_selected.xml"

    selected_execs = [
        per_input_exec[p] for p in selected
        if p in per_input_exec and per_input_exec[p].exists()
    ]

    totals: dict = {}
    if selected_execs and merge_execs(jacoco_cli, selected_execs, merged_exec):
        if exec_to_xml(classfiles, jacoco_cli, merged_exec, merged_xml):
            totals = get_coverage_totals(merged_xml)

    if totals:
        def _pct(d: dict) -> str:
            return (f"{d.get('covered', 0):,}/{d.get('count', 0):,} "
                    f"({d.get('percent', 0.0):.1f}%)")
        print()
        print("RADDOSE-3D coverage (se.raddo.raddose3D.* only):")
        for key in ("instruction", "line", "branch", "method", "class"):
            if key in totals:
                print(f"  {key.capitalize():12s}: {_pct(totals[key])}")

    # -------------------------------------------------------------------------
    # Phase 4: copy selected inputs to --out-dir
    # -------------------------------------------------------------------------
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        for i, p in enumerate(selected, 1):
            dest = args.out_dir / f"{i:03d}_{p.name}"
            _shutil.copy2(p, dest)
        print(f"\nCopied {len(selected)} input(s) to {args.out_dir}")

    # Save merged exec (--save-exec, or auto into --out-dir)
    effective_exec_dest = args.save_exec
    if effective_exec_dest is None and args.out_dir is not None:
        effective_exec_dest = args.out_dir / "coverage.exec"
    if effective_exec_dest and merged_exec.exists():
        _shutil.copy2(merged_exec, effective_exec_dest)
        print(f"Exec file saved to {effective_exec_dest}")
        print(f"  Use with: java -jar {jacoco_cli} report {effective_exec_dest} "
              f"--classfiles {args.java_jar} --html <dir>")

    # -------------------------------------------------------------------------
    # Phase 5: JSON report
    # -------------------------------------------------------------------------
    if args.report:
        report_data: dict = {
            "jacoco_version":   JACOCO_VERSION,
            "inputs_analysed":  len(inputs),
            "inputs_skipped":   len(skipped),
            "inputs_selected":  len(selected),
            "total_lines":      total_lines,
            "covered_lines":    len(selected_lines),
            "coverage_pct":     round(100 * len(selected_lines) / max(total_lines, 1), 2),
            "totals":           totals,
            "selected": [
                {
                    "path":         str(p),
                    "wall_time":    wall_times.get(p, 0.0),
                    "lines_covered": len(bitmaps[p]),
                    "marginal_gain": len(bitmaps[p] - set().union(
                        *([bitmaps[q] for q in selected[:i]] or [set()])
                    )),
                }
                for i, p in enumerate(selected)
            ],
        }
        args.report.write_text(json.dumps(report_data, indent=2))
        print(f"Report written to {args.report}")

    # -------------------------------------------------------------------------
    # Phase 6: HTML report (optional)
    # -------------------------------------------------------------------------
    effective_html_dir = (
        args.out_dir / "html" if args.html and args.out_dir
        else Path("html") if args.html
        else None
    )

    if effective_html_dir and merged_exec.exists():
        effective_html_dir.mkdir(parents=True, exist_ok=True)
        html_cmd = [
            "java", "-jar", str(jacoco_cli),
            "report", str(merged_exec),
            "--classfiles", str(classfiles),
            "--html", str(effective_html_dir),
            "--name", "RADDOSE-3D Java coverage",
        ]
        if args.source_dir.exists():
            html_cmd += ["--sourcefiles", str(args.source_dir)]
        else:
            print(f"Note: --source-dir {args.source_dir} not found; "
                  "HTML report will lack source annotations.")
        html_result = subprocess.run(html_cmd, capture_output=True, text=True)
        if html_result.returncode == 0:
            print(f"HTML report written to {effective_html_dir / 'index.html'}")
        else:
            print(f"WARNING: HTML generation failed:\n{html_result.stderr[:400]}",
                  file=sys.stderr)

    work_dir_obj.cleanup()


if __name__ == "__main__":
    main()
