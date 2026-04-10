#!/usr/bin/env bash
# Build the coverage-instrumented binary.
# Output: target/coverage/release/raddose3d  (separate from target/release/)
set -euo pipefail
cd "$(dirname "$0")"
RUSTFLAGS="-C instrument-coverage" cargo build --release --target-dir target/coverage "$@"
echo
echo "Coverage binary: $(pwd)/target/coverage/release/raddose3d"
