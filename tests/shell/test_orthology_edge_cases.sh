#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "SKIP: Rscript is unavailable; orthology edge-case test"
  exit 0
fi

Rscript "$root_dir/tests/test_orthology_empty_matches.R" "$root_dir"
