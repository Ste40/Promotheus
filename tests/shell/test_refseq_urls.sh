#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "SKIP: Rscript is unavailable; RefSeq URL unit test"
  exit 0
fi

Rscript "$root_dir/tests/test_refseq_urls.R" "$root_dir"
