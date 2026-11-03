#!/usr/bin/env bash
set -euo pipefail
root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/promotheus-modular-container.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT
python3 "$root_dir/bin/prepare_gse_inputs.py" \
  --input-dir "$root_dir/tests/fixtures/toy_input" --output-dir "$test_dir/shards"
mkdir "$test_dir/output"
docker run --rm --platform linux/amd64 --cpus 2 --memory 4g --network none \
  -v "$root_dir/Scripts:/scripts:ro" \
  -v "$root_dir/bin:/adapters:ro" \
  -v "$root_dir/tests/test_modular_container.R:/test_modular_container.R:ro" \
  -v "$test_dir/shards/GSE000002:/input:ro" \
  -v "$test_dir/output:/output" \
  --entrypoint Rscript \
  ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7 \
  /test_modular_container.R
