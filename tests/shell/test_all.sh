#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

for executable in python3 Rscript nextflow; do
  if ! command -v "$executable" >/dev/null 2>&1; then
    echo "ERROR: required test executable is unavailable: $executable" >&2
    exit 1
  fi
done
if ! Rscript -e 'quit(status = if (requireNamespace("openxlsx", quietly = TRUE)) 0 else 1)' >/dev/null 2>&1; then
  echo "ERROR: required R test package is unavailable: openxlsx" >&2
  exit 1
fi

PYTHONPYCACHEPREFIX="${TMPDIR:-/tmp}/promotheus-pycache" python3 -m py_compile \
  "$root_dir/bin/validate_promotheus_inputs.py" \
  "$root_dir/bin/verify_promotheus_outputs.py" \
  "$root_dir/bin/repair_promotheus_xlsx.py" \
  "$root_dir/bin/export_promotheus_tsv.py" \
  "$root_dir/bin/build_checkpoint_manifest.py" \
  "$root_dir/bin/xlsx_integrity.py" \
  "$root_dir/tests/test_xlsx_integrity.py" \
  "$root_dir/tests/test_pipeline_artifacts.py"

bash -n \
  "$root_dir/bin/run_promotheus_docker.sh" \
  "$root_dir/tests/bin/docker" \
  "$root_dir/tests/bin/run_mock_promotheus.sh" \
  "$root_dir/tests/shell/test_docker_runner.sh" \
  "$root_dir/tests/shell/test_orthology_edge_cases.sh" \
  "$root_dir/tests/shell/test_refseq_urls.sh" \
  "$root_dir/tests/shell/test_spreadsheet_compatibility.sh" \
  "$root_dir/tests/shell/test_validation.sh" \
  "$root_dir/tests/shell/test_nextflow.sh"

bash "$root_dir/tests/shell/test_validation.sh"
python3 "$root_dir/tests/test_xlsx_integrity.py"
python3 "$root_dir/tests/test_pipeline_artifacts.py"
python3 "$root_dir/tests/test_strict_output_validation.py"
python3 "$root_dir/tests/test_fastme_deterministic.py"
bash -n "$root_dir/bin/run_orthofinder.sh"
Rscript "$root_dir/tests/test_fastq_retry.R" "$root_dir"
Rscript "$root_dir/tests/test_rnaseq_performance.R" "$root_dir"
Rscript "$root_dir/tests/test_rnaseq_identity.R" "$root_dir"
Rscript "$root_dir/tests/test_operons_disabled.R" "$root_dir"
Rscript "$root_dir/tests/test_rockhopper_groups.R" "$root_dir"
Rscript "$root_dir/tests/test_interexp_performance.R" "$root_dir"
Rscript "$root_dir/tests/test_integration_samples.R" "$root_dir"
Rscript "$root_dir/tests/test_orthology_merge.R" "$root_dir"
bash "$root_dir/tests/shell/test_spreadsheet_compatibility.sh"
bash "$root_dir/tests/shell/test_docker_runner.sh"
bash "$root_dir/tests/shell/test_orthology_edge_cases.sh"
bash "$root_dir/tests/shell/test_refseq_urls.sh"
bash "$root_dir/tests/shell/test_nextflow.sh"
python3 "$root_dir/tests/test_modular_workflow.py"
python3 "$root_dir/tests/test_parallel_orthology.py"
python3 "$root_dir/tests/test_source_cache.py"
python3 "$root_dir/tests/test_resource_reservation.py"

echo "PASS: all Promotheus Nextflow tests"
