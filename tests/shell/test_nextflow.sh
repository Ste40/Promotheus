#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
temp_root="${TMPDIR:-/tmp}"
test_dir="$(mktemp -d "${temp_root%/}/promotheus-nextflow-test.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT

input_dir="$test_dir/input"
output_dir="$test_dir/output"
mkdir -p "$input_dir" "$output_dir"
cp "$root_dir/tests/fixtures/toy_input/GEO_Datasets_Search_Microarrays_Toy.xlsx" "$input_dir/"
cp "$root_dir/tests/fixtures/toy_input/GEO_Datasets_Search_RNA_Seq_Toy.xlsx" "$input_dir/"

cd "$test_dir"
NXF_OFFLINE=true nextflow run "$root_dir/main.nf" \
  -c "$root_dir/nextflow.config" \
  -profile test \
  -work-dir "$test_dir/work" \
  --input_dir "$input_dir" \
  --output_dir "$output_dir" \
  --asmcode ASM317683v1 \
  --threads 2 \
  --memory 1024 \
  --correction_de BH \
  --pvalue 0.01 \
  --operons_detection false \
  --amino_diff 5 \
  --similarity 40 \
  --operons_threshold 2 \
  --multitaxa false \
  -ansi-log false

test -s "$output_dir/results/final_table.xlsx"
test -s "$output_dir/results/final_table.tsv"
test -s "$output_dir/orthologous_genes.txt"
test -s "$output_dir/.promotheus-nextflow/input_validation.json"
test -s "$output_dir/.promotheus-nextflow/run_record.tsv"
test -s "$output_dir/.promotheus-nextflow/checkpoint_manifest.tsv"
test -s "$output_dir/.promotheus-nextflow/output_validation.json"

report_file="$(find "$output_dir/pipeline_info" -maxdepth 1 -name 'execution_report_*.html' -type f -print -quit)"
timeline_file="$(find "$output_dir/pipeline_info" -maxdepth 1 -name 'execution_timeline_*.html' -type f -print -quit)"
trace_file="$(find "$output_dir/pipeline_info" -maxdepth 1 -name 'execution_trace_*.txt' -type f -print -quit)"
dag_file="$(find "$output_dir/pipeline_info" -maxdepth 1 -name 'pipeline_dag_*.html' -type f -print -quit)"
test -s "$report_file"
test -s "$timeline_file"
test -s "$trace_file"
test -s "$dag_file"

cmp "$root_dir/tests/fixtures/toy_expected/final_table.xlsx" "$output_dir/results/final_table.xlsx"
grep -q $'engine\tmock' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'threads\t2' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'correction_de\tBH' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'asmcode\tASM317683v1' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'operons_detection\tfalse' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'multitaxa\tfalse' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'max_gse_parallel\t2' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'export_tsv\ttrue' "$output_dir/.promotheus-nextflow/run_record.tsv"
grep -q $'GSE000001\tcomplete' "$output_dir/.promotheus-nextflow/checkpoint_manifest.tsv"
grep -q $'GSE000002\tcomplete' "$output_dir/.promotheus-nextflow/checkpoint_manifest.tsv"
grep -q '"status": "pass"' "$output_dir/.promotheus-nextflow/input_validation.json"
grep -q '"status": "pass"' "$output_dir/.promotheus-nextflow/output_validation.json"
grep -q $'name\tstatus\texit\trealtime\t%cpu\tpeak_rss\tpeak_vmem' "$trace_file"
grep -q 'VALIDATE_INPUTS' "$trace_file"
grep -q 'ORTHOLOGY' "$trace_file"
grep -q 'MICROARRAY_GSE' "$trace_file"
grep -q 'RNA_DOWNLOAD_GSE' "$trace_file"
grep -q 'RNA_ALIGN_GSE' "$trace_file"
grep -q 'RNA_ANALYZE' "$trace_file"
grep -q 'INTEGRATE_RESULTS' "$trace_file"
grep -q 'VERIFY_OUTPUTS' "$trace_file"

echo "PASS: Nextflow mock end-to-end test"
