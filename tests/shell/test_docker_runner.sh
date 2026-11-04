#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
temp_root="${TMPDIR:-/tmp}"
test_dir="$(mktemp -d "${temp_root%/}/promotheus-runner-test.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT

input_dir="$test_dir/input"
output_dir="$test_dir/output"
manifest="$test_dir/input_manifest.tsv"
record="$test_dir/run_record.tsv"
docker_log="$test_dir/docker_args.txt"
cache_dir="$test_dir/cache"
mkdir -p "$input_dir" "$output_dir" "$cache_dir"
input_dir="$(cd "$input_dir" && pwd -P)"
output_dir="$(cd "$output_dir" && pwd -P)"
printf 'kind\tfile\trows\tsha256\taccessions\n' > "$manifest"

PATH="$root_dir/tests/bin:$PATH" PROMOTHEUS_DOCKER_LOG="$docker_log" \
  bash "$root_dir/bin/run_promotheus_docker.sh" \
    --input-dir "$input_dir" \
    --output-dir "$output_dir" \
    --input-manifest "$manifest" \
    --record "$record" \
    --image 'ste40/promotheus@sha256:test' \
    --platform linux/amd64 \
    --threads 6 \
    --memory 8192 \
    --correction-de BH \
    --pvalue 0.01 \
    --operons-detection false \
    --amino-diff 4 \
    --similarity 55 \
    --asmcode ASM317683v1 \
    --operons-threshold 2 \
    --multitaxa true \
    --max-gse-parallel 3 \
    --cache-dir "$cache_dir" \
    --keep-cache true \
    --export-tsv true \
    --orthology-script "$root_dir/Scripts/Microarray/ortologhi_new.R" \
    --refseq-download-script "$root_dir/Scripts/Microarray/dowload_file_assembly.R" \
    --microarray-metadata-script "$root_dir/Scripts/Microarray/metadata_retrieval.R" \
    --microarray-intra-script "$root_dir/Scripts/Microarray/intra_experiment_analysis.R" \
    --rna-data-retrieval-script "$root_dir/Scripts/RNA-Seq/Data_retrieval.R" \
    --rna-bam-parsing-script "$root_dir/Scripts/RNA-Seq/BAM_parsing.R" \
    --rna-nde-detection-script "$root_dir/Scripts/RNA-Seq/NDE_detection.R" \
    --rna-pipeline-script "$root_dir/Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R" \
    --rna-operons-script "$root_dir/Scripts/RNA-Seq/Operons_detection.R" \
    --interexp-analysis-script "$root_dir/Scripts/Microarray/interexp_analysis.R"

for expected in \
  'run' \
  '--rm' \
  '-i' \
  '--platform' \
  'linux/amd64' \
  'THREADS=6' \
  'MEMORY=8192' \
  'CORRECTION_DE=BH' \
  'PVALUE=0.01' \
  'OPERONS_DETECTION=false' \
  'AMINO_DIFF=4' \
  'SIMILARITY=55' \
  'ASMcode=ASM317683v1' \
  'OPERONS_THRESHOLD=2' \
  'MULTITAXA=true' \
  'MAX_GSE_PARALLEL=3' \
  'PROMOTHEUS_KEEP_CACHE=true' \
  'PROMOTHEUS_EXPORT_TSV=true' \
  "$input_dir:/input:ro" \
  "$output_dir:/output" \
  "$root_dir/Scripts/Microarray/ortologhi_new.R:/scripts/Microarray/ortologhi_new.R:ro" \
  "$root_dir/Scripts/Microarray/dowload_file_assembly.R:/scripts/Microarray/dowload_file_assembly.R:ro" \
  "$root_dir/Scripts/Microarray/metadata_retrieval.R:/scripts/Microarray/metadata_retrieval.R:ro" \
  "$root_dir/Scripts/Microarray/intra_experiment_analysis.R:/scripts/Microarray/intra_experiment_analysis.R:ro" \
  "$root_dir/Scripts/RNA-Seq/Data_retrieval.R:/scripts/RNA-Seq/Data_retrieval.R:ro" \
  "$root_dir/Scripts/RNA-Seq/BAM_parsing.R:/scripts/RNA-Seq/BAM_parsing.R:ro" \
  "$root_dir/Scripts/RNA-Seq/NDE_detection.R:/scripts/RNA-Seq/NDE_detection.R:ro" \
  "$root_dir/Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:/scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:ro" \
  "$root_dir/Scripts/RNA-Seq/Operons_detection.R:/scripts/RNA-Seq/Operons_detection.R:ro" \
  "$root_dir/Scripts/Microarray/interexp_analysis.R:/scripts/Microarray/interexp_analysis.R:ro" \
  'ste40/promotheus@sha256:test'; do
  grep -Fxq -- "$expected" "$docker_log"
done

grep -q $'engine\tdocker' "$record"
grep -q $'asmcode\tASM317683v1' "$record"
grep -q $'max_gse_parallel\t3' "$record"
grep -q $'keep_cache\ttrue' "$record"
grep -q $'export_tsv\ttrue' "$record"
cache_work_dir="$(awk -F '\t' '$1 == "cache_work_dir" { print $2 }' "$record")"
grep -Fxq -- "$cache_work_dir:/output/tmp" "$docker_log"
echo "PASS: Docker runner compatibility contract"
