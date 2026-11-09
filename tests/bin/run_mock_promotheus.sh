#!/usr/bin/env bash
set -euo pipefail

input_dir=""
output_dir=""
input_manifest=""
record=""
fixture_output=""
image=""
platform=""
threads=""
memory=""
correction_de=""
pvalue=""
operons_detection=""
amino_diff=""
similarity=""
asmcode=""
operons_threshold=""
multitaxa=""
max_gse_parallel=""
cache_dir=""
keep_cache=""
export_tsv=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) input_dir="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --input-manifest) input_manifest="$2"; shift 2 ;;
    --record) record="$2"; shift 2 ;;
    --fixture-output) fixture_output="$2"; shift 2 ;;
    --image) image="$2"; shift 2 ;;
    --platform) platform="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --memory) memory="$2"; shift 2 ;;
    --correction-de) correction_de="$2"; shift 2 ;;
    --pvalue) pvalue="$2"; shift 2 ;;
    --operons-detection) operons_detection="$2"; shift 2 ;;
    --amino-diff) amino_diff="$2"; shift 2 ;;
    --similarity) similarity="$2"; shift 2 ;;
    --asmcode) asmcode="$2"; shift 2 ;;
    --operons-threshold) operons_threshold="$2"; shift 2 ;;
    --multitaxa) multitaxa="$2"; shift 2 ;;
    --max-gse-parallel) max_gse_parallel="$2"; shift 2 ;;
    --cache-dir) cache_dir="$2"; shift 2 ;;
    --keep-cache) keep_cache="$2"; shift 2 ;;
    --export-tsv) export_tsv="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

for required in input_dir output_dir input_manifest record fixture_output asmcode; do
  if [[ -z "${!required}" ]]; then
    echo "Missing mock runner argument: $required" >&2
    exit 2
  fi
done

mkdir -p "$output_dir/results/fulltab" "$output_dir/RNA_Seq/results"
cp "$fixture_output" "$output_dir/results/final_table.xlsx"
printf 'geneA\tTOY_0001\n geneB\tTOY_0002\n' > "$output_dir/orthologous_genes.txt"
printf 'locus_tag\tmeans\tdiff_signif\nTOY_0001\t10\tNO\n' > "$output_dir/results/fulltab/GSE000001_full_tab.txt"
printf 'locus_tag\tmeans\tdiff_signif\nTOY_0001\t11\tYES\n' > "$output_dir/RNA_Seq/results/GSE000002_Toy_full_tab.txt"

{
  printf 'field\tvalue\n'
  printf 'engine\tmock\n'
  printf 'image\t%s\n' "$image"
  printf 'platform\t%s\n' "$platform"
  printf 'input_dir\t%s\n' "$input_dir"
  printf 'output_dir\t%s\n' "$output_dir"
  printf 'threads\t%s\n' "$threads"
  printf 'memory\t%s\n' "$memory"
  printf 'correction_de\t%s\n' "$correction_de"
  printf 'pvalue\t%s\n' "$pvalue"
  printf 'operons_detection\t%s\n' "$operons_detection"
  printf 'amino_diff\t%s\n' "$amino_diff"
  printf 'similarity\t%s\n' "$similarity"
  printf 'asmcode\t%s\n' "$asmcode"
  printf 'operons_threshold\t%s\n' "$operons_threshold"
  printf 'multitaxa\t%s\n' "$multitaxa"
  printf 'max_gse_parallel\t%s\n' "$max_gse_parallel"
  printf 'cache_dir\t%s\n' "$cache_dir"
  printf 'keep_cache\t%s\n' "$keep_cache"
  printf 'export_tsv\t%s\n' "$export_tsv"
} > "$record"
