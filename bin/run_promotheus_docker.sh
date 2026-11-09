#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --input-dir DIR --output-dir DIR --input-manifest FILE --record FILE --asmcode ASM... [parameters]" >&2
}

input_dir=""
output_dir=""
input_manifest=""
record=""
image="ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7"
platform="linux/amd64"
threads="8"
memory="4000"
correction_de="none"
pvalue="0.05"
operons_detection="true"
amino_diff="3"
similarity="30"
asmcode=""
operons_threshold="3"
multitaxa="false"
max_gse_parallel="2"
cache_dir=""
keep_cache="true"
export_tsv="true"
orthology_script=""
refseq_download_script=""
microarray_metadata_script=""
microarray_intra_script=""
rna_data_retrieval_script=""
rna_metadata_script=""
rna_groups_script=""
rna_bam_parsing_script=""
rna_nde_detection_script=""
rna_pipeline_script=""
rna_operons_script=""
interexp_analysis_script=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) input_dir="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --input-manifest) input_manifest="$2"; shift 2 ;;
    --record) record="$2"; shift 2 ;;
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
    --orthology-script) orthology_script="$2"; shift 2 ;;
    --refseq-download-script) refseq_download_script="$2"; shift 2 ;;
    --microarray-metadata-script) microarray_metadata_script="$2"; shift 2 ;;
    --microarray-intra-script) microarray_intra_script="$2"; shift 2 ;;
    --rna-data-retrieval-script) rna_data_retrieval_script="$2"; shift 2 ;;
    --rna-metadata-script) rna_metadata_script="$2"; shift 2 ;;
    --rna-groups-script) rna_groups_script="$2"; shift 2 ;;
    --rna-bam-parsing-script) rna_bam_parsing_script="$2"; shift 2 ;;
    --rna-nde-detection-script) rna_nde_detection_script="$2"; shift 2 ;;
    --rna-pipeline-script) rna_pipeline_script="$2"; shift 2 ;;
    --rna-operons-script) rna_operons_script="$2"; shift 2 ;;
    --interexp-analysis-script) interexp_analysis_script="$2"; shift 2 ;;
    --fixture-output) shift 2 ;;
    *) usage; echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$input_dir" || -z "$output_dir" || -z "$input_manifest" || -z "$record" || -z "$asmcode" ]]; then
  usage
  exit 2
fi
if [[ ! "$max_gse_parallel" =~ ^[1-9][0-9]*$ ]]; then
  echo "--max-gse-parallel must be a positive integer" >&2
  exit 2
fi
for boolean_value in "$keep_cache" "$export_tsv"; do
  if [[ "$boolean_value" != "true" && "$boolean_value" != "false" ]]; then
    echo "--keep-cache and --export-tsv must be true or false" >&2
    exit 2
  fi
done

runner_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
orthology_script="${orthology_script:-$runner_root/Scripts/Microarray/ortologhi_new.R}"
refseq_download_script="${refseq_download_script:-$runner_root/Scripts/Microarray/dowload_file_assembly.R}"
microarray_metadata_script="${microarray_metadata_script:-$runner_root/Scripts/Microarray/metadata_retrieval.R}"
microarray_intra_script="${microarray_intra_script:-$runner_root/Scripts/Microarray/intra_experiment_analysis.R}"
rna_data_retrieval_script="${rna_data_retrieval_script:-$runner_root/Scripts/RNA-Seq/Data_retrieval.R}"
rna_metadata_script="${rna_metadata_script:-$runner_root/Scripts/RNA-Seq/MetaData_retrieval.R}"
rna_groups_script="${rna_groups_script:-$runner_root/Scripts/RNA-Seq/Groups_retrieval.R}"
rna_bam_parsing_script="${rna_bam_parsing_script:-$runner_root/Scripts/RNA-Seq/BAM_parsing.R}"
rna_nde_detection_script="${rna_nde_detection_script:-$runner_root/Scripts/RNA-Seq/NDE_detection.R}"
rna_pipeline_script="${rna_pipeline_script:-$runner_root/Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R}"
rna_operons_script="${rna_operons_script:-$runner_root/Scripts/RNA-Seq/Operons_detection.R}"
interexp_analysis_script="${interexp_analysis_script:-$runner_root/Scripts/Microarray/interexp_analysis.R}"
for script_path in "$orthology_script" "$refseq_download_script" "$microarray_metadata_script" "$microarray_intra_script" "$rna_data_retrieval_script" "$rna_metadata_script" "$rna_groups_script" "$rna_bam_parsing_script" "$rna_nde_detection_script" "$rna_pipeline_script" "$rna_operons_script" "$interexp_analysis_script"; do
  if [[ ! -f "$script_path" ]]; then
    echo "Promotheus overlay script does not exist: $script_path" >&2
    exit 2
  fi
done
orthology_script="$(realpath "$orthology_script")"
refseq_download_script="$(realpath "$refseq_download_script")"
microarray_metadata_script="$(realpath "$microarray_metadata_script")"
microarray_intra_script="$(realpath "$microarray_intra_script")"
rna_data_retrieval_script="$(realpath "$rna_data_retrieval_script")"
rna_metadata_script="$(realpath "$rna_metadata_script")"
rna_groups_script="$(realpath "$rna_groups_script")"
rna_bam_parsing_script="$(realpath "$rna_bam_parsing_script")"
rna_nde_detection_script="$(realpath "$rna_nde_detection_script")"
rna_pipeline_script="$(realpath "$rna_pipeline_script")"
rna_operons_script="$(realpath "$rna_operons_script")"
interexp_analysis_script="$(realpath "$interexp_analysis_script")"
if [[ ! -d "$input_dir" ]]; then
  echo "Input directory does not exist: $input_dir" >&2
  exit 2
fi
if [[ ! -f "$input_manifest" ]]; then
  echo "Input manifest does not exist: $input_manifest" >&2
  exit 2
fi
if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is required for the standard profile." >&2
  exit 127
fi
if ! docker info >/dev/null 2>&1; then
  echo "Docker is installed but its daemon is not running." >&2
  exit 1
fi

mkdir -p "$output_dir"
input_dir="$(cd "$input_dir" && pwd -P)"
output_dir="$(cd "$output_dir" && pwd -P)"

cache_mount=()
cache_key="disabled"
cache_work_dir=""
if [[ -n "$cache_dir" ]]; then
  mkdir -p "$cache_dir"
  cache_dir="$(cd "$cache_dir" && pwd -P)"
  overlay_hashes=""
  if command -v sha256sum >/dev/null 2>&1; then
    manifest_hash="$(sha256sum "$input_manifest" | awk '{print $1}')"
    for script_path in "$orthology_script" "$refseq_download_script" "$microarray_metadata_script" "$microarray_intra_script" "$rna_data_retrieval_script" "$rna_bam_parsing_script" "$rna_nde_detection_script" "$rna_pipeline_script" "$rna_operons_script" "$interexp_analysis_script"; do
      overlay_hashes+="$(sha256sum "$script_path" | awk '{print $1}')"
    done
    cache_key="$(printf '%s' "$manifest_hash|$overlay_hashes|$asmcode|$amino_diff|$similarity|$correction_de|$pvalue" | sha256sum | awk '{print substr($1,1,20)}')"
  else
    manifest_hash="$(shasum -a 256 "$input_manifest" | awk '{print $1}')"
    for script_path in "$orthology_script" "$refseq_download_script" "$microarray_metadata_script" "$microarray_intra_script" "$rna_data_retrieval_script" "$rna_bam_parsing_script" "$rna_nde_detection_script" "$rna_pipeline_script" "$rna_operons_script" "$interexp_analysis_script"; do
      overlay_hashes+="$(shasum -a 256 "$script_path" | awk '{print $1}')"
    done
    cache_key="$(printf '%s' "$manifest_hash|$overlay_hashes|$asmcode|$amino_diff|$similarity|$correction_de|$pvalue" | shasum -a 256 | awk '{print substr($1,1,20)}')"
  fi
  cache_work_dir="$cache_dir/$cache_key/tmp"
  mkdir -p "$cache_work_dir"
  cache_mount=(-v "$cache_work_dir:/output/tmp")
fi

docker_args=(run --rm -i)
if [[ -n "$platform" && "$platform" != "auto" ]]; then
  docker_args+=(--platform "$platform")
fi
docker_args+=(
  -e "THREADS=$threads"
  -e "MEMORY=$memory"
  -e "CORRECTION_DE=$correction_de"
  -e "PVALUE=$pvalue"
  -e "OPERONS_DETECTION=$operons_detection"
  -e "AMINO_DIFF=$amino_diff"
  -e "SIMILARITY=$similarity"
  -e "ASMcode=$asmcode"
  -e "OPERONS_THRESHOLD=$operons_threshold"
  -e "MULTITAXA=$multitaxa"
  -e "MAX_GSE_PARALLEL=$max_gse_parallel"
  -e "PROMOTHEUS_KEEP_CACHE=$keep_cache"
  -e "PROMOTHEUS_EXPORT_TSV=$export_tsv"
  -v "$input_dir:/input:ro"
  -v "$output_dir:/output"
  "${cache_mount[@]}"
  -v "$orthology_script:/scripts/Microarray/ortologhi_new.R:ro"
  -v "$refseq_download_script:/scripts/Microarray/dowload_file_assembly.R:ro"
  -v "$microarray_metadata_script:/scripts/Microarray/metadata_retrieval.R:ro"
  -v "$microarray_intra_script:/scripts/Microarray/intra_experiment_analysis.R:ro"
  -v "$rna_data_retrieval_script:/scripts/RNA-Seq/Data_retrieval.R:ro"
  -v "$rna_bam_parsing_script:/scripts/RNA-Seq/BAM_parsing.R:ro"
  -v "$rna_metadata_script:/scripts/RNA-Seq/MetaData_retrieval.R:ro"
  -v "$rna_groups_script:/scripts/RNA-Seq/Groups_retrieval.R:ro"
  -v "$rna_nde_detection_script:/scripts/RNA-Seq/NDE_detection.R:ro"
  -v "$rna_pipeline_script:/scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:ro"
  -v "$rna_operons_script:/scripts/RNA-Seq/Operons_detection.R:ro"
  -v "$interexp_analysis_script:/scripts/Microarray/interexp_analysis.R:ro"
  "$image"
)

# Both legacy branches can request confirmation after inferring experimental
# groups. Supplying one line for each prompt preserves the defaults while
# keeping scheduled Nextflow runs non-interactive.
printf '\n\n' | docker "${docker_args[@]}"

{
  printf 'field\tvalue\n'
  printf 'engine\tdocker\n'
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
  printf 'cache_key\t%s\n' "$cache_key"
  printf 'cache_work_dir\t%s\n' "$cache_work_dir"
  printf 'keep_cache\t%s\n' "$keep_cache"
  printf 'export_tsv\t%s\n' "$export_tsv"
} > "$record"
