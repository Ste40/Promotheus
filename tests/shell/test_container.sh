#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
image="${PROMOTHEUS_IMAGE:-ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7}"

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: Docker is not installed" >&2
  exit 1
fi
if ! docker info >/dev/null 2>&1; then
  echo "ERROR: Docker daemon is not running" >&2
  exit 1
fi
if ! docker image inspect "$image" >/dev/null 2>&1; then
  echo "ERROR: image is not local; run: docker pull $image" >&2
  exit 1
fi

metadata="$(docker image inspect "$image" --format '{{.Os}}/{{.Architecture}} {{json .Config.Entrypoint}}')"
if [[ "$metadata" != 'linux/amd64 ["sh","/scripts/Introduction.sh"]' ]]; then
  echo "ERROR: unexpected container metadata: $metadata" >&2
  exit 1
fi

docker run --rm --platform linux/amd64 --entrypoint bash "$image" -lc '
  set -euo pipefail
  test -f /scripts/Microarray/main_pipeline.R
  test -f /scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R
  test -f /scripts/RNA-Seq/Rockhopper.jar
  test -x /OrthoFinder/orthofinder
  command -v Rscript
  command -v samtools
  command -v fasterq-dump
  command -v pysradb
  Rscript -e "stopifnot(all(vapply(c(\"openxlsx\",\"GEOquery\",\"Rsubread\",\"DESeq2\",\"sva\",\"Biostrings\"), requireNamespace, logical(1), quietly=TRUE)))"
'

files=(
  Scripts/Introduction.sh
  Scripts/Microarray/main_pipeline.R
)
container_files=(
  /scripts/Introduction.sh
  /scripts/Microarray/main_pipeline.R
)

for index in "${!files[@]}"; do
  local_hash="$(shasum -a 256 "$root_dir/${files[$index]}" | awk '{print $1}')"
  container_hash="$(docker run --rm --platform linux/amd64 --entrypoint sha256sum "$image" "${container_files[$index]}" | awk '{print $1}')"
  if [[ "$local_hash" != "$container_hash" ]]; then
    echo "ERROR: repository/container mismatch for ${files[$index]}" >&2
    exit 1
  fi
done

docker run --rm --platform linux/amd64 \
  -v "$root_dir/Scripts/Microarray/ortologhi_new.R:/scripts/Microarray/ortologhi_new.R:ro" \
  -v "$root_dir/Scripts/Microarray/dowload_file_assembly.R:/scripts/Microarray/dowload_file_assembly.R:ro" \
  -v "$root_dir/Scripts/Microarray/metadata_retrieval.R:/scripts/Microarray/metadata_retrieval.R:ro" \
  -v "$root_dir/Scripts/Microarray/intra_experiment_analysis.R:/scripts/Microarray/intra_experiment_analysis.R:ro" \
  -v "$root_dir/Scripts/RNA-Seq/Data_retrieval.R:/scripts/RNA-Seq/Data_retrieval.R:ro" \
  -v "$root_dir/Scripts/RNA-Seq/BAM_parsing.R:/scripts/RNA-Seq/BAM_parsing.R:ro" \
  -v "$root_dir/Scripts/RNA-Seq/NDE_detection.R:/scripts/RNA-Seq/NDE_detection.R:ro" \
  -v "$root_dir/Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:/scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:ro" \
  -v "$root_dir/Scripts/RNA-Seq/Operons_detection.R:/scripts/RNA-Seq/Operons_detection.R:ro" \
  -v "$root_dir/Scripts/Microarray/interexp_analysis.R:/scripts/Microarray/interexp_analysis.R:ro" \
  --entrypoint Rscript \
  "$image" \
  -e '
    scripts <- c(
      "/scripts/Microarray/dowload_file_assembly.R",
      "/scripts/RNA-Seq/Data_retrieval.R",
      "/scripts/Microarray/interexp_analysis.R"
    )
    expected <- paste0(
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/176/835/",
      "GCF_003176835.1_ASM317683v1/",
      "GCF_003176835.1_ASM317683v1_genomic.gbff.gz"
    )
    for (script in scripts) {
      environment <- new.env(parent = globalenv())
      sys.source(script, envir = environment)
      actual <- environment$refseq_file_location(
        paste0(
          "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/176/835/",
          "GCF_003176835.1_ASM317683v1/"
        ),
        "_genomic.gbff.gz"
      )$url
      stopifnot(identical(actual, expected))
    }
    for (script in c(
      "/scripts/RNA-Seq/BAM_parsing.R",
      "/scripts/RNA-Seq/NDE_detection.R",
      "/scripts/RNA-Seq/Operons_detection.R",
      "/scripts/Microarray/intra_experiment_analysis.R",
      "/scripts/Microarray/metadata_retrieval.R"
    )) {
      sys.source(script, envir = new.env(parent = globalenv()))
    }
    invisible(parse(file = "/scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R"))
  '

docker run --rm --platform linux/amd64 \
  -v "$root_dir/Scripts/Microarray/ortologhi_new.R:/workspace/Scripts/Microarray/ortologhi_new.R:ro" \
  -v "$root_dir/tests/test_orthology_empty_matches.R:/workspace/tests/test_orthology_empty_matches.R:ro" \
  --entrypoint Rscript \
  "$image" \
  /workspace/tests/test_orthology_empty_matches.R /workspace

echo "PASS: pinned Promotheus container, script overlay, and critical script identity"
