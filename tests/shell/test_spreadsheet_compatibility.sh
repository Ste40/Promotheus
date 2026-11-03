#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "SKIP: Rscript is unavailable; Python/Nextflow workbook validation still runs"
  exit 0
fi
if ! Rscript -e 'quit(status = if (requireNamespace("openxlsx", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1; then
  echo "SKIP: R package openxlsx is unavailable; Python/Nextflow workbook validation still runs"
  exit 0
fi

Rscript -e '
  micro <- openxlsx::read.xlsx(file.path(commandArgs(TRUE)[1], "tests/fixtures/toy_input/GEO_Datasets_Search_Microarrays_Toy.xlsx"))
  rna <- openxlsx::read.xlsx(file.path(commandArgs(TRUE)[1], "tests/fixtures/toy_input/GEO_Datasets_Search_RNA_Seq_Toy.xlsx"))
  out <- openxlsx::read.xlsx(file.path(commandArgs(TRUE)[1], "tests/fixtures/toy_expected/final_table.xlsx"))
  required <- c("Accession", "taxon", "reference_genome")
  stopifnot(
    is.data.frame(micro), is.data.frame(rna), is.data.frame(out),
    identical(nrow(micro), 1L), identical(nrow(rna), 1L), identical(nrow(out), 4L),
    all(required %in% names(micro)), all(required %in% names(rna)),
    all(nchar(out$Sequence) == 300)
  )
' "$root_dir"

echo "PASS: toy workbooks are compatible with Promotheus/openxlsx"
