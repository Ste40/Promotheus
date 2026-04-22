#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

if ! command -v nextflow >/dev/null 2>&1; then
  echo "ERROR: nextflow binary not found in PATH" >&2
  exit 1
fi

INPUT_DIR="$TMP_DIR/input"
OUTPUT_DIR="$TMP_DIR/output"
mkdir -p "$INPUT_DIR" "$OUTPUT_DIR"

touch "$INPUT_DIR/Test_Microarray.xlsx"

cd "$ROOT_DIR"
nextflow run main.nf \
  -profile test \
  --input_dir "$INPUT_DIR" \
  --output_dir "$OUTPUT_DIR" \
  --asmcode ASM12345v1 \
  -ansi-log false

if ! find .nextflow -type f -name '.session*' >/dev/null 2>&1; then
  echo "ERROR: nextflow session metadata not generated" >&2
  exit 1
fi

echo "Nextflow smoke test completed successfully"
