#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
temp_root="${TMPDIR:-/tmp}"
test_dir="$(mktemp -d "${temp_root%/}/promotheus-validation-test.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT

valid_dir="$test_dir/valid"
missing_rna_dir="$test_dir/missing-rna"
mkdir -p "$valid_dir" "$missing_rna_dir"
cp "$root_dir/tests/fixtures/toy_input/GEO_Datasets_Search_Microarrays_Toy.xlsx" "$valid_dir/"
cp "$root_dir/tests/fixtures/toy_input/GEO_Datasets_Search_RNA_Seq_Toy.xlsx" "$valid_dir/"
cp "$root_dir/tests/fixtures/toy_input/GEO_Datasets_Search_Microarrays_Toy.xlsx" "$missing_rna_dir/"

python3 "$root_dir/bin/validate_promotheus_inputs.py" \
  --input-dir "$valid_dir" \
  --asmcode ASM317683v1 \
  --manifest "$test_dir/valid.tsv" \
  --report "$test_dir/valid.json"
grep -q '"status": "pass"' "$test_dir/valid.json"

if python3 "$root_dir/bin/validate_promotheus_inputs.py" \
  --input-dir "$valid_dir" \
  --asmcode bad-code \
  --manifest "$test_dir/bad-asm.tsv" \
  --report "$test_dir/bad-asm.json"; then
  echo "ERROR: invalid ASM code was accepted" >&2
  exit 1
fi
grep -q 'invalid --asmcode' "$test_dir/bad-asm.json"

if python3 "$root_dir/bin/validate_promotheus_inputs.py" \
  --input-dir "$missing_rna_dir" \
  --asmcode ASM317683v1 \
  --manifest "$test_dir/missing.tsv" \
  --report "$test_dir/missing.json"; then
  echo "ERROR: missing RNA-seq workbook was accepted" >&2
  exit 1
fi
grep -q 'expected exactly one rnaseq workbook' "$test_dir/missing.json"

echo "PASS: input validation tests"
