#!/usr/bin/env bash
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
for format in svg png; do
  dot -T"$format" "$repo/docs/diagrams/nextflow.dot" -o "$repo/docs/images/promotheus-nextflow-dag.$format"
done
