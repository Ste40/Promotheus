#!/usr/bin/env bash
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
python_bin="${PYTHON:-python3}"
"$python_bin" "$repo/docs/diagrams/render_nextflow_dag.py"
svg="$repo/docs/images/promotheus-nextflow-dag.svg"
png="$repo/docs/images/promotheus-nextflow-dag.png"
if command -v rsvg-convert >/dev/null 2>&1; then
  rsvg-convert "$svg" -o "$png"
elif "$python_bin" -c 'import resvg_py' 2>/dev/null; then
  "$python_bin" - "$svg" "$png" <<'PY'
import sys
from pathlib import Path
import resvg_py
Path(sys.argv[2]).write_bytes(resvg_py.svg_to_bytes(svg_path=sys.argv[1]))
PY
else
  echo 'SVG generated. PNG not updated: install rsvg-convert or resvg-py==0.3.2.' >&2
  exit 1
fi
