#!/usr/bin/env python3
"""Export the validated Promotheus workbook to a portable TSV table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from validate_promotheus_inputs import first_sheet_rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = first_sheet_rows(args.workbook)
    if len(rows) < 2:
        raise SystemExit(f"workbook has no data rows: {args.workbook}")

    width = max(len(row) for row in rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(row + [""] * (width - len(row)) for row in rows)

    print(f"Exported {len(rows) - 1} Promotheus rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
