#!/usr/bin/env python3
"""Repair and validate final_table.xlsx produced by the legacy container."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from xlsx_integrity import repair_xlsx_package, validate_xlsx_package


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True, type=Path)
    args = parser.parse_args()

    if not args.workbook.is_file():
        print(f"ERROR: workbook does not exist: {args.workbook}", file=sys.stderr)
        return 1
    try:
        repairs = repair_xlsx_package(args.workbook)
        errors = validate_xlsx_package(args.workbook)
    except Exception as exc:
        print(f"ERROR: cannot repair {args.workbook}: {exc}", file=sys.stderr)
        return 1
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    for repair in repairs:
        print(f"Repaired final_table.xlsx: {repair}")
    if not repairs:
        print("final_table.xlsx required no structural repairs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
