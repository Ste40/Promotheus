#!/usr/bin/env python3
"""Verify the principal artifacts produced by Promotheus."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

from validate_promotheus_inputs import first_sheet_rows
from xlsx_integrity import validate_xlsx_package


REQUIRED_FINAL_COLUMNS = (
    "Gene_name",
    "GeneID",
    "Mean_Expression",
    "CV",
    "Mean_Percentile",
    "Percentile_SD",
    "DE",
    "Sequence",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--run-record", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--checkpoint-manifest", required=True, type=Path)
    parser.add_argument("--manifest", type=Path, default=Path("output_manifest.tsv"))
    parser.add_argument("--report", type=Path, default=Path("output_validation.json"))
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []
    final_table = args.output_dir / "results" / "final_table.xlsx"
    final_tsv = args.output_dir / "results" / "final_table.tsv"
    orthologues = args.output_dir / "orthologous_genes.txt"

    run_settings: dict[str, str] = {}
    if not args.run_record.is_file() or args.run_record.stat().st_size == 0:
        errors.append("run record is missing or empty")
    else:
        try:
            with args.run_record.open(encoding="utf-8", newline="") as handle:
                run_settings = {
                    row["field"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
        except (OSError, csv.Error, KeyError) as exc:
            errors.append(f"cannot read run record: {exc}")
    expected_accessions: set[str] = set()
    if not args.input_manifest.is_file() or args.input_manifest.stat().st_size == 0:
        errors.append("input manifest is missing or empty")
    else:
        try:
            with args.input_manifest.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    expected_accessions.update(
                        accession.strip()
                        for accession in row.get("accessions", "").split(",")
                        if accession.strip()
                    )
        except (OSError, csv.Error) as exc:
            errors.append(f"cannot read input manifest: {exc}")
        if not expected_accessions:
            errors.append("input manifest contains no expected GSE accessions")
    if not final_table.is_file():
        errors.append(f"missing final result: {final_table}")
    if not orthologues.is_file() or orthologues.stat().st_size == 0:
        errors.append(f"missing or empty orthology table: {orthologues}")
    checkpoint_status: dict[str, str] = {}
    if not args.checkpoint_manifest.is_file() or args.checkpoint_manifest.stat().st_size == 0:
        errors.append("checkpoint manifest is missing or empty")
    else:
        try:
            with args.checkpoint_manifest.open(encoding="utf-8", newline="") as handle:
                checkpoint_status = {
                    row["accession"]: row["status"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
        except (OSError, csv.Error, KeyError) as exc:
            errors.append(f"cannot read checkpoint manifest: {exc}")

    headers: list[str] = []
    data_rows: list[list[str]] = []
    expression_accessions: set[str] = set()
    percentile_accessions: set[str] = set()
    if final_table.is_file():
        package_errors = validate_xlsx_package(final_table)
        errors.extend(f"invalid final_table.xlsx package: {error}" for error in package_errors)
        try:
            rows = first_sheet_rows(final_table)
            headers = rows[0] if rows else []
            data_rows = [row for row in rows[1:] if any(value for value in row)]
        except Exception as exc:  # The report should preserve the concrete workbook error.
            errors.append(f"cannot read final_table.xlsx: {exc}")
        missing = [column for column in REQUIRED_FINAL_COLUMNS if column not in headers]
        if missing:
            errors.append(f"final_table.xlsx is missing column(s): {', '.join(missing)}")
        if not any(column.startswith("Expression_") for column in headers):
            errors.append("final_table.xlsx has no Expression_<GSE> column")
        if not any(column.startswith("Percentile_") and column != "Percentile_SD" for column in headers):
            errors.append("final_table.xlsx has no Percentile_<GSE> column")
        for header in headers:
            expression_match = re.match(r"^Expression_(GSE\d+)(?:\.|$)", header)
            percentile_match = re.match(r"^Percentile_(GSE\d+)(?:\.|$)", header)
            if expression_match:
                expression_accessions.add(expression_match.group(1))
            if percentile_match:
                percentile_accessions.add(percentile_match.group(1))
        missing_expression = sorted(expected_accessions - expression_accessions)
        missing_percentile = sorted(expected_accessions - percentile_accessions)
        if missing_expression:
            errors.append(
                "final_table.xlsx is missing expression results for expected GSE accession(s): "
                + ", ".join(missing_expression)
            )
        if missing_percentile:
            errors.append(
                "final_table.xlsx is missing percentile results for expected GSE accession(s): "
                + ", ".join(missing_percentile)
            )
        if not data_rows:
            errors.append("final_table.xlsx contains no genes")
        if "Sequence" in headers and data_rows:
            sequence_index = headers.index("Sequence")
            invalid_lengths = 0
            for row in data_rows:
                sequence = row[sequence_index] if sequence_index < len(row) else ""
                if len(sequence) != 300:
                    invalid_lengths += 1
            if invalid_lengths:
                warnings.append(f"{invalid_lengths} promoter sequence(s) are not 300 bp")

    if run_settings.get("export_tsv", "true") == "true":
        if not final_tsv.is_file() or final_tsv.stat().st_size == 0:
            errors.append(f"missing or empty TSV result: {final_tsv}")
        else:
            try:
                with final_tsv.open(encoding="utf-8", newline="") as handle:
                    tsv_rows = list(csv.reader(handle, delimiter="\t"))
                if not tsv_rows or tsv_rows[0] != headers:
                    errors.append("final_table.tsv headers differ from final_table.xlsx")
                if len(tsv_rows[1:]) != len(data_rows):
                    errors.append("final_table.tsv row count differs from final_table.xlsx")
            except (OSError, csv.Error) as exc:
                errors.append(f"cannot read final_table.tsv: {exc}")

    missing_checkpoints = sorted(
        accession for accession in expected_accessions
        if checkpoint_status.get(accession) != "complete"
    )
    if missing_checkpoints:
        warnings.append(
            "checkpoint manifest has no reusable full-tab for: "
            + ", ".join(missing_checkpoints)
        )

    artifacts = [
        path
        for path in (
            final_table,
            final_tsv,
            orthologues,
            args.run_record,
            args.checkpoint_manifest,
        )
        if path.is_file()
    ]
    with args.manifest.open("w", encoding="utf-8") as handle:
        handle.write("artifact\tsize_bytes\tsha256\n")
        for artifact in artifacts:
            handle.write(f"{artifact}\t{artifact.stat().st_size}\t{sha256(artifact)}\n")

    report = {
        "status": "fail" if errors else "pass",
        "output_dir": str(args.output_dir),
        "final_table": str(final_table),
        "genes": len(data_rows),
        "columns": headers,
        "expected_accessions": sorted(expected_accessions),
        "expression_accessions": sorted(expression_accessions),
        "percentile_accessions": sorted(percentile_accessions),
        "missing_accessions": sorted(
            expected_accessions - (expression_accessions & percentile_accessions)
        ),
        "checkpoints": checkpoint_status,
        "errors": errors,
        "warnings": warnings,
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"Verified Promotheus output: {len(data_rows)} genes, {len(headers)} columns.")
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
