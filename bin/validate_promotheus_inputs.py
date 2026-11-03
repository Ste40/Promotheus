#!/usr/bin/env python3
"""Validate the two Excel input tables consumed by the original Promotheus image."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree as ET


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
REQUIRED_COLUMNS = ("Accession", "taxon", "reference_genome")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def column_index(reference: str) -> int:
    letters = "".join(character for character in reference if character.isalpha())
    result = 0
    for character in letters.upper():
        result = result * 26 + ord(character) - ord("A") + 1
    return result - 1


def first_sheet_rows(path: Path) -> list[list[str]]:
    with zipfile.ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared_root.findall(f"{{{MAIN_NS}}}si"):
                shared_strings.append("".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t")))

        workbook_root = ET.fromstring(archive.read("xl/workbook.xml"))
        first_sheet = workbook_root.find(f".//{{{MAIN_NS}}}sheet")
        if first_sheet is None:
            raise ValueError("workbook contains no worksheets")
        relationship_id = first_sheet.attrib[f"{{{REL_NS}}}id"]
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = None
        for relationship in relationships.findall(f"{{{PKG_REL_NS}}}Relationship"):
            if relationship.attrib.get("Id") == relationship_id:
                target = relationship.attrib["Target"]
                break
        if target is None:
            raise ValueError("first worksheet relationship is missing")
        worksheet_name = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
        worksheet_name = str(PurePosixPath(worksheet_name))
        worksheet_root = ET.fromstring(archive.read(worksheet_name))

        rows: list[list[str]] = []
        for row in worksheet_root.findall(f".//{{{MAIN_NS}}}row"):
            values: dict[int, str] = {}
            for cell in row.findall(f"{{{MAIN_NS}}}c"):
                index = column_index(cell.attrib.get("r", "A1"))
                cell_type = cell.attrib.get("t")
                if cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{MAIN_NS}}}t"))
                else:
                    value_node = cell.find(f"{{{MAIN_NS}}}v")
                    value = value_node.text if value_node is not None and value_node.text is not None else ""
                    if cell_type == "s" and value:
                        value = shared_strings[int(value)]
                values[index] = value.strip()
            if values:
                width = max(values) + 1
                rows.append([values.get(index, "") for index in range(width)])
        return rows


def classify(path: Path) -> str | None:
    name = path.name.lower()
    if "microarray" in name:
        return "microarray"
    if "rna_seq" in name or "rna-seq" in name or "rnaseq" in name:
        return "rnaseq"
    return None


def validate_workbook(path: Path) -> tuple[dict[str, object], list[str]]:
    errors: list[str] = []
    kind = classify(path)
    if kind is None:
        errors.append(f"{path.name}: filename must identify Microarray or RNA_Seq data")
    try:
        rows = first_sheet_rows(path)
    except (OSError, KeyError, ValueError, zipfile.BadZipFile, ET.ParseError) as exc:
        return (
            {"file": path.name, "kind": kind, "rows": 0, "sha256": sha256(path)},
            [f"{path.name}: invalid .xlsx file ({exc})"],
        )

    if not rows:
        errors.append(f"{path.name}: first worksheet is empty")
        headers: list[str] = []
        data_rows: list[list[str]] = []
    else:
        headers = rows[0]
        data_rows = [row for row in rows[1:] if any(value for value in row)]

    missing = [column for column in REQUIRED_COLUMNS if column not in headers]
    if missing:
        errors.append(f"{path.name}: missing required column(s): {', '.join(missing)}")
    if not data_rows:
        errors.append(f"{path.name}: at least one data row is required")

    accessions: list[str] = []
    if not missing and data_rows:
        indices = {column: headers.index(column) for column in REQUIRED_COLUMNS}
        for number, row in enumerate(data_rows, start=2):
            values = {column: row[index] if index < len(row) else "" for column, index in indices.items()}
            accession = values["Accession"]
            accessions.append(accession)
            if not re.fullmatch(r"GSE\d+", accession):
                errors.append(f"{path.name}: row {number} has invalid Accession '{accession}'")
            if not values["taxon"]:
                errors.append(f"{path.name}: row {number} has an empty taxon")
            if not re.fullmatch(r"ASM\d+v\d+", values["reference_genome"]):
                errors.append(
                    f"{path.name}: row {number} has invalid reference_genome '{values['reference_genome']}'"
                )

    metadata: dict[str, object] = {
        "file": path.name,
        "kind": kind,
        "rows": len(data_rows),
        "columns": headers,
        "accessions": accessions,
        "sha256": sha256(path),
    }
    return metadata, errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--asmcode", required=True)
    parser.add_argument("--manifest", type=Path, default=Path("input_manifest.tsv"))
    parser.add_argument("--report", type=Path, default=Path("input_validation.json"))
    args = parser.parse_args()

    errors: list[str] = []
    if not re.fullmatch(r"ASM\d+v\d+", args.asmcode):
        errors.append(f"invalid --asmcode '{args.asmcode}'; expected ASM<digits>v<digits>")
    if not args.input_dir.is_dir():
        errors.append(f"input directory does not exist: {args.input_dir}")

    workbooks = sorted(args.input_dir.glob("*.xlsx")) if args.input_dir.is_dir() else []
    if not workbooks:
        errors.append("input directory contains no .xlsx files")

    metadata: list[dict[str, object]] = []
    for workbook in workbooks:
        workbook_metadata, workbook_errors = validate_workbook(workbook)
        metadata.append(workbook_metadata)
        errors.extend(workbook_errors)

    kinds = [item.get("kind") for item in metadata]
    for required_kind in ("microarray", "rnaseq"):
        count = kinds.count(required_kind)
        if count != 1:
            errors.append(f"expected exactly one {required_kind} workbook, found {count}")

    all_accessions = [accession for item in metadata for accession in item.get("accessions", [])]
    duplicates = sorted(
        accession
        for accession, count in Counter(all_accessions).items()
        if count > 1
    )
    if duplicates:
        errors.append(f"duplicate Accession values across inputs: {', '.join(duplicates)}")

    with args.manifest.open("w", encoding="utf-8") as handle:
        handle.write("kind\tfile\trows\tsha256\taccessions\n")
        for item in metadata:
            handle.write(
                f"{item.get('kind') or 'unknown'}\t{item['file']}\t{item['rows']}\t{item['sha256']}\t"
                f"{','.join(item.get('accessions', []))}\n"
            )

    report = {
        "status": "fail" if errors else "pass",
        "asmcode": args.asmcode,
        "workbooks": metadata,
        "errors": errors,
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"Validated {len(metadata)} Promotheus input workbooks ({len(all_accessions)} accessions).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
