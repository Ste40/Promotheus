#!/usr/bin/env python3
"""Split validated workbook rows into deterministic, isolated Nextflow inputs."""
import argparse
import csv
import hashlib
import json
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

from validate_promotheus_inputs import first_sheet_rows, classify


def write_xlsx(path, rows):
    """Small dependency-free workbook writer; preserve input values as strings."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    def column(n):
        result = ""
        while n:
            n, remainder = divmod(n - 1, 26)
            result = chr(65 + remainder) + result
        return result
    xml_rows = []
    for i, row in enumerate(rows, 1):
        cells = []
        for j, value in enumerate(row, 1):
            reference = f"{column(j)}{i}"
            if isinstance(value, (int, float)):
                cells.append(f'<c r="{reference}"><v>{value}</v></c>')
            else:
                # The pinned image's old openxlsx parser treats attributes on <t>
                # as literal cell text. Plain <t> works with both legacy and modern readers.
                cells.append(f'<c r="{reference}" t="inlineStr"><is><t>{escape(str(value))}</t></is></c>')
        xml_rows.append(f'<row r="{i}">{"".join(cells)}</row>')
    files = {
        "[Content_Types].xml": '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>',
        "_rels/.rels": '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        "xl/workbook.xml": '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>',
        "xl/_rels/workbook.xml.rels": '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>',
        "xl/worksheets/sheet1.xml": '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + ''.join(xml_rows) + '</sheetData></worksheet>',
    }
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, contents in files.items():
            entry = zipfile.ZipInfo(name, (2000, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' + contents)


def prepare(input_dir, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    organisms = {}
    for workbook in sorted(Path(input_dir).glob("*.xlsx")):
        kind = classify(workbook)
        rows = first_sheet_rows(workbook)
        header = rows[0]
        for ordinal, row in enumerate(rows[1:]):
            if not any(row):
                continue
            values = dict(zip(header, row))
            accession = values["Accession"]
            # The legacy R path replaces punctuation with spaces. Reject ambiguous paths.
            organism = re.sub(r"[^\w\s]|_", " ", values["taxon"])
            if kind == "rnaseq":
                identity = (values["taxon"], values["reference_genome"])
                if organism in organisms and organisms[organism] != identity:
                    raise ValueError("RNA-Seq organism has conflicting taxa or reference assemblies: " + organism)
                organisms[organism] = identity
            group = hashlib.sha256(organism.encode()).hexdigest()[:16]
            shard = output_dir / accession
            shard.mkdir()
            # Legacy metadata scripts consume the numeric GSE column, not Accession.
            if "GSE" in header:
                row[header.index("GSE")] = accession[3:]
            name = "Microarrays.xlsx" if kind == "microarray" else "RNA_Seq.xlsx"
            write_xlsx(shard / name, [header, row])
            descriptor = dict(kind=kind, accession=accession, group=group,
                              organism=organism, taxon=values["taxon"],
                              assembly=values["reference_genome"], ordinal=ordinal)
            (shard / "descriptor.json").write_text(json.dumps(descriptor, sort_keys=True) + "\n")
            records.append(descriptor)
    with (output_dir / "datasets.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, ["kind", "accession", "group", "ordinal"], delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    prepare(args.input_dir, args.output_dir)
