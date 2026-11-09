#!/usr/bin/env python3
"""Regression tests for malformed workbooks produced by legacy openxlsx."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from validate_promotheus_inputs import first_sheet_rows  # noqa: E402
from xlsx_integrity import (  # noqa: E402
    repair_xlsx_package,
    validate_xlsx_package,
)


def make_legacy_malformed(source: Path, destination: Path) -> None:
    with zipfile.ZipFile(source) as archive:
        entries = [(item, archive.read(item.filename)) for item in archive.infolist()]

    rewritten: dict[str, bytes] = {}
    worksheet_name = "xl/worksheets/sheet1.xml"
    entry_data = dict((item.filename, data) for item, data in entries)
    worksheet = entry_data[worksheet_name].decode("utf-8")
    if "<dimension " in worksheet:
        worksheet = worksheet.replace(worksheet[worksheet.index("<dimension ") : worksheet.index("/>", worksheet.index("<dimension ")) + 2], '<dimension ref="A1"/>', 1)
    else:
        worksheet_start = worksheet.index("<worksheet")
        worksheet_start_end = worksheet.index(">", worksheet_start)
        worksheet = worksheet[: worksheet_start_end + 1] + '<dimension ref="A1"/>' + worksheet[worksheet_start_end + 1 :]
    rewritten[worksheet_name] = worksheet.encode("utf-8")

    relationships_name = "xl/worksheets/_rels/sheet1.xml.rels"
    relationships = entry_data[relationships_name].decode("utf-8")
    relationships = relationships.replace(
        "</Relationships>",
        '<Relationship Id="rIdMissingDrawing" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/drawing" Target="../drawings/missing.xml"/></Relationships>',
    )
    rewritten[relationships_name] = relationships.encode("utf-8")

    content_types = entry_data["[Content_Types].xml"].decode("utf-8")
    content_types = content_types.replace(
        "</Types>",
        '<Override PartName="/xl/drawings/missing.xml" ContentType="application/vnd.openxmlformats-officedocument.drawing+xml"/></Types>',
    )
    rewritten["[Content_Types].xml"] = content_types.encode("utf-8")

    with zipfile.ZipFile(destination, "w") as output:
        for item, data in entries:
            output.writestr(item, rewritten.get(item.filename, data))


class XlsxIntegrityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="promotheus-xlsx-test.")
        self.directory = Path(self.temporary.name)
        self.fixture = ROOT / "tests/fixtures/toy_expected/final_table.xlsx"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_input_manifest(self, accessions: str = "GSE000001,GSE000002") -> Path:
        manifest = self.directory / "input_manifest.tsv"
        manifest.write_text(
            "kind\tfile\trows\tsha256\taccessions\n"
            f"microarray\tinput.xlsx\t1\tunused\t{accessions}\n",
            encoding="utf-8",
        )
        return manifest

    def run_output_verifier(
        self, output: Path, run_record: Path, input_manifest: Path
    ) -> tuple[subprocess.CompletedProcess[str], dict[str, object]]:
        report = self.directory / "report.json"
        checkpoint_manifest = self.directory / "checkpoint_manifest.tsv"
        accessions = input_manifest.read_text(encoding="utf-8").splitlines()[1].split("\t")[-1].split(",")
        checkpoint_manifest.write_text(
            "kind\taccession\tstatus\tartifacts\tsha256\n"
            + "".join(
                f"test\t{accession}\tcomplete\tfixture\tunused\n"
                for accession in accessions
            ),
            encoding="utf-8",
        )
        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/verify_promotheus_outputs.py"),
                "--output-dir",
                str(output),
                "--run-record",
                str(run_record),
                "--input-manifest",
                str(input_manifest),
                "--checkpoint-manifest",
                str(checkpoint_manifest),
                "--manifest",
                str(self.directory / "manifest.tsv"),
                "--report",
                str(report),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        return completed, json.loads(report.read_text(encoding="utf-8"))

    def test_detects_and_repairs_legacy_openxlsx_package(self) -> None:
        malformed = self.directory / "malformed.xlsx"
        make_legacy_malformed(self.fixture, malformed)

        errors = validate_xlsx_package(malformed)
        self.assertTrue(any("dangling relationship" in error for error in errors))
        self.assertTrue(any("content type references missing part" in error for error in errors))
        self.assertTrue(any("does not cover populated range" in error for error in errors))

        repairs = repair_xlsx_package(malformed)
        self.assertGreaterEqual(len(repairs), 3)
        self.assertEqual(validate_xlsx_package(malformed), [])
        self.assertEqual(first_sheet_rows(malformed), first_sheet_rows(self.fixture))

    def test_output_verifier_rejects_malformed_package(self) -> None:
        output = self.directory / "output"
        results = output / "results"
        results.mkdir(parents=True)
        malformed = results / "final_table.xlsx"
        make_legacy_malformed(self.fixture, malformed)
        (output / "orthologous_genes.txt").write_text("geneA\tgeneB\n", encoding="utf-8")
        run_record = self.directory / "run_record.tsv"
        run_record.write_text(
            "field\tvalue\nengine\tmock\nexport_tsv\tfalse\n", encoding="utf-8"
        )
        completed, validation = self.run_output_verifier(
            output, run_record, self.write_input_manifest()
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertEqual(validation["status"], "fail")
        self.assertTrue(any("invalid final_table.xlsx package" in error for error in validation["errors"]))

    def test_output_verifier_rejects_missing_expected_accession(self) -> None:
        output = self.directory / "coverage-output"
        results = output / "results"
        results.mkdir(parents=True)
        shutil.copyfile(self.fixture, results / "final_table.xlsx")
        (output / "orthologous_genes.txt").write_text("geneA\tgeneB\n", encoding="utf-8")
        run_record = self.directory / "coverage-run-record.tsv"
        run_record.write_text(
            "field\tvalue\nengine\tmock\nexport_tsv\tfalse\n", encoding="utf-8"
        )

        completed, validation = self.run_output_verifier(
            output,
            run_record,
            self.write_input_manifest("GSE000001,GSE000002,GSE999999"),
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertEqual(validation["status"], "fail")
        self.assertEqual(validation["missing_accessions"], ["GSE999999"])
        self.assertTrue(
            any("missing expression results" in error for error in validation["errors"])
        )


if __name__ == "__main__":
    unittest.main()
