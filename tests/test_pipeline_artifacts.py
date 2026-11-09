#!/usr/bin/env python3
"""Tests for portable TSV export and accession-level checkpoints."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PipelineArtifactTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="promotheus-artifacts.")
        self.directory = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_exports_tsv_and_records_complete_and_missing_checkpoints(self) -> None:
        output = self.directory / "output"
        results = output / "results"
        fulltab = results / "fulltab"
        fulltab.mkdir(parents=True)
        workbook = ROOT / "tests/fixtures/toy_expected/final_table.xlsx"
        shutil.copyfile(workbook, results / "final_table.xlsx")
        (fulltab / "GSE000001_full_tab.txt").write_text(
            "locus_tag\tmeans\tdiff_signif\nTOY_0001\t10\tNO\n",
            encoding="utf-8",
        )
        input_manifest = self.directory / "input_manifest.tsv"
        input_manifest.write_text(
            "kind\tfile\trows\tsha256\taccessions\n"
            "microarray\tmicro.xlsx\t1\tunused\tGSE000001\n"
            "rnaseq\trna.xlsx\t1\tunused\tGSE000002\n",
            encoding="utf-8",
        )

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/export_promotheus_tsv.py"),
                "--workbook",
                str(results / "final_table.xlsx"),
                "--output",
                str(results / "final_table.tsv"),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        checkpoint = self.directory / "checkpoint.tsv"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/build_checkpoint_manifest.py"),
                "--output-dir",
                str(output),
                "--input-manifest",
                str(input_manifest),
                "--manifest",
                str(checkpoint),
            ],
            check=True,
        )

        with (results / "final_table.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 5)
        self.assertIn("Sequence", rows[0])

        with checkpoint.open(encoding="utf-8", newline="") as handle:
            statuses = {
                row["accession"]: row["status"]
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(statuses, {"GSE000001": "complete", "GSE000002": "missing"})


if __name__ == "__main__":
    unittest.main()
