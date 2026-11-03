#!/usr/bin/env python3
"""Content changes must invalidate all scientific tasks despite preserved metadata."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SourceCacheTests(unittest.TestCase):
    def test_content_changes_and_unchanged_resume(self):
        with tempfile.TemporaryDirectory(prefix="promotheus-source-cache-") as temp:
            folder = Path(temp)
            project = folder / "project"
            project.mkdir()
            for name in ("main.nf", "nextflow.config"):
                shutil.copy2(ROOT / name, project / name)
            for name in ("bin", "Scripts", "tests/fixtures"):
                shutil.copytree(ROOT / name, project / name, ignore=shutil.ignore_patterns("__pycache__"))
            command = ["nextflow", "run", str(project / "main.nf"), "-profile", "test",
                       "--input_dir", str(project / "tests/fixtures/toy_input"),
                       "--output_dir", str(folder / "output"), "--asmcode", "ASM317683v1",
                       "--threads", "4", "--memory", "2048", "--operons_detection", "false",
                       "-ansi-log", "false"]

            def execute(name, resume=False):
                trace = folder / (name + ".tsv")
                result = subprocess.run(command + ["-with-trace", str(trace)] + (["-resume"] if resume else []),
                                        cwd=folder, env=dict(os.environ, NXF_OFFLINE="true"),
                                        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
                self.assertEqual(result.returncode, 0, result.stdout)
                with trace.open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(len(rows), 12)
                return rows

            def edit_preserving_metadata(path):
                before = path.stat()
                parent = path.parent.stat()
                text = path.read_text()
                # Same size and timestamps: the content hash must detect this.
                position = text.index("#") + 1
                replacement = "X" if text[position] != "X" else "Y"
                path.write_text(text[:position] + replacement + text[position + 1:])
                os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
                os.utime(path.parent, ns=(parent.st_atime_ns, parent.st_mtime_ns))

            self.assertTrue(all(r["status"] == "COMPLETED" for r in execute("initial")))
            self.assertTrue(all(r["status"] == "CACHED" for r in execute("cached", True)))
            edit_preserving_metadata(project / "Scripts/RNA-Seq/NDE_detection.R")
            changed = execute("scientific_changed", True)
            scientific = [r for r in changed if r["name"].split(" (")[0] not in ("VALIDATE_INPUTS", "VERIFY_OUTPUTS")]
            self.assertTrue(all(r["status"] == "COMPLETED" for r in scientific), changed)
            self.assertTrue(all(r["status"] == "CACHED" for r in changed if "VALIDATE_INPUTS" in r["name"]))
            edit_preserving_metadata(project / "bin/xlsx_integrity.py")
            changed = execute("helper_changed", True)
            self.assertTrue(all(r["status"] == "COMPLETED" for r in changed), changed)


if __name__ == "__main__":
    unittest.main()
