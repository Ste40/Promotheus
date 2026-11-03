#!/usr/bin/env python3
"""Exercise per-genome scheduling, failure recovery and source invalidation."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class OrthologySchedulingTests(unittest.TestCase):
    def test_parallel_pairs_resume_and_code_changes(self):
        with tempfile.TemporaryDirectory(prefix="promotheus-orthology-") as temp:
            root = Path(temp)
            project = root / "project"
            project.mkdir()
            for name in ("main.nf", "nextflow.config"):
                shutil.copy2(ROOT / name, project)
            for name in ("bin", "Scripts", "tests/fixtures"):
                shutil.copytree(ROOT / name, project / name, ignore=shutil.ignore_patterns("__pycache__"))
            base = ["nextflow", "run", str(project / "main.nf"), "-entry", "ORTHOLOGY_ONLY",
                    "-profile", "test", "--input_dir", str(project / "tests/fixtures/toy_input"),
                    "--output_dir", str(root / "output"), "--asmcode", "ASM317683v1",
                    "--threads", "4", "--memory", "2048", "--mock_delay", "1", "-ansi-log", "false"]
            def execute(name, options, success=True):
                trace = root / f"{name}.tsv"
                result = subprocess.run(base + ["-with-trace", str(trace)] + options, cwd=root,
                                        env=dict(os.environ, NXF_OFFLINE="true"), text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
                self.assertEqual(result.returncode == 0, success, result.stdout)
                with trace.open() as handle:
                    return list(csv.DictReader(handle, delimiter="\t"))
            failed = execute("failed", ["--max_orthology_parallel", "1", "--mock_fail_stage",
                                        "orthology_pair", "--mock_fail_gse", "genome_000002"], False)
            self.assertFalse(any("ORTHOLOGY_MERGE" in row["name"] for row in failed))
            resumed = execute("resumed", ["--max_orthology_parallel", "2", "-resume"])
            self.assertEqual(len(resumed), 6)
            completed_names = {row["name"] for row in failed if row["status"] == "COMPLETED"}
            self.assertTrue(completed_names)
            self.assertTrue(all(row["status"] == "CACHED" for row in resumed if row["name"] in completed_names))
            # Submission order is not guaranteed, even with maxForks=1.
            # Use a fresh successful run to measure concurrency independently
            # from how many tasks happened to finish before injected failure.
            parallel = execute("parallel", ["--max_orthology_parallel", "2"])
            intervals = []
            for row in parallel:
                if "ORTHOLOGY_PAIR" in row["name"] and row["status"] != "CACHED":
                    record = next(Path(row["workdir"]).glob("pair_*/task_record.json"))
                    intervals.append(json.loads(record.read_text()))
            self.assertEqual(len(intervals), 3)
            self.assertTrue(any(max(a["started"], b["started"]) < min(a["finished"], b["finished"])
                                for i, a in enumerate(intervals) for b in intervals[i + 1:]))
            self.assertTrue((root / "output/orthology/orthologous_genes.txt").is_file())
            cached = execute("cached", ["--max_orthology_parallel", "2", "-resume"])
            self.assertTrue(all(row["status"] == "CACHED" for row in cached))
            # Editing an existing source leaves the parent directory's metadata
            # unchanged; explicit content fingerprints must still invalidate it.
            script = project / "Scripts/Microarray/ortologhi_new.R"
            script.write_text(script.read_text() + "\n# source fingerprint regression\n")
            changed = execute("changed", ["--max_orthology_parallel", "2", "-resume"])
            self.assertTrue(all(row["status"] == "COMPLETED" for row in changed if ":ORTHOLOGY:" in row["name"]))
            self.assertTrue(all(row["status"] == "CACHED" for row in changed if "VALIDATE_INPUTS" in row["name"]))
            # A snapshot containing only the preferred genome has no pair jobs.
            # Alter only this test copy's fixture generator to exercise that DAG.
            runner = project / "bin/run_modular_stage.py"
            fixture_source = runner.read_text()
            self.assertEqual(fixture_source.count("range(1, 4)"), 1)
            runner.write_text(fixture_source.replace("range(1, 4)", "range(1, 1)"))
            reference_only = execute("reference_only", ["-resume"])
            self.assertEqual(len(reference_only), 3)
            self.assertTrue(any("ORTHOLOGY_MERGE" in row["name"] and row["status"] == "COMPLETED"
                                for row in reference_only))


if __name__ == "__main__":
    unittest.main()
