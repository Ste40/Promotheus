#!/usr/bin/env python3
"""Offline integration: real Nextflow scheduling, failure propagation and resume."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from prepare_gse_inputs import prepare, write_xlsx
from validate_promotheus_inputs import first_sheet_rows
from xlsx_integrity import validate_xlsx_package
from run_modular_stage import assemble, digest_file, run


def make_inputs(folder):
    folder.mkdir(parents=True, exist_ok=True)
    for kind, accessions in (("Microarrays", ["GSE000001", "GSE000004"]),
                              ("RNA_Seq", ["GSE000002", "GSE000003", "GSE000005"])):
        source = ROOT / f"tests/fixtures/toy_input/GEO_Datasets_Search_{kind}_Toy.xlsx"
        original = first_sheet_rows(source)
        header = original[0]
        rows = [header]
        for accession in accessions:
            row = list(original[1])
            row[header.index("Accession")] = accession
            row[header.index("GSE")] = accession[3:]
            rows.append(row)
        write_xlsx(folder / source.name, rows)


class InputAndIsolationTests(unittest.TestCase):
    def test_organism_names_reach_integration_without_raw_data(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            dependency = directory / "aligned"
            organism = dependency / "RNA_Seq/Toy"
            organism.mkdir(parents=True)
            (organism / "raw.bam").write_text("large intermediate")
            (dependency / "descriptor.json").write_text(json.dumps({
                "organism": "Toy", "ordinal": 0, "accession": "GSE000002"}))
            (dependency / "state.rds").write_text("state")
            previous = Path.cwd()
            try:
                os.chdir(directory)
                run({"stage": "rna_analyze", "key": "Toy", "output": "analysis",
                     "inputs": [], "dependencies": [str(dependency)],
                     "engine": "mock", "cpus": 1})
            finally:
                os.chdir(previous)
            retained = directory / "analysis/RNA_Seq/Toy"
            self.assertTrue(retained.is_dir(), "Operon integration discovers organisms from directory names")
            self.assertEqual(list(retained.iterdir()), [])
            output = directory / "integrated"
            output.mkdir()
            assemble("integrate", [directory / "analysis"], output)
            self.assertTrue((output / "RNA_Seq/Toy").is_dir())
            self.assertFalse(list(output.rglob("*.bam")))
            self.assertTrue(list((output / "RNA_Seq/results").glob("*full_tab.txt")))

    def test_shards_and_ambiguous_references(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            make_inputs(directory / "input")
            prepare(directory / "input", directory / "shards")
            shards = list((directory / "shards").glob("*/*.xlsx"))
            self.assertEqual(len(shards), 5)
            for shard in shards:
                self.assertEqual(validate_xlsx_package(shard), [])
                self.assertEqual(len(first_sheet_rows(shard)), 2)
            rna = next((directory / "input").glob("*RNA_Seq*"))
            rows = first_sheet_rows(rna)
            rows[2][rows[0].index("reference_genome")] = "ASM999999v1"
            write_xlsx(rna, rows)
            with self.assertRaisesRegex(ValueError, "conflicting"):
                prepare(directory / "input", directory / "invalid")

    def test_shared_genomes_copied_once_and_conflicts_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            dependencies = []
            for i in range(2):
                dependency = directory / str(i)
                organism = dependency / "RNA_Seq/Toy"
                organism.mkdir(parents=True)
                (organism / "genome.fna").write_text("ACGT")
                (dependency / "descriptor.json").write_text(json.dumps({"organism": "Toy", "ordinal": i, "accession": f"GSE{i}"}))
                (dependency / "state.rds").write_text(str(i))
                dependencies.append(dependency)
            output = directory / "output"
            output.mkdir()
            assemble("rna_analyze", list(reversed(dependencies)), output)
            self.assertEqual((output / "states/00000000.rds").read_text(), "0")
            self.assertEqual((output / "states/00000001.rds").read_text(), "1")
            self.assertEqual((output / "RNA_Seq/Toy/genome.fna").read_text(), "ACGT")
            (dependencies[1] / "RNA_Seq/Toy/genome.fna").write_text("DIFFERENT")
            (directory / "conflict").mkdir()
            with self.assertRaisesRegex(ValueError, "Conflicting shared reference"):
                assemble("rna_analyze", dependencies, directory / "conflict")


class NextflowTests(unittest.TestCase):
    def test_parallel_failure_resume_and_operon_barrier(self):
        configured = os.environ.get("PROMOTHEUS_MODULAR_TEST_DIR")
        temporary = None if configured else tempfile.TemporaryDirectory(prefix="promotheus-modular-")
        directory = Path(configured or temporary.name).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        make_inputs(directory / "input")
        base = ["nextflow", "run", str(ROOT / "main.nf"), "-profile", "test",
                "--input_dir", str(directory / "input"), "--output_dir", str(directory / "output"),
                "--asmcode", "ASM317683v1", "--threads", "4", "--memory", "2048",
                "--max_gse_parallel", "2", "--mock_delay", "0.7", "-ansi-log", "false"]
        environment = dict(os.environ, NXF_OFFLINE="true")
        def execute(name, extra, expect_success=True):
            # The per-run explicit trace also avoids timestamp collisions between fast launches.
            command = base + ["-with-trace", str(directory / f"{name}.tsv")] + extra
            result = subprocess.run(command, cwd=directory, env=environment, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
            (directory / f"{name}.log").write_text(result.stdout)
            self.assertEqual(result.returncode == 0, expect_success, result.stdout)
            with (directory / f"{name}.tsv").open() as handle:
                return list(csv.DictReader(handle, delimiter="\t"))
        failure = execute("failure", ["--operons_detection", "false", "--mock_fail_stage", "rna_align",
                                       "--mock_fail_gse", "GSE000003"], False)
        self.assertTrue(any(row["status"] == "FAILED" and "GSE000003" in row["name"] for row in failure))
        self.assertFalse(any("INTEGRATE_RESULTS" in row["name"] for row in failure))
        completed = execute("resumed", ["--operons_detection", "false", "-resume"])
        self.assertEqual(len(completed), 19)
        self.assertTrue(any(row["status"] == "CACHED" for row in completed))
        self.assertEqual(sum("RNA_ANALYZE" in row["name"] for row in completed), 3)
        report = json.loads((directory / "output/.promotheus-nextflow/output_validation.json").read_text())
        self.assertEqual(report["status"], "pass")
        self.assertEqual(len(report["expected_accessions"]), 5)
        self.assertEqual(report["missing_accessions"], [])
        for artifact in ("results/final_table.xlsx", "results/final_table.tsv", "orthologous_genes.txt",
                         ".promotheus-nextflow/checkpoint_manifest.tsv", ".promotheus-nextflow/run_record.tsv",
                         ".promotheus-nextflow/input_validation.json"):
            self.assertTrue((directory / "output" / artifact).is_file(), artifact)
        with (directory / "output/.promotheus-nextflow/checkpoint_manifest.tsv").open() as handle:
            for checkpoint in csv.DictReader(handle, delimiter="\t"):
                for artifact, checksum in zip(checkpoint["artifacts"].split(","), checkpoint["sha256"].split(",")):
                    published = directory / "output" / artifact
                    self.assertTrue(published.is_file(), artifact)
                    self.assertEqual(digest_file(published).hex(), checksum)
        # Inspect actual task execution intervals, including cached successful tasks from run one.
        intervals = []
        for row in completed:
            if "RNA_DOWNLOAD_GSE" not in row["name"]:
                continue
            work = Path(row["workdir"])
            for record in work.glob("*/task_record.json"):
                data = json.loads(record.read_text())
                if data["stage"] == "rna_download":
                    intervals.append(data)
        self.assertEqual(len(intervals), 3)
        # A second resume must reuse every task without executing scientific work again.
        cached = execute("cached", ["--operons_detection", "false", "-resume"])
        self.assertTrue(all(row["status"] == "CACHED" for row in cached), cached)
        operons = execute("operons", ["--operons_detection", "true", "-resume"])
        self.assertEqual(sum("RNA_ANALYZE" in row["name"] for row in operons), 1)
        self.assertTrue(all(row["status"] == "CACHED" for row in operons if "RNA_DOWNLOAD_GSE" in row["name"] or "RNA_ALIGN_GSE" in row["name"]))
        all_records = [json.loads(p.read_text()) for p in (directory / "work").glob("*/*/*/task_record.json")]
        tasks = [r for r in all_records if r["stage"] in ("orthology_pair", "rna_download", "rna_align", "microarray")]
        self.assertTrue(any(a["pid"] != b["pid"] and max(a["started"], b["started"]) < min(a["finished"], b["finished"])
                            for i, a in enumerate(tasks) for b in tasks[i+1:]), "No concurrent task execution observed")
        print("PASS: five GSEs; concurrent task intervals; failure isolation; resume; operon grouping")
        if temporary:
            temporary.cleanup()


if __name__ == "__main__":
    unittest.main()
