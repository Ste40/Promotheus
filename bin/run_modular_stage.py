#!/usr/bin/env python3
"""Run one isolated scientific task; all dependencies are explicit Nextflow paths."""
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from prepare_gse_inputs import write_xlsx
from validate_promotheus_inputs import first_sheet_rows


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.digest()


def copy_tree(source, target, verify_existing=False):
    source, target = Path(source), Path(target)
    if not source.exists():
        return
    for item in sorted(source.rglob("*")):
        dest = target / item.relative_to(source)
        if item.is_dir():
            dest.mkdir(parents=True, exist_ok=True)
        elif item.is_file():
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists() and verify_existing:
                if digest_file(item) != digest_file(dest):
                    raise ValueError(f"Conflicting shared reference artifact: {dest}")
            else:
                shutil.copy2(item, dest)


def descriptors(dependencies):
    result = []
    for dependency in dependencies:
        if (dependency / "descriptor.json").exists():
            result.append((json.loads((dependency / "descriptor.json").read_text()), dependency))
    return sorted(result, key=lambda item: (item[0]["ordinal"], item[0]["accession"]))


def assemble(stage, dependencies, output):
    if stage in ("orthology_prepare", "orthology_merge") and dependencies:
        reference = dependencies[0]
        summary = reference / "assembly_summary.txt"
        genomes = reference / "tmp/orthologous_genome"
        if not summary.is_file() or not list(genomes.glob("*.gb")):
            raise ValueError("Orthology snapshot requires assembly_summary.txt and tmp/orthologous_genome/*.gb")
        shutil.copy2(summary, output)
        destination = output / "tmp/orthologous_genome"
        destination.mkdir(parents=True, exist_ok=True)
        for genome in genomes.glob("*.gb"):
            shutil.copy2(genome, destination)
        if stage == "orthology_merge":
            copy_tree(reference / "shared", output / "shared")
            shutil.copy2(reference / "pairs.tsv", output)
            mappings = output / "mappings"
            mappings.mkdir()
            for pair in dependencies[1:]:
                info = json.loads((pair / "pair.json").read_text())
                key = info["key"]
                if not key.startswith("genome_") or not key[7:].isdigit():
                    raise ValueError("Invalid orthology pair key")
                target = mappings / f"{key}.rds"
                if target.exists():
                    raise ValueError(f"Duplicate orthology pair: {key}")
                shutil.copy2(pair / "mapping.rds", target)
        return
    if stage == "orthology_pair":
        copy_tree(dependencies[0], output / "shared")
        return
    for dependency in dependencies:
        if (dependency / "orthologous_genes.txt").exists():
            shutil.copy2(dependency / "orthologous_genes.txt", output)
            if (dependency / "assembly_summary.txt").exists():
                shutil.copy2(dependency / "assembly_summary.txt", output)
            copy_tree(dependency / "tmp/orthologous_genome", output / "tmp/orthologous_genome")
    if stage == "rna_align":
        copy_tree(dependencies[0], output)
    elif stage == "rna_analyze":
        states = output / "states"
        states.mkdir()
        for i, (descriptor, dependency) in enumerate(descriptors(dependencies)):
            # Only merge the organism tree, never the per-GSE state or group workbook.
            copy_tree(dependency / "RNA_Seq" / descriptor["organism"],
                      output / "RNA_Seq" / descriptor["organism"], verify_existing=True)
            if (dependency / "state.rds").exists():
                shutil.copy2(dependency / "state.rds", states / f"{i:08d}.rds")
    elif stage == "integrate":
        for dependency in dependencies:
            copy_tree(dependency / "results/fulltab", output / "tmp/microarray")
            copy_tree(dependency / "RNA_Seq/results", output / "RNA_Seq/results")
            # The R operon filter discovers organisms by directory name.
            # Carry the names forward, without staging reads or BAM files.
            for organism in (dependency / "RNA_Seq").glob("*"):
                if organism.is_dir() and organism.name != "results":
                    (output / "RNA_Seq" / organism.name).mkdir(parents=True, exist_ok=True)


def mock(config, inputs, dependencies, output):
    stage = config["stage"]
    time.sleep(float(config.get("mock_delay", 0)))
    if config.get("mock_fail_stage") == stage and (
            not config.get("mock_fail_gse") or config["key"] == config["mock_fail_gse"]):
        raise RuntimeError("Deliberate mock task failure")
    if stage == "orthology_prepare":
        (output / "shared").mkdir()
        genomes = output / "tmp/orthologous_genome"
        genomes.mkdir(parents=True, exist_ok=True)
        (genomes / "reference.gb").write_text("mock reference\n")
        (output / "assembly_summary.txt").write_text("mock summary\n")
        rows = ["ordinal\tkey"]
        for i in range(1, 4):
            key = f"genome_{i:06d}"
            pair = output / "pairs" / key
            pair.mkdir(parents=True)
            (pair / "pair.json").write_text(json.dumps({"key": key, "ordinal": i}))
            rows.append(f"{i}\t{key}")
        (output / "pairs.tsv").write_text("\n".join(rows) + "\n")
    elif stage == "orthology_pair":
        shutil.copy2(inputs / "pair.json", output)
        (output / "mapping.rds").write_text(config["key"])
    elif stage == "orthology_merge":
        with (output / "pairs.tsv").open() as handle:
            expected = [row["key"] for row in csv.DictReader(handle, delimiter="\t")]
        actual = sorted(p.read_text() for p in (output / "mappings").glob("*.rds"))
        if sorted(expected) != actual:
            raise ValueError("Missing orthology pair results")
        (output / "orthologous_genes.txt").write_text("geneA\tTOY_0001\ngeneB\tTOY_0002\n")
    elif stage == "orthology":
        (output / "orthologous_genes.txt").write_text("geneA\tTOY_0001\ngeneB\tTOY_0002\n")
    elif stage in ("rna_download", "rna_align"):
        (output / f"{stage}.complete").write_text(config["key"] + "\n")
    elif stage in ("microarray", "rna_analyze"):
        records = [json.loads((inputs / "descriptor.json").read_text())] if stage == "microarray" else [d for d, _ in descriptors(dependencies)]
        folder = output / ("results/fulltab" if stage == "microarray" else "RNA_Seq/results")
        folder.mkdir(parents=True, exist_ok=True)
        for record in records:
            (folder / f'{record["accession"]}_full_tab.txt').write_text(
                "GeneID\tlocus_tag\tmeans\tsd\tdiff_signif\n"
                "geneA\tTOY_0001\t10\t1\tNO\n")
    elif stage == "integrate":
        accessions = sorted(p.name.split("_")[0] for folder in (output / "tmp/microarray", output / "RNA_Seq/results") for p in folder.glob("*full_tab.txt"))
        expected = []
        for workbook in inputs.glob("*.xlsx"):
            rows = first_sheet_rows(workbook)
            expected += [r[rows[0].index("Accession")] for r in rows[1:] if any(r)]
        if sorted(expected) != accessions:
            raise ValueError("Integration received missing or duplicate GSE results")
        results = output / "results"
        results.mkdir(exist_ok=True)
        # Original two-GSE fixture remains byte-identical for legacy regression checks.
        fixture = Path(config["fixture"])
        if accessions == ["GSE000001", "GSE000002"]:
            shutil.copy2(fixture, results / "final_table.xlsx")
        else:
            fixture_rows = first_sheet_rows(fixture)
            headers = ["Gene_name", "GeneID"] + ["Expression_" + a for a in accessions] + ["Mean_Expression", "CV"] + ["Percentile_" + a for a in accessions] + ["Mean_Percentile", "Percentile_SD", "DE", "Sequence"]
            rows = [headers]
            for i, row in enumerate(fixture_rows[1:], 1):
                rows.append([row[0], row[1]] + [i * 10] * len(accessions) + [i * 10, 0] + [i / 4] * len(accessions) + [i / 4, 0, "NO", row[-1]])
            write_xlsx(results / "final_table.xlsx", rows)
        copy_tree(output / "tmp/microarray", results / "fulltab")


def docker(config, inputs, output):
    args = ["docker", "run", "--rm", "-i"]
    platform = config["container_platform"]
    if platform != "auto":
        args += ["--platform", platform]
    # Resource limits are enforced both by the local executor and the container.
    args += ["--cpus", str(config["cpus"]), "--memory", f'{config["memory_mb"]}m']
    variables = {"THREADS": config["cpus"], "MEMORY": config["memory_mb"],
        "ASMcode": config["asmcode"], "CORRECTION_DE": config["correction_de"],
        "PVALUE": config["pvalue"], "OPERONS_DETECTION": config["operons_detection"],
        "OPERONS_THRESHOLD": config["operons_threshold"], "AMINO_DIFF": config["amino_diff"],
        "SIMILARITY": config["similarity"], "MULTITAXA": config["multitaxa"],
        "MAX_GSE_PARALLEL": 1, "PROMOTHEUS_KEEP_CACHE": "TRUE",
        "OMP_NUM_THREADS": config["cpus"], "OPENBLAS_NUM_THREADS": 1}
    for name, value in variables.items():
        args += ["-e", f'{name}={str(value).upper() if isinstance(value, bool) else value}']
    args += ["-v", f"{inputs.resolve()}:/input:ro", "-v", f"{output.resolve()}:/output",
             "-v", f'{Path(config["scripts"]).resolve()}:/scripts:ro',
             "-v", f"{Path(__file__).resolve().parent}:/adapters:ro",
             "--entrypoint", "Rscript", config["container_image"], "/adapters/run_stage.R", config["stage"]]
    subprocess.run(args, input="\n\n", text=True, check=True)


def run(config):
    output = Path(config["output"])
    output.mkdir()
    inputs = Path("task_input")
    inputs.mkdir()
    for item in config["inputs"]:
        item = Path(item)
        if item.is_dir():
            copy_tree(item, inputs)
        else:
            shutil.copy2(item, inputs / item.name)
    dependencies = [Path(p) for p in config["dependencies"]]
    assemble(config["stage"], dependencies, output)
    descriptor = inputs / "descriptor.json"
    if descriptor.exists():
        shutil.copy2(descriptor, output)
    started = time.time()
    if config["engine"] == "mock":
        mock(config, inputs, dependencies, output)
    else:
        docker(config, inputs, output)
    stage = config["stage"]
    if stage == "integrate":
        helpers = Path(__file__).resolve().parent
        workbook = output / "results/final_table.xlsx"
        subprocess.run([sys.executable, str(helpers / "repair_promotheus_xlsx.py"), "--workbook", str(workbook)], check=True)
        if config["export_tsv"]:
            subprocess.run([sys.executable, str(helpers / "export_promotheus_tsv.py"), "--workbook", str(workbook), "--output", str(output / "results/final_table.tsv")], check=True)
        metadata = output / ".promotheus-nextflow"
        metadata.mkdir(exist_ok=True)
        shutil.copy2(config["manifest"], metadata / "input_manifest.tsv")
        # Provenance must reference only published artifacts, not private integration copies.
        # Reusable orthology/MA data remain in their upstream Nextflow task directories.
        shutil.rmtree(output / "tmp", ignore_errors=True)
        subprocess.run([sys.executable, str(helpers / "build_checkpoint_manifest.py"), "--output-dir", str(output), "--input-manifest", str(metadata / "input_manifest.tsv"), "--manifest", str(metadata / "checkpoint_manifest.tsv")], check=True)
        with (metadata / "run_record.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["field", "value"])
            for key, value in config.items():
                if key not in ("inputs", "dependencies"):
                    writer.writerow([key, str(value).lower() if isinstance(value, bool) else value])
    # Do not forward downloaded reads/BAMs to integration; keep them in their task cache.
    if stage in ("microarray", "rna_analyze"):
        keep = "results" if stage == "microarray" else "RNA_Seq/results"
        retained = output / keep
        if not retained.is_dir() or not list(retained.rglob("*full_tab.txt")):
            raise ValueError("Missing GSE results")
        if stage == "rna_analyze":
            for child in (output / "RNA_Seq").iterdir():
                if child.name != "results":
                    if child.is_dir():
                        shutil.rmtree(child)
                        # Preserve the organism identity used by integration.
                        child.mkdir()
                    else:
                        child.unlink()
        else:
            shutil.rmtree(output / "tmp", ignore_errors=True)
            (output / "orthologous_genes.txt").unlink(missing_ok=True)
    (output / "task_record.json").write_text(json.dumps({"stage": stage, "key": config["key"], "started": started, "finished": time.time(), "pid": os.getpid(), "cpus": config["cpus"]}, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(json.loads(Path(args.config).read_text()))
