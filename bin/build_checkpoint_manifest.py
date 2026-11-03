#!/usr/bin/env python3
"""Record which requested Promotheus GSE accessions have reusable full-tab outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_accessions(path: Path) -> list[tuple[str, str]]:
    accessions: list[tuple[str, str]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            kind = row.get("kind", "unknown")
            accessions.extend(
                (kind, accession.strip())
                for accession in row.get("accessions", "").split(",")
                if accession.strip()
            )
    return accessions


def full_tabs(output_dir: Path) -> dict[str, list[Path]]:
    by_accession: dict[str, list[Path]] = {}
    candidates = list((output_dir / "results" / "fulltab").glob("*full_tab.txt"))
    candidates += list((output_dir / "RNA_Seq" / "results").glob("*full_tab.txt"))
    candidates += list((output_dir / "tmp" / "microarray").glob("*full_tab.txt"))
    for path in candidates:
        match = re.search(r"GSE\d+", path.name)
        if match and path.is_file() and path.stat().st_size > 0:
            by_accession.setdefault(match.group(0), []).append(path)
    return by_accession


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()

    available = full_tabs(args.output_dir)
    with args.manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("kind", "accession", "status", "artifacts", "sha256"))
        for kind, accession in expected_accessions(args.input_manifest):
            paths = sorted(available.get(accession, []))
            writer.writerow(
                (
                    kind,
                    accession,
                    "complete" if paths else "missing",
                    ",".join(str(path.relative_to(args.output_dir)) for path in paths),
                    ",".join(sha256(path) for path in paths),
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
