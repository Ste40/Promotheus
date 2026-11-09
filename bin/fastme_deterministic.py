#!/usr/bin/env python3
"""Stabilize OrthoFinder 2.5.5's unordered distance matrices before FastME."""
import os
from pathlib import Path
import sys


def canonical_matrix(text):
    lines = text.splitlines()
    size = int(lines[0].strip())
    rows = [line.split() for line in lines[1:] if line.strip()]
    if size < 1 or len(rows) != size or any(len(row) != size + 1 for row in rows):
        raise ValueError("Expected one full square PHYLIP distance matrix")
    labels = [row[0] for row in rows]
    if len(set(labels)) != size:
        raise ValueError("Duplicate taxon labels in distance matrix")
    # Permute BOTH axes; preserve the original decimal strings without rounding.
    order = sorted(range(size), key=lambda i: labels[i])
    return str(size) + "\n" + "".join(
        labels[i] + " " + " ".join(rows[i][j + 1] for j in order) + "\n"
        for i in order)


def main():
    args = sys.argv[1:]
    if "-i" in args:
        matrix = Path(args[args.index("-i") + 1])
        matrix.write_text(canonical_matrix(matrix.read_text()))
    # The image's original executable and all inference options stay unchanged.
    os.execv("/OrthoFinder/bin/fastme", ["fastme"] + args)


if __name__ == "__main__":
    main()
