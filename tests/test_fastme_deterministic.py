#!/usr/bin/env python3
"""Distance ordering must preserve labelled values and reject ambiguous inputs."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location(
    "fastme_adapter", Path(__file__).resolve().parents[1] / "bin/fastme_deterministic.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


class DistanceOrderingTests(unittest.TestCase):
    def test_permuted_axes_and_precision(self):
        # Includes tied distances and a nonzero exponent: no value conversion.
        first = "3\nz 0.0 1e-09 2.500000\na 1e-09 0 2.500000\nm 2.500000 2.500000 0\n"
        permuted = "3\nm 0 2.500000 2.500000\nz 2.500000 0.0 1e-09\na 2.500000 1e-09 0\n"
        expected = "3\na 0 2.500000 1e-09\nm 2.500000 0 2.500000\nz 1e-09 2.500000 0.0\n"
        self.assertEqual(adapter.canonical_matrix(first), expected)
        self.assertEqual(adapter.canonical_matrix(permuted), expected)
        self.assertEqual(adapter.canonical_matrix(expected), expected)

    def test_rejects_missing_columns_duplicate_labels_and_multiple_matrices(self):
        for text in ("2\na 0 1\nb 1\n", "2\na 0 1\na 1 0\n", "1\na 0\n1\nb 0\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                adapter.canonical_matrix(text)


if __name__ == "__main__":
    unittest.main()
