# Nextflow performance measurements

See the [per-task report](TASK_PERFORMANCE.md) for durations, concurrency,
resource requests and all 136 tasks of the complete run.

## Full workflow: 3 November 2026

The `standard,balanced` profile ran the modular workflow at commit `5cf286b`
with 8 CPUs and a 20,000 MB local budget. Docker reported 8 CPUs and about
22 GiB RAM. Scientific tasks requested 2 CPUs / 8000 MB; orthology pairs
requested 1 CPU / 3000 MB. The requested six pairs were capped at four to
reserve resources for RNA-Seq.

Inputs: GSE8182 (six Microarray samples), GSE116952 (four RNA-Seq samples),
reference ASM317683v1, 128 archived genomes, operons disabled. GEO/SRA data
were downloaded afresh; genome download time was excluded.

| Measurement | Result |
| --- | --- |
| Completed tasks | 136; none cached |
| Genome comparisons | 127 |
| Workflow task span | 3 h 42 min 30 s |
| Orthology comparison span | 3 h 39 min 7 s |
| RNA-Seq download | 33 min 17 s |
| Download-to-alignment waiting time | 1.135 s |
| RNA-Seq alignment | 54 min 26 s |
| RNA-Seq analysis | 4 min 34 s |
| Final table | 352 genes × 28 columns; all ten samples |

Stages overlap: their times must not be added. All 179 frozen source/input
files passed checksum verification. The final orthology map matched the
independently validated snapshot (SHA-256
`3f1112596754007d2d94935ae843048f64bfad68ef5494e29cd8ae6705c579fc`).
The final table had the same genes, columns and values as the preceding full
run: text compared exactly, numeric tolerance `rtol=1e-8`, `atol=1e-10`.

The preceding resumed run waited 2 h 13 min 10 s before alignment and its
newly executed tasks spanned 3 h 16 min 39 s. It reused 23 tasks, including
21 orthology pairs. These are **not comparable cold-run totals**. The shorter
queue wait demonstrates improved overlap; the measurements do not establish
an overall workflow speedup. Alignment itself also took longer under the new
concurrent workload (54 min versus 27 min previously).

## Controlled comparison: 4 November 2026

Fresh Nextflow runs used reference ASM317683v1 and archived A12 / AF116-22
genomes, with the same scientific code as `5cf286b`. Each run completed five
tasks (validation, preparation, two comparisons and merge), without resume.

| Measurement | Serial | Parallel |
| --- | --- | --- |
| Concurrent comparisons observed | 1 | 2 |
| Workflow task span | 525.784 s (8 min 46 s) | 297.766 s (4 min 58 s) |
| Pair comparison span | 495.903 s | 271.376 s |
| Whole command, including Nextflow startup/shutdown | 534.119 s | 305.799 s |

The workflow task span improved **1.77×**, a **43.4% reduction**. Both individual
mapping files and the final orthology table were byte-identical. Source/input
checksums were unchanged throughout. This demonstrates parallelism on these
two real comparisons, not an end-to-end RNA-Seq/Microarray speedup.

The [machine-readable result](benchmarks/2026-11-04-orthology.json) includes
hashes, timings, budgets and limitations. This is a single serial-first pair
of runs, not a repeated statistical benchmark.

## Repeatable orthology comparison

Use the same local input and genome snapshot for two fresh runs. The harness
runs `ORTHOLOGY_ONLY` first with one pair at a time, then with two concurrent
pairs, keeping the total budget (4 CPUs / 12000 MB) and per-task requests
(1 CPU / 3000 MB) fixed:

```bash
python3 tests/benchmark_orthology.py \
  --input /absolute/path/to/input \
  --genomes /absolute/path/to/archived-genomes \
  --output /absolute/path/to/new-benchmark-directory \
  --asmcode ASM317683v1
```

The output directory must not exist. The genome snapshot must contain the
preferred reference and at least two other genomes in the layout accepted by
`--orthology_genomes`. The harness checks completed task coverage, actual
concurrency, source/input checksums, byte-identical per-pair mappings and the
final orthology table. It writes `benchmark.json`, `status.json`, commands,
logs and the Nextflow reports for both runs. Failures return a nonzero exit
code and a failed status.

This is an orthology-only test. Neither run uses `-resume`, but OS/Docker
filesystem caches are not flushed. One fixed-order serial/parallel comparison
cannot quantify run-to-run variance or remove host-load/order effects. Nextflow
trace CPU/RSS fields do not reliably describe the nested Docker processes on
macOS; no CPU or peak-memory claim is made from those fields.
