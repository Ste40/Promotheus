# Orthology repeatability correction — 28 October 2026

The modular orthology adapter now uses a stable order for distance matrices
and candidate genes. It retains the pinned OrthoFinder 2.5.5/FastME 2.1.4
executables and the existing scientific thresholds. This correction applies
to `main.nf`, including `ORTHOLOGY_ONLY`; `main_legacy.nf` and the direct Docker
entrypoint retain their previous behavior.

## Diagnosis

The original repeated A12/ASM317683v1 runs had identical sequence searches,
MCL results and initial orthogroups, but different gene trees and hierarchical
orthogroups. Setting `PYTHONHASHSEED=0` did not eliminate this in the bundled
executable. Three further trials with FastME's seed fixed also varied.

Capturing FastME's input matrices showed the cause: all 37 captured matrices
(36 gene trees plus the dependency check) were identical after jointly
permuting rows and columns by taxon ID. Their values had not changed; their
input order had. The
[OrthoFinder 2.5.5 MCL reader](https://github.com/davidemms/OrthoFinder/blob/2.5.5/scripts_of/mcl.py)
returns gene sets. The
[orthologue code](https://github.com/davidemms/OrthoFinder/blob/2.5.5/scripts_of/orthologues.py)
iterates those sets to construct sequence lists and matrices, without sorting
the genes first. The resulting order is not stable across processes.

There was a second ordering issue downstream: lists of genes within `N0.tsv`
also varied, even once group membership and trees were stable. Promotheus's
greedy matching selects the first candidate among equal alignment scores.
Therefore these list permutations can change the selected paralog without
changing any scores or group membership.

## Correction

- `bin/fastme_deterministic.py` jointly orders both matrix axes lexicographically
  by taxon ID before invoking the original FastME binary. It preserves each
  distance's original decimal representation. Missing columns and duplicate
  labels fail explicitly.
- `bin/run_orthofinder.sh` provides a private runtime next to that adapter.
  OrthoFinder prepends its own `bin/` to the search path, so merely prepending
  a wrapper to the outer shell's path would not reliably intercept FastME.
  The original executables in the image are not overwritten.
- `canonical_orthogroups` in `bin/orthology_stages.R` sorts candidate gene IDs
  and group rows before calling the existing matching function. Equal-score
  decisions thus follow a defined order. It does not discard groups or change
  scores, thresholds or group membership.

Each pair retains original `N0.tsv`, the normalized `N0.canonical.tsv` used by
matching, and the ordering policy in `pair.json`. Both new executable adapters
are included in the orthology source fingerprint, so `-resume` invalidates old
orthology results after this code change. Source-content fingerprints have
subsequently been extended to every process; see the
[cache policy](MODULAR_NEXTFLOW.md#restart-and-compatibility).

This selects a reproducible resolution of ambiguous cases; it does not prove
that one tied paralog is biologically preferable to another. A corrected run
can differ from an earlier run that used an unstable order.

## Validation and scope

Three independent corrected standalone runs on the real A12/reference
proteomes produced identical 36 gene trees and identical membership of all
1976 hierarchical groups. Raw group files can still differ in gene-list
serialization order, which the matching adapter normalizes explicitly.

Two fresh Nextflow runs, each comparing A12 and AF116-22 with ASM317683v1,
completed all five tasks without cache reuse. Both pair mappings and normalized
group tables were byte-identical between runs. The published merged tables were
also byte-identical: **2565 rows × 9 columns**, SHA-256
`69e2584fecfe64c6ee774666b7ee25c6277e75be4e32fe6bd5efbf06ac738f7a`.
The complete local unit/integration suite passed after the correction. These
runs establish repeatability on this snapshot, not a new performance benchmark.

Local evidence and reproducible commands are under
`../real-runs/2026-10-28-orthology/diagnosis/`: `run_fixed.sh` (seed-only control),
`run_canonical.sh` (three corrected standalone runs), `run_pipeline.sh`
(independent Nextflow runs), and `audit_determinism.py`. The audit asserts
completed, uncached tasks and exact equality of pair mappings, normalized group
tables and published merged tables; it writes `determinism_report.json` only
after those checks pass.

Automated regressions cover joint matrix permutations, preservation of numeric
strings, malformed matrices, and permuted candidate lists with tied alignment
scores. Run them as part of `bash tests/shell/test_all.sh`.

### Full snapshot verification — 29 October 2026

Two independent runs of **all 127 comparisons across 128 archived genomes**
completed successfully against ASM317683v1. Each run completed 130 Nextflow
tasks, with **zero cached tasks**, using separate launch/cache/work directories.
Up to six comparisons ran concurrently, with one CPU / 3000 MB per pair and
an overall budget of eight CPUs / 20,000 MB.

The automatic audit passed:

- All 127 pair mappings and normalized group tables were byte-identical.
- Group membership matched for every comparison.
- The final orthology tables were byte-identical: **2565 rows × 339 columns**,
  with **zero differing cells**.
- All 179 archived input/source files retained their recorded SHA-256 checksums.

Final table SHA-256:
`3f1112596754007d2d94935ae843048f64bfad68ef5494e29cd8ae6705c579fc`.

Evidence is archived under `../real-runs/2026-10-29-orthology-full/`, including
`full_validation.json`, `audit.py`, `run_all.sh`, `parameters.json`,
`snapshot_sha256.json`, the frozen project/input/genome snapshot, both traces,
and both published tables. The audit completed at 15:47 CET on 29 October 2026.
This extends the initial two-comparison check to the entire archived genome set;
it verifies the orthology branch, not a repeated downstream transcriptomics run.

Validation is limited to the archived genomes and pinned image. It does not
establish repeatability across every genome, changed input snapshots, tool
versions or numerical platforms, nor full biological equivalence with legacy
results.
