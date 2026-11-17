# Two GSEs per species: Nextflow validation

The small real-data test completed on 17 November 2026 with `--multitaxa true --operons_detection false`, using scientific source `b319fc0`. Test-only hooks staged existing inputs and recorded the actual ComBat calls without changing their calculations.

## Inputs and scope

| Species | GSE | Technology | Samples |
| --- | --- | --- | ---: |
| Bacillus subtilis | GSE54577 | Microarray | 6 |
| Bacillus subtilis | GSE64918 | Two-channel microarray | 6 |
| Lactococcus lactis | GSE116952 | Paired-end RNA-Seq | 4 |
| Lactococcus lactis | GSE74554 | Single-end RNA-Seq | 6 |

RNA inputs contain the first 100,000 reads per FASTQ, preserving paired-mate order. References ASM317683v1 and ASM904v1 were staged from an existing snapshot; ASM317683v1 was the preferred reference. GSE54577 reused its previously validated real full-table checkpoint. GSE64918 was analyzed from its six original raw files.

GEOquery initially reported no supplementary files for GSE64918 despite an available public RAW archive. The first attempt failed explicitly rather than silently dropping the GSE. The archive was downloaded directly from NCBI and decompressed into the test seed, then Nextflow was resumed. This is a documented test workaround, not a production downloader fix. It does not validate automatic recovery from this GEO listing failure.

## Results

- Final output: **246 genes, 22 expression columns, 52 total columns**.
- Exactly two actual ComBat calls: 12 Bacillus samples (6 + 6), and 10 Lactococcus samples (4 + 6). Each call included only the two GSEs belonging to that species.
- Batch labels matched the corresponding input column GSE identifiers; all corrected values were finite. Maximum absolute expression changes were approximately 1.2612 and 3.4819, confirming correction was applied.
- Final expression values and CVs were finite, final gene means were positive, and CV matched sample SD / mean to tolerance `1e-10`.
- The production output validator passed. All 146 recorded project, input and seed file hashes matched the frozen manifest.

These checks establish execution and GSE-level batch assignment on this fixture. They do not establish biological accuracy, individual sample identity independently of the input metadata, or validity of cross-technology expression comparisons. The two species use different technologies, so species and technology are confounded. CV is computed on the integrated expression scale. Operon detection was disabled for this test and is covered separately in [the operon validation report](VALIDATION_OPERONS_CV_DISK.md).

## Task timings

Durations are wall-clock task times, not CPU times. The successful resumed attempt ran new tasks over **227.5 seconds (3 min 48 s)**, measured from its first newly started task to the last completion. Eight tasks were recovered from the previous attempt and six ran anew. This is not a cold-run or full-download benchmark; do not sum overlapping durations or count cached timings as new work.

| Task | Duration (seconds) | Final attempt |
| --- | ---: | --- |
| Validate inputs | 0.473 | Cached |
| Prepare orthology | 21.571 | Cached |
| Prepare RNA inputs: GSE116952 | 31.144 | Cached |
| Prepare RNA inputs: GSE74554 | 36.155 | Cached |
| Orthology pair | 225.407 | Cached |
| Merge orthology | 14.646 | Cached |
| Align GSE116952 | 150.948 | Cached |
| Align GSE74554 | 182.093 | Executed |
| Analyze GSE116952 | 18.609 | Cached |
| Analyze GSE74554 | 18.913 | Executed |
| Microarray GSE54577 | 25.815 | Executed; checkpoint restoration |
| Microarray GSE64918 | 32.533 | Executed; raw data analysis |
| Integrate results | 26.050 | Executed |
| Verify outputs | 0.146 | Executed |

The runner used up to 6 configured CPUs and 16,000 MB globally, 2 CPUs / 6,000 MB per analysis task, and 1 CPU / 3,000 MB for the single orthology pair. Docker ran the pinned amd64 image under emulation on the ARM host. Container CPU/RSS metrics were unavailable through the nested Docker launcher and are not treated as zero. Bulk reference/SRA downloads and the separately downloaded GSE64918 archive are excluded; live metadata requests are included.

## Local evidence

The preserved directory `real-runs/2026-11-17-multitaxa-four-gse/` outside the repository contains `full_validation.json`, `combat_audit.rds`, `source_integrity.json`, `snapshot_sha256.json`, `TASK_TIMINGS.md`, Nextflow traces, inputs, final outputs and the failed-attempt evidence in `attempt_1/`.

Repeat the independent checks on these retained artifacts with:

```bash
Rscript real-runs/2026-11-17-multitaxa-four-gse/audit.R real-runs/2026-11-17-multitaxa-four-gse
python3 real-runs/2026-11-17-multitaxa-four-gse/task_timings.py
```

## Follow-up: GEO downloader correction

The production microarray downloader now falls back to the canonical NCBI HTTPS RAW archive when GEOquery fails or does not produce a valid archive. Downloads use a temporary file and are promoted only after archive validation. Exhausted retries raise an explicit error identifying the GSE and URL; a success on the third attempt is now accepted. Extraction failures also stop the task.

`Rscript tests/test_geo_raw_download.R .` covers empty/failed listings, the normal GEOquery path, existing valid archives, invalid responses, third-attempt success, nonzero download status, and partial-file cleanup. The test is included in the regression suite.

In the pinned production container, a fresh live GSE64918 download passed. A separate test simulated only an empty listing and performed the real HTTPS fallback: its archive matched the GEOquery-downloaded archive byte-for-byte by digest and contained the six original GSM files. Evidence is retained locally in `real-runs/2026-11-17-geo-download-fix/`. This verifies the download component; the four-GSE scientific run above used the earlier manual staging workaround.

The subsequent [clean-download, interruption and reproducibility validation](VALIDATION_CLEAN_RESUME_REPRO.md) exercised both microarrays without checkpoints, found and corrected a separate GSM metadata issue, and compared two complete results.
