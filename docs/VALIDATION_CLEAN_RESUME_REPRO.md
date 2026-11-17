# Clean GEO downloads, interruption recovery and reproducibility

## Result

Two real four-GSE Nextflow executions passed on 17 November 2026. Both produced **246 genes × 52 columns**, including **22 expression columns**. After matching GeneID and column names, every numeric value matched exactly (maximum absolute difference **0**, acceptance tolerance `1e-10`); nonnumeric annotations were identical. Both runs passed the production output validator and the independent ComBat/CV audit.

Inputs were the same as [the four-GSE fixture](VALIDATION_MULTITAXA_FOUR_GSE.md): GSE54577 and GSE64918 (6 microarray samples each, Bacillus subtilis), GSE116952 and GSE74554 (4 paired-end and 6 single-end RNA samples, Lactococcus lactis). Both microarrays were downloaded and analyzed afresh in each run, without raw-file seeds or full-table checkpoints. RNA used the first 100,000 reads per FASTQ and archived reference genomes. Thus this is a clean Nextflow/GEO-microarray test, **not a full raw RNA download benchmark**. Operons were disabled; `multitaxa` was enabled.

The source was commit `9c158e8` plus the family-GSM metadata correction described below. Test-only hooks staged subsampled RNA inputs, labeled Docker containers, and recorded actual ComBat calls without changing the calculations. All 131 frozen project/input/genome/seed hashes passed, and the production metadata source matched the tested snapshot.

## Failures discovered and correction

An early attempt failed to parse a GEO RNA series matrix. An unchanged retry succeeded; automatic RNA metadata retry was not added by this change. This remains a transient-failure limitation.

A later attempt successfully recovered both microarray RAW archives through the new HTTPS fallback, but failed while retrieving individual GSM metadata. Inspection showed an HTML GEO error page cached under a `.soft` filename. The complete GSE family response already contained the same sample headers.

The correction saves those headers and uses them, in requested GSM order, for experimental grouping. Missing samples fail explicitly. Legacy callers without the family snapshot retain individual GSM retrieval, now removing invalid cached SOFT responses between attempts and accepting success on attempt three. Focused tests and the complete regression suite passed. The earlier sources, logs and outputs remain separately preserved; neither failed attempt is represented as a successful benchmark.

Both final runs completed microarray analysis without individual-GSM metadata requests. In the second run GSE64918 again exercised the real automatic RAW fallback. A nonfatal legacy cleanup warning about trying to remove a directory when there are no individual SOFT files was observed; it did not affect output or remove the data directory.

## Interruption and partial downloads

Run A received SIGINT while tasks were active, after RNA input preparation had completed. Only containers labeled for this test were stopped. Its subsequent `-resume` recovered four tasks (input validation, orthology preparation and both RNA input tasks); all previously completed scientific tasks were reused. Incomplete tasks ran again and the full workflow passed.

A separate component check injected a connection error after writing a partial RAW response, then used the real HTTPS downloader on retry. The second attempt produced a valid archive and removed its temporary file. This tests ordinary transfer-error recovery; it does not claim cleanup after SIGKILL or a host crash during transfer.

Run B used a separate work directory and cache, with all 14 tasks executed anew. Genome and subsampled RNA inputs were intentionally identical to Run A.

## Timing and space

Run A's resumed invocation took **390.94 s (6 min 31 s)**; this excludes its initial interrupted attempt. Run B's uninterrupted invocation took **430.03 s (7 min 10 s)**. These are measured launcher wall times. Task durations below come from Nextflow start/completion timestamps; overlapping tasks must not be summed. Cached durations describe work done before interruption.

| Task | Run A seconds | A status | Run B seconds |
| --- | ---: | --- | ---: |
| VALIDATE_INPUTS (validate and split GSE inputs) | 0.138 | CACHED | 0.148 |
| ORTHOLOGY:ORTHOLOGY_PREPARE (download and prepare reference proteomes) | 18.165 | CACHED | 19.717 |
| RNA_DOWNLOAD_GSE (GSE74554) | 37.108 | CACHED | 65.862 |
| RNA_DOWNLOAD_GSE (GSE116952) | 31.990 | CACHED | 48.851 |
| RNA_ALIGN_GSE (GSE74554) | 197.686 | COMPLETED | 170.068 |
| ORTHOLOGY:ORTHOLOGY_PAIR (genome_000001) | 224.552 | COMPLETED | 200.159 |
| RNA_ANALYZE (GSE74554) | 18.507 | COMPLETED | 17.465 |
| ORTHOLOGY:ORTHOLOGY_MERGE (merge orthology pairs in reference order) | 14.138 | COMPLETED | 14.436 |
| MICROARRAY_GSE (GSE64918) | 33.416 | COMPLETED | 32.713 |
| MICROARRAY_GSE (GSE54577) | 40.135 | COMPLETED | 40.410 |
| RNA_ALIGN_GSE (GSE116952) | 137.108 | COMPLETED | 133.527 |
| RNA_ANALYZE (GSE116952) | 16.582 | COMPLETED | 15.461 |
| INTEGRATE_RESULTS (integrate all GSE results) | 26.099 | COMPLETED | 29.858 |
| VERIFY_OUTPUTS (verify integrated results) | 0.138 | COMPLETED | 0.186 |

The configuration allowed 6 CPUs / 16,000 MB globally, 2 CPUs / 6,000 MB per analysis task, and 1 CPU / 3,000 MB for orthology. The pinned amd64 Docker image ran under emulation on an ARM host. Nested-container CPU/RSS metrics were unavailable and were not interpreted as zero. These two runs do not establish a general performance improvement.

Host free space was sampled every 3 seconds. Run A started with 38.67 GiB free and reached a sampled minimum of 35.20 GiB; Run B started with 35.21 GiB and reached 31.91 GiB. Retained allocated footprints were approximately 3.04 GiB and 2.96 GiB respectively. These host readings include other activity and Docker storage behavior; retained directory size is not peak scratch space. The earlier failed investigation is excluded from these two figures and is still retained.

## Evidence

Local evidence is under `real-runs/2026-11-17-clean-resume-repro/`, outside the repository: `comparison.json`, `resume_validation.json`, `source_integrity.json`, `disk_summary.json`, partial-download results, the regression log, and each run's `full_validation.json`, `TASK_TIMINGS.md`, traces, interruption record and final table. The controller finished with status `complete`.

```bash
Rscript real-runs/2026-11-17-clean-resume-repro/compare.R real-runs/2026-11-17-clean-resume-repro
```
