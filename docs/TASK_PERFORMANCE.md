# Per-task performance report — complete Nextflow run

## What was measured

This report describes the **3 November 2026 full workflow**, not the smaller two-comparison benchmark. It used GSE8182 (six Microarray samples), GSE116952 (four RNA-Seq samples), reference ASM317683v1 and 128 archived genomes (127 comparisons). Operons were disabled. All **136 tasks completed without Nextflow cache reuse**; the final table contained 352 genes × 28 columns and passed the output/comparison audits.

The host used Docker Desktop with 8 CPUs and about 22 GiB RAM; Nextflow had an 8 CPU / 20,000 MB budget. The `standard,balanced` profile allowed four orthology pairs (each requesting 1 CPU / 3000 MB) and reserved capacity for a scientific task (2 CPUs / 8000 MB). The pinned image uses `linux/amd64`; these measurements are machine- and workload-specific, not native-Linux performance guarantees.

**Total task span: 3 h 42 min 29.9 s**, from 10:17:43 to 14:00:13 CET. Genome downloads were excluded because the archived genome snapshot was supplied. GEO/SRA retrieval, task staging, container startup, scientific calculations and output writing were included in task durations. Nextflow startup before the first task and the external final audit after the last task are excluded.

All durations are **wall-clock time = trace `complete − start`**, including waiting inside a running task. They are not CPU time. The `%cpu` and peak-RSS fields are unavailable for these nested Docker tasks on macOS; no measured CPU utilization or peak-memory claim is made. Requested resources are limits, not measured use.

## Duration by process

The span is from the first start to the last completion for a process. The sum adds individual task durations even when they overlap; neither sums nor process spans may be added to infer the workflow runtime.

| Process | Tasks | Minimum | Median | Maximum | Span | Sum of task durations | Peak simultaneous |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `VALIDATE_INPUTS` | 1 | 0.69 s | 0.69 s | 0.69 s | 0.69 s | 0.69 s | 1 |
| `ORTHOLOGY_PREPARE` | 1 | 1 min 57.3 s | 1 min 57.3 s | 1 min 57.3 s | 1 min 57.3 s | 1 min 57.3 s | 1 |
| `RNA_DOWNLOAD_GSE` | 1 | 33 min 16.6 s | 33 min 16.6 s | 33 min 16.6 s | 33 min 16.6 s | 33 min 16.6 s | 1 |
| `ORTHOLOGY_PAIR` | 127 | 1 min 47.9 s | 6 min 50.6 s | 10 min 57.8 s | 3 h 39 min 07.4 s | 14 h 31 min 20.2 s | 4 |
| `RNA_ALIGN_GSE` | 1 | 54 min 26.0 s | 54 min 26.0 s | 54 min 26.0 s | 54 min 26.0 s | 54 min 26.0 s | 1 |
| `RNA_ANALYZE` | 1 | 4 min 34.0 s | 4 min 34.0 s | 4 min 34.0 s | 4 min 34.0 s | 4 min 34.0 s | 1 |
| `ORTHOLOGY_MERGE` | 1 | 17.90 s | 17.90 s | 17.90 s | 17.90 s | 17.90 s | 1 |
| `MICROARRAY_GSE` | 1 | 35.13 s | 35.13 s | 35.13 s | 35.13 s | 35.13 s | 1 |
| `INTEGRATE_RESULTS` | 1 | 30.77 s | 30.77 s | 30.77 s | 30.77 s | 30.77 s | 1 |
| `VERIFY_OUTPUTS` | 1 | 0.11 s | 0.11 s | 0.11 s | 0.11 s | 0.11 s | 1 |

## Interpretation

- Orthology was the longest branch: its 127 pair tasks overlapped up to four at a time. Microarray analysis could only start after the orthology merge; final integration waited for both branches.
- RNA-Seq download finished at 10:51:00.748 and alignment started at 10:51:01.883: **1.135 seconds of downstream waiting**. RNA-Seq analysis finished at 11:50:02.729 while orthology continued. This overlap prevents adding the RNA and orthology times together.
- One aggregate alignment task covers all four RNA-Seq samples. This trace does not split time into individual samples, indexing, mapping, filtering or compression. A finer breakdown would require additional instrumentation.
- Individual genome comparisons have different workloads; a longer task is not evidence of a regression. The median and range describe this run, not variability across repeated runs.
- In the separate **4 November two-comparison test**, the orthology-only workflow took 525.784 s serially and 297.766 s with two concurrent pairs (1.77×). It used three archived genomes and five tasks per run, not all 128 genomes or the transcriptomics branches. See [method and limitations](PERFORMANCE.md#controlled-comparison-4-november-2026).

## Every task

Times are CET on 3 November. “Submit → start” measures only the interval recorded by the executor; it does **not** include time waiting for dependencies or for tasks to become eligible for submission. Genome names identify the other genome compared with ASM317683v1. Durations include milliseconds, while displayed clock times are rounded down to whole seconds.

<details>
<summary>All 136 completed tasks, ordered by start time</summary>

| ID | Process | Genome / GSE | Start | Complete | Duration | Submit → start | Requested CPU / MB |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: |
| 1 | `VALIDATE_INPUTS` | validate and split GSE inputs | 10:17:43 | 10:17:43 | 0.69 s | 0.073 s | 1 / 512 |
| 2 | `ORTHOLOGY_PREPARE` | download and prepare reference proteomes | 10:17:44 | 10:19:41 | 1 min 57.3 s | 0.030 s | 2 / 8000 |
| 3 | `RNA_DOWNLOAD_GSE` | GSE116952 | 10:17:44 | 10:51:00 | 33 min 16.6 s | 0.103 s | 2 / 8000 |
| 4 | `ORTHOLOGY_PAIR` | genome_000001: A12 | 10:19:41 | 10:26:32 | 6 min 50.6 s | 0.073 s | 1 / 3000 |
| 7 | `ORTHOLOGY_PAIR` | genome_000004: AF116-28_contig_1_np1 | 10:19:41 | 10:26:46 | 7 min 05.5 s | 0.061 s | 1 / 3000 |
| 8 | `ORTHOLOGY_PAIR` | genome_000005: AF116-29_contig_12_np1 | 10:19:41 | 10:26:42 | 7 min 00.9 s | 0.135 s | 1 / 3000 |
| 10 | `ORTHOLOGY_PAIR` | genome_000007: AF116-32_contig_2_np1 | 10:19:41 | 10:26:48 | 7 min 06.9 s | 0.145 s | 1 / 3000 |
| 9 | `ORTHOLOGY_PAIR` | genome_000006: AF116-30_contig_3_np1 | 10:26:32 | 10:32:43 | 6 min 11.8 s | 0.081 s | 1 / 3000 |
| 5 | `ORTHOLOGY_PAIR` | genome_000002: AF116-22 | 10:26:42 | 10:32:52 | 6 min 10.2 s | 0.091 s | 1 / 3000 |
| 6 | `ORTHOLOGY_PAIR` | genome_000003: AF116-24 | 10:26:46 | 10:33:02 | 6 min 15.4 s | 0.076 s | 1 / 3000 |
| 12 | `ORTHOLOGY_PAIR` | genome_000009: AF120-91-WCA | 10:26:48 | 10:32:48 | 6 min 00.4 s | 0.085 s | 1 / 3000 |
| 13 | `ORTHOLOGY_PAIR` | genome_000010: ASM1339501v3 | 10:32:43 | 10:39:10 | 6 min 27.0 s | 0.069 s | 1 / 3000 |
| 11 | `ORTHOLOGY_PAIR` | genome_000008: AF117-82-BHI | 10:32:49 | 10:39:03 | 6 min 14.3 s | 0.090 s | 1 / 3000 |
| 14 | `ORTHOLOGY_PAIR` | genome_000011: ASM1433471v1 | 10:32:52 | 10:39:04 | 6 min 12.1 s | 0.087 s | 1 / 3000 |
| 21 | `ORTHOLOGY_PAIR` | genome_000018: ASM1692099v2 | 10:33:02 | 10:39:08 | 6 min 06.5 s | 0.087 s | 1 / 3000 |
| 17 | `ORTHOLOGY_PAIR` | genome_000014: ASM1602881v1 | 10:39:03 | 10:47:05 | 8 min 02.4 s | 0.052 s | 1 / 3000 |
| 22 | `ORTHOLOGY_PAIR` | genome_000019: ASM1692101v2 | 10:39:05 | 10:47:21 | 8 min 16.8 s | 0.091 s | 1 / 3000 |
| 15 | `ORTHOLOGY_PAIR` | genome_000012: ASM1488460v1 | 10:39:09 | 10:47:24 | 8 min 15.4 s | 0.090 s | 1 / 3000 |
| 16 | `ORTHOLOGY_PAIR` | genome_000013: ASM1602797v1 | 10:39:11 | 10:47:09 | 7 min 58.4 s | 0.087 s | 1 / 3000 |
| 23 | `ORTHOLOGY_PAIR` | genome_000020: ASM1695295v2 | 10:47:05 | 10:56:16 | 9 min 11.0 s | 0.063 s | 1 / 3000 |
| 18 | `ORTHOLOGY_PAIR` | genome_000015: ASM1602883v1 | 10:47:09 | 10:56:23 | 9 min 13.6 s | 0.245 s | 1 / 3000 |
| 19 | `ORTHOLOGY_PAIR` | genome_000016: ASM1612713v1 | 10:47:21 | 10:56:43 | 9 min 21.1 s | 0.069 s | 1 / 3000 |
| 24 | `ORTHOLOGY_PAIR` | genome_000021: ASM1782119v2 | 10:47:24 | 10:57:37 | 10 min 13.0 s | 0.093 s | 1 / 3000 |
| 131 | `RNA_ALIGN_GSE` | GSE116952 | 10:51:01 | 11:45:27 | 54 min 26.0 s | 0.112 s | 2 / 8000 |
| 25 | `ORTHOLOGY_PAIR` | genome_000022: ASM1786941v1 | 10:56:17 | 11:04:53 | 8 min 36.5 s | 0.087 s | 1 / 3000 |
| 20 | `ORTHOLOGY_PAIR` | genome_000017: ASM1664919v2 | 10:56:23 | 11:06:17 | 9 min 54.6 s | 0.074 s | 1 / 3000 |
| 27 | `ORTHOLOGY_PAIR` | genome_000024: ASM2022175v1 | 10:56:43 | 11:05:37 | 8 min 54.4 s | 0.077 s | 1 / 3000 |
| 26 | `ORTHOLOGY_PAIR` | genome_000023: ASM19270v1 | 10:57:37 | 11:06:56 | 9 min 18.7 s | 0.075 s | 1 / 3000 |
| 32 | `ORTHOLOGY_PAIR` | genome_000029: ASM207861v4 | 11:04:53 | 11:14:00 | 9 min 06.8 s | 0.065 s | 1 / 3000 |
| 31 | `ORTHOLOGY_PAIR` | genome_000028: ASM207849v2 | 11:05:37 | 11:14:45 | 9 min 08.0 s | 0.074 s | 1 / 3000 |
| 28 | `ORTHOLOGY_PAIR` | genome_000025: ASM2046375v1 | 11:06:18 | 11:15:31 | 9 min 13.0 s | 0.067 s | 1 / 3000 |
| 33 | `ORTHOLOGY_PAIR` | genome_000030: ASM207885v1 | 11:06:56 | 11:15:22 | 8 min 26.2 s | 0.065 s | 1 / 3000 |
| 30 | `ORTHOLOGY_PAIR` | genome_000027: ASM207847v4 | 11:14:00 | 11:23:14 | 9 min 14.1 s | 0.124 s | 1 / 3000 |
| 29 | `ORTHOLOGY_PAIR` | genome_000026: ASM207843v1 | 11:14:45 | 11:25:14 | 10 min 28.9 s | 0.056 s | 1 / 3000 |
| 36 | `ORTHOLOGY_PAIR` | genome_000033: ASM214821v1 | 11:15:22 | 11:25:47 | 10 min 24.6 s | 0.077 s | 1 / 3000 |
| 34 | `ORTHOLOGY_PAIR` | genome_000031: ASM207897v5 | 11:15:31 | 11:26:03 | 10 min 32.3 s | 0.027 s | 1 / 3000 |
| 37 | `ORTHOLOGY_PAIR` | genome_000034: ASM2249490v1 | 11:23:15 | 11:26:46 | 3 min 31.4 s | 0.231 s | 1 / 3000 |
| 41 | `ORTHOLOGY_PAIR` | genome_000038: ASM2320675v1 | 11:25:14 | 11:35:26 | 10 min 11.8 s | 0.072 s | 1 / 3000 |
| 35 | `ORTHOLOGY_PAIR` | genome_000032: ASM207899v4 | 11:25:47 | 11:36:39 | 10 min 52.1 s | 0.086 s | 1 / 3000 |
| 38 | `ORTHOLOGY_PAIR` | genome_000035: ASM2267081v1 | 11:26:03 | 11:36:47 | 10 min 43.5 s | 0.061 s | 1 / 3000 |
| 40 | `ORTHOLOGY_PAIR` | genome_000037: ASM2307851v1 | 11:26:46 | 11:37:44 | 10 min 57.8 s | 0.072 s | 1 / 3000 |
| 39 | `ORTHOLOGY_PAIR` | genome_000036: ASM2280975v1 | 11:35:26 | 11:45:28 | 10 min 01.7 s | 0.121 s | 1 / 3000 |
| 48 | `ORTHOLOGY_PAIR` | genome_000045: ASM2382220v1 | 11:36:39 | 11:47:02 | 10 min 23.2 s | 0.081 s | 1 / 3000 |
| 43 | `ORTHOLOGY_PAIR` | genome_000040: ASM2361188v1 | 11:36:47 | 11:46:48 | 10 min 00.9 s | 0.096 s | 1 / 3000 |
| 47 | `ORTHOLOGY_PAIR` | genome_000044: ASM2382218v1 | 11:37:44 | 11:48:36 | 10 min 51.9 s | 0.070 s | 1 / 3000 |
| 42 | `ORTHOLOGY_PAIR` | genome_000039: ASM2361186v3 | 11:45:28 | 11:53:31 | 8 min 02.5 s | 0.133 s | 1 / 3000 |
| 132 | `RNA_ANALYZE` | GSE116952 | 11:45:28 | 11:50:02 | 4 min 34.0 s | 0.105 s | 2 / 8000 |
| 46 | `ORTHOLOGY_PAIR` | genome_000043: ASM2380974v1 | 11:46:48 | 11:54:31 | 7 min 42.9 s | 0.122 s | 1 / 3000 |
| 45 | `ORTHOLOGY_PAIR` | genome_000042: ASM2373475v1 | 11:47:02 | 11:54:44 | 7 min 41.0 s | 0.070 s | 1 / 3000 |
| 44 | `ORTHOLOGY_PAIR` | genome_000041: ASM2373473v1 | 11:48:36 | 11:55:57 | 7 min 20.8 s | 0.051 s | 1 / 3000 |
| 49 | `ORTHOLOGY_PAIR` | genome_000046: ASM2420534v1 | 11:53:31 | 12:00:11 | 6 min 40.1 s | 0.071 s | 1 / 3000 |
| 54 | `ORTHOLOGY_PAIR` | genome_000051: ASM2759476v1 | 11:54:31 | 12:01:27 | 6 min 56.1 s | 0.066 s | 1 / 3000 |
| 50 | `ORTHOLOGY_PAIR` | genome_000047: ASM2466563v1 | 11:54:44 | 11:56:49 | 2 min 05.6 s | 0.094 s | 1 / 3000 |
| 51 | `ORTHOLOGY_PAIR` | genome_000048: ASM2504v1 | 11:55:57 | 12:03:06 | 7 min 09.3 s | 0.080 s | 1 / 3000 |
| 55 | `ORTHOLOGY_PAIR` | genome_000052: ASM280418v1 | 11:56:49 | 12:03:43 | 6 min 53.2 s | 0.101 s | 1 / 3000 |
| 53 | `ORTHOLOGY_PAIR` | genome_000050: ASM2556295v3 | 12:00:11 | 12:02:20 | 2 min 09.4 s | 0.080 s | 1 / 3000 |
| 52 | `ORTHOLOGY_PAIR` | genome_000049: ASM2556289v2 | 12:01:27 | 12:03:35 | 2 min 08.5 s | 0.081 s | 1 / 3000 |
| 56 | `ORTHOLOGY_PAIR` | genome_000053: ASM280428v1 | 12:02:21 | 12:09:10 | 6 min 49.0 s | 0.069 s | 1 / 3000 |
| 57 | `ORTHOLOGY_PAIR` | genome_000054: ASM289522v1 | 12:03:06 | 12:10:31 | 7 min 24.4 s | 0.069 s | 1 / 3000 |
| 58 | `ORTHOLOGY_PAIR` | genome_000055: ASM2899423v1 | 12:03:36 | 12:11:24 | 7 min 48.5 s | 0.087 s | 1 / 3000 |
| 59 | `ORTHOLOGY_PAIR` | genome_000056: ASM2902386v1 | 12:03:43 | 12:11:13 | 7 min 30.3 s | 0.098 s | 1 / 3000 |
| 64 | `ORTHOLOGY_PAIR` | genome_000061: ASM3048966v1 | 12:09:10 | 12:17:06 | 7 min 56.2 s | 0.046 s | 1 / 3000 |
| 61 | `ORTHOLOGY_PAIR` | genome_000058: ASM2958111v1 | 12:10:31 | 12:18:20 | 7 min 49.0 s | 0.045 s | 1 / 3000 |
| 60 | `ORTHOLOGY_PAIR` | genome_000057: ASM2954222v1 | 12:11:13 | 12:13:29 | 2 min 15.6 s | 0.078 s | 1 / 3000 |
| 62 | `ORTHOLOGY_PAIR` | genome_000059: ASM3043797v1 | 12:11:24 | 12:19:32 | 8 min 07.4 s | 0.098 s | 1 / 3000 |
| 63 | `ORTHOLOGY_PAIR` | genome_000060: ASM3048046v1 | 12:13:29 | 12:21:28 | 7 min 59.1 s | 0.154 s | 1 / 3000 |
| 65 | `ORTHOLOGY_PAIR` | genome_000062: ASM3287881v1 | 12:17:06 | 12:23:58 | 6 min 52.2 s | 0.065 s | 1 / 3000 |
| 66 | `ORTHOLOGY_PAIR` | genome_000063: ASM3363120v1 | 12:18:20 | 12:20:11 | 1 min 50.8 s | 0.070 s | 1 / 3000 |
| 72 | `ORTHOLOGY_PAIR` | genome_000069: ASM372227v1 | 12:19:32 | 12:25:34 | 6 min 02.7 s | 0.070 s | 1 / 3000 |
| 71 | `ORTHOLOGY_PAIR` | genome_000068: ASM372225v1 | 12:20:11 | 12:26:06 | 5 min 54.9 s | 0.086 s | 1 / 3000 |
| 68 | `ORTHOLOGY_PAIR` | genome_000065: ASM34457v1 | 12:21:28 | 12:27:42 | 6 min 14.1 s | 0.079 s | 1 / 3000 |
| 73 | `ORTHOLOGY_PAIR` | genome_000070: ASM3872123v1 | 12:23:58 | 12:25:58 | 1 min 59.3 s | 0.081 s | 1 / 3000 |
| 67 | `ORTHOLOGY_PAIR` | genome_000064: ASM3387991v1 | 12:25:34 | 12:32:22 | 6 min 47.6 s | 0.090 s | 1 / 3000 |
| 69 | `ORTHOLOGY_PAIR` | genome_000066: ASM3598596v1 | 12:25:58 | 12:32:46 | 6 min 48.2 s | 0.091 s | 1 / 3000 |
| 70 | `ORTHOLOGY_PAIR` | genome_000067: ASM362739v2 | 12:26:06 | 12:33:02 | 6 min 55.7 s | 0.091 s | 1 / 3000 |
| 79 | `ORTHOLOGY_PAIR` | genome_000076: ASM4355164v1 | 12:27:42 | 12:29:47 | 2 min 04.4 s | 0.088 s | 1 / 3000 |
| 78 | `ORTHOLOGY_PAIR` | genome_000075: ASM4355085v1 | 12:29:47 | 12:37:00 | 7 min 13.2 s | 0.162 s | 1 / 3000 |
| 76 | `ORTHOLOGY_PAIR` | genome_000073: ASM4154546v1 | 12:32:22 | 12:39:23 | 7 min 00.4 s | 0.050 s | 1 / 3000 |
| 75 | `ORTHOLOGY_PAIR` | genome_000072: ASM4051602v1 | 12:32:46 | 12:39:51 | 7 min 05.5 s | 0.084 s | 1 / 3000 |
| 77 | `ORTHOLOGY_PAIR` | genome_000074: ASM4187943v1 | 12:33:02 | 12:35:01 | 1 min 59.8 s | 0.089 s | 1 / 3000 |
| 74 | `ORTHOLOGY_PAIR` | genome_000071: ASM4044831v1 | 12:35:02 | 12:42:13 | 7 min 11.7 s | 0.106 s | 1 / 3000 |
| 81 | `ORTHOLOGY_PAIR` | genome_000078: ASM4433675v1 | 12:37:00 | 12:39:06 | 2 min 05.7 s | 0.069 s | 1 / 3000 |
| 84 | `ORTHOLOGY_PAIR` | genome_000081: ASM4436055v1 | 12:39:06 | 12:41:04 | 1 min 58.2 s | 0.078 s | 1 / 3000 |
| 80 | `ORTHOLOGY_PAIR` | genome_000077: ASM4355239v1 | 12:39:23 | 12:42:19 | 2 min 56.2 s | 0.095 s | 1 / 3000 |
| 83 | `ORTHOLOGY_PAIR` | genome_000080: ASM4435830v1 | 12:39:52 | 12:46:45 | 6 min 53.4 s | 0.091 s | 1 / 3000 |
| 82 | `ORTHOLOGY_PAIR` | genome_000079: ASM4433739v1 | 12:41:04 | 12:48:12 | 7 min 07.9 s | 0.081 s | 1 / 3000 |
| 87 | `ORTHOLOGY_PAIR` | genome_000084: ASM4486679v1 | 12:42:14 | 12:49:04 | 6 min 50.9 s | 0.078 s | 1 / 3000 |
| 88 | `ORTHOLOGY_PAIR` | genome_000085: ASM4515977v1 | 12:42:19 | 12:49:09 | 6 min 49.9 s | 0.090 s | 1 / 3000 |
| 90 | `ORTHOLOGY_PAIR` | genome_000087: ASM4534840v1 | 12:46:45 | 12:54:06 | 7 min 20.6 s | 0.087 s | 1 / 3000 |
| 85 | `ORTHOLOGY_PAIR` | genome_000082: ASM4436367v1 | 12:48:12 | 12:55:03 | 6 min 51.0 s | 0.082 s | 1 / 3000 |
| 86 | `ORTHOLOGY_PAIR` | genome_000083: ASM4436571v1 | 12:49:05 | 12:56:20 | 7 min 15.5 s | 0.071 s | 1 / 3000 |
| 91 | `ORTHOLOGY_PAIR` | genome_000088: ASM4534841v1 | 12:49:09 | 12:57:20 | 8 min 11.3 s | 0.096 s | 1 / 3000 |
| 89 | `ORTHOLOGY_PAIR` | genome_000086: ASM4534839v1 | 12:54:06 | 13:01:35 | 7 min 28.9 s | 0.084 s | 1 / 3000 |
| 93 | `ORTHOLOGY_PAIR` | genome_000090: ASM4704114v1 | 12:55:03 | 13:02:21 | 7 min 17.1 s | 0.089 s | 1 / 3000 |
| 94 | `ORTHOLOGY_PAIR` | genome_000091: ASM47825v2 | 12:56:20 | 13:03:30 | 7 min 09.6 s | 0.072 s | 1 / 3000 |
| 95 | `ORTHOLOGY_PAIR` | genome_000092: ASM5036884v1 | 12:57:20 | 13:04:27 | 7 min 06.8 s | 0.089 s | 1 / 3000 |
| 92 | `ORTHOLOGY_PAIR` | genome_000089: ASM4623761v1 | 13:01:35 | 13:08:26 | 6 min 50.8 s | 0.082 s | 1 / 3000 |
| 98 | `ORTHOLOGY_PAIR` | genome_000095: ASM5082238v1 | 13:02:21 | 13:08:49 | 6 min 28.7 s | 0.059 s | 1 / 3000 |
| 101 | `ORTHOLOGY_PAIR` | genome_000098: ASM5381886v1 | 13:03:30 | 13:09:51 | 6 min 21.2 s | 0.086 s | 1 / 3000 |
| 96 | `ORTHOLOGY_PAIR` | genome_000093: ASM5042958v1 | 13:04:27 | 13:10:49 | 6 min 21.8 s | 0.077 s | 1 / 3000 |
| 97 | `ORTHOLOGY_PAIR` | genome_000094: ASM5051145v1 | 13:08:26 | 13:15:07 | 6 min 41.6 s | 0.082 s | 1 / 3000 |
| 99 | `ORTHOLOGY_PAIR` | genome_000096: ASM5082239v1 | 13:08:49 | 13:15:02 | 6 min 12.4 s | 0.093 s | 1 / 3000 |
| 102 | `ORTHOLOGY_PAIR` | genome_000099: ASM5517987v1 | 13:09:51 | 13:16:12 | 6 min 21.5 s | 0.085 s | 1 / 3000 |
| 105 | `ORTHOLOGY_PAIR` | genome_000102: ASM5518025v1 | 13:10:49 | 13:12:42 | 1 min 53.0 s | 0.091 s | 1 / 3000 |
| 100 | `ORTHOLOGY_PAIR` | genome_000097: ASM5113114v1 | 13:12:42 | 13:18:58 | 6 min 16.1 s | 0.185 s | 1 / 3000 |
| 103 | `ORTHOLOGY_PAIR` | genome_000100: ASM5517988v1 | 13:15:02 | 13:21:14 | 6 min 12.5 s | 0.070 s | 1 / 3000 |
| 104 | `ORTHOLOGY_PAIR` | genome_000101: ASM5518009v1 | 13:15:07 | 13:21:38 | 6 min 30.6 s | 0.098 s | 1 / 3000 |
| 111 | `ORTHOLOGY_PAIR` | genome_000108: ASM5538869v1 | 13:16:13 | 13:22:31 | 6 min 18.0 s | 0.073 s | 1 / 3000 |
| 110 | `ORTHOLOGY_PAIR` | genome_000107: ASM5518059v1 | 13:18:59 | 13:25:07 | 6 min 08.4 s | 0.076 s | 1 / 3000 |
| 108 | `ORTHOLOGY_PAIR` | genome_000105: ASM5518047v1 | 13:21:15 | 13:27:22 | 6 min 07.5 s | 0.085 s | 1 / 3000 |
| 107 | `ORTHOLOGY_PAIR` | genome_000104: ASM5518035v1 | 13:21:38 | 13:27:49 | 6 min 11.1 s | 0.095 s | 1 / 3000 |
| 112 | `ORTHOLOGY_PAIR` | genome_000109: ASM5651745v1 | 13:22:31 | 13:29:01 | 6 min 30.2 s | 0.082 s | 1 / 3000 |
| 106 | `ORTHOLOGY_PAIR` | genome_000103: ASM5518031v1 | 13:25:07 | 13:31:37 | 6 min 29.8 s | 0.092 s | 1 / 3000 |
| 109 | `ORTHOLOGY_PAIR` | genome_000106: ASM5518057v1 | 13:27:22 | 13:34:00 | 6 min 37.6 s | 0.087 s | 1 / 3000 |
| 116 | `ORTHOLOGY_PAIR` | genome_000113: ASM6030721v1 | 13:27:49 | 13:34:04 | 6 min 15.1 s | 0.092 s | 1 / 3000 |
| 113 | `ORTHOLOGY_PAIR` | genome_000110: ASM5744397v1 | 13:29:01 | 13:35:33 | 6 min 32.0 s | 0.082 s | 1 / 3000 |
| 114 | `ORTHOLOGY_PAIR` | genome_000111: ASM5835280v1 | 13:31:37 | 13:38:03 | 6 min 26.4 s | 0.043 s | 1 / 3000 |
| 115 | `ORTHOLOGY_PAIR` | genome_000112: ASM5835312v1 | 13:34:00 | 13:40:34 | 6 min 33.8 s | 0.064 s | 1 / 3000 |
| 117 | `ORTHOLOGY_PAIR` | genome_000114: ASM6034408v1 | 13:34:05 | 13:40:28 | 6 min 23.6 s | 0.092 s | 1 / 3000 |
| 119 | `ORTHOLOGY_PAIR` | genome_000116: ASM6132581v1 | 13:35:33 | 13:41:52 | 6 min 19.0 s | 0.081 s | 1 / 3000 |
| 118 | `ORTHOLOGY_PAIR` | genome_000115: ASM6058836v1 | 13:38:04 | 13:44:23 | 6 min 19.3 s | 0.078 s | 1 / 3000 |
| 122 | `ORTHOLOGY_PAIR` | genome_000119: ASM795476v1 | 13:40:28 | 13:46:46 | 6 min 17.8 s | 0.074 s | 1 / 3000 |
| 120 | `ORTHOLOGY_PAIR` | genome_000117: ASM686v1 | 13:40:34 | 13:46:29 | 5 min 54.8 s | 0.091 s | 1 / 3000 |
| 125 | `ORTHOLOGY_PAIR` | genome_000122: LAC-D-17 | 13:41:52 | 13:43:40 | 1 min 47.9 s | 0.089 s | 1 / 3000 |
| 126 | `ORTHOLOGY_PAIR` | genome_000123: Lactococcus_lactis_CIRM-BIA2553 | 13:43:40 | 13:49:50 | 6 min 10.1 s | 0.086 s | 1 / 3000 |
| 127 | `ORTHOLOGY_PAIR` | genome_000124: LlIMDOWA12L8_assembly | 13:44:23 | 13:50:58 | 6 min 35.3 s | 0.082 s | 1 / 3000 |
| 121 | `ORTHOLOGY_PAIR` | genome_000118: ASM76111v1 | 13:46:29 | 13:52:43 | 6 min 14.1 s | 0.083 s | 1 / 3000 |
| 123 | `ORTHOLOGY_PAIR` | genome_000120: ASM80737v1 | 13:46:46 | 13:53:31 | 6 min 45.3 s | 0.068 s | 1 / 3000 |
| 124 | `ORTHOLOGY_PAIR` | genome_000121: ASM872887v1 | 13:49:50 | 13:56:07 | 6 min 16.0 s | 0.079 s | 1 / 3000 |
| 129 | `ORTHOLOGY_PAIR` | genome_000126: PacBio_sequel_II | 13:50:58 | 13:57:08 | 6 min 09.8 s | 0.080 s | 1 / 3000 |
| 128 | `ORTHOLOGY_PAIR` | genome_000125: Lla1834_v1_all | 13:52:43 | 13:58:15 | 5 min 32.2 s | 0.078 s | 1 / 3000 |
| 130 | `ORTHOLOGY_PAIR` | genome_000127: Varani_LAB_06 | 13:53:32 | 13:58:48 | 5 min 16.7 s | 0.072 s | 1 / 3000 |
| 133 | `ORTHOLOGY_MERGE` | merge orthology pairs in reference order | 13:58:48 | 13:59:06 | 17.90 s | 0.032 s | 2 / 8000 |
| 134 | `MICROARRAY_GSE` | GSE8182 | 13:59:06 | 13:59:41 | 35.13 s | 0.092 s | 2 / 8000 |
| 135 | `INTEGRATE_RESULTS` | integrate all GSE results | 13:59:42 | 14:00:12 | 30.77 s | 0.074 s | 2 / 8000 |
| 136 | `VERIFY_OUTPUTS` | verify integrated results | 14:00:12 | 14:00:13 | 0.11 s | 0.090 s | 1 / 512 |

</details>

The [machine-readable report](benchmarks/2026-11-03-task-timings.json) retains exact timestamps, durations, per-stage statistics and the SHA-256 of the original Nextflow trace. Scientific code: commit `5cf286b`.
