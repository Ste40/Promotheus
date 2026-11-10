# Real operon detection test — 10 November 2026

The updated `Operons_detection` function completed successfully on **GSE116952**, *Lactococcus lactis*, using four paired-end runs: SRR7511294–SRR7511297. Runs 1294/1295 form group 1; runs 1296/1297 form group 2.

| Measurement | Result |
| --- | --- |
| Elapsed time for detection and output checks | 789.0594 seconds (13 min 9 s) |
| Predicted multi-gene operons | 490 |
| Gene pairs predicted in the same operon | 904 |
| Resources | Docker limit: 2 CPUs, 8,000 MB; Java heap: 6,000 MB |
| Image | `ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7` |
| Platform | linux/amd64 under emulation on an ARM host |

The test reused downloaded reads and reference artifacts, but Rockhopper performed its own read alignment. Timing excludes original data downloads, Nextflow scheduling and the rest of the pipeline. The nonempty operon table, required columns and organism-level mode file were checked. Local evidence is retained in `real-runs/2026-11-10-operons/` outside the Git repository.

The first attempt exposed incorrect replicate grouping and was stopped. The corrected run explicitly passed two comma-separated replicate groups to Rockhopper. The adapter now matches reads to accession IDs, validates read pairing and group metadata, checks external command exit codes, and fails if operon detection is incomplete. Tests cover paired/single-end reads, noncontiguous groups, missing mates and invalid metadata; the full local regression suite passed.

Only one GSE was available, so the requested consensus threshold of 3 was reduced to 1 by the existing implementation. **Consensus across multiple GSEs, end-to-end Nextflow execution with operons enabled, and biological accuracy against a curated operon reference remain unvalidated.** These predictions are not a biological ground truth.
