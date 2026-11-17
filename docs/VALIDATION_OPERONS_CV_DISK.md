# Focused validation: operons, CV and RNA-Seq disk use

## Two-GSE operon detection

Real reads from GSE116952 (four paired-end samples) and GSE74554 / SRP065562 (six single-end samples) were truncated to the **first 100,000 reads per FASTQ**; paired mates were retained in order. The second study has three early-stationary and three late-stationary replicates, verified against ENA/GEO metadata. Both used the existing ASM317683v1 reference, 2 CPUs and an 8,000 MB container limit (6,000 MB Java heap).

Rockhopper completed in 44.4 s and 64.1 s respectively. It produced 480 and 525 operon rows. With threshold 2, the mode file was absent after the first GSE and present after the second: 593 distinct gene lists, of which 412 occurred in both GSEs.

**Threshold semantics:** the parameter controls how many GSEs must be available before producing a modal table. The table is a union of gene lists with modal coordinates, not an intersection or a requirement for support by two studies. The existing policy was retained.

This is a component smoke test on real subsampled reads, not a complete Nextflow run with operons enabled. Truncation and mapping different strains to a shared reference limit biological interpretation. It does not establish prediction accuracy.

Targeted regression tests reproduced and corrected:

- Earlier waiting expression tables were annotated using the latest GSE's operons. Each table now uses its own GSE accession.
- Consensus contained comma-separated gene lists but was compared directly with single gene IDs. Both immediate and deferred annotation now identify the strand-aware first gene.
- An existing mode file caused an early return before copying later GSE operon tables.
- Empty operon predictions now retain character gene IDs and use zero-safe loops. Reverse-strand annotation is covered by regression tests.
- More cached GSEs than the threshold prevented mode generation (`==`); reaching or exceeding the threshold now works.

Run `Rscript tests/test_operons_multigse.R .` and `bash tests/shell/test_all.sh` for regression checks.

## CV after normalization

The real interspecies matrix contained 345 genes and 10 samples. All means were positive (minimum 3.7511), all expression values and CVs were finite, and the largest CV was 0.8736. A separate numerical stress test applied mean-only ComBat across the two GSEs; it also produced no nonpositive means or non-finite values. This cross-species/technology stress test is not a recommended biological correction design.

Artificial zero/negative means, NA and infinite values were rejected explicitly. Positive scaling preserved CV, and constant positive expression produced CV zero. CV remains SD/mean on the integrated expression scale, not on reconstructed unlogged abundance.

## RNA-Seq temporary-space measurement

One full paired-end accession, SRR7511294 (10,658,636 spots), was extracted and compressed with the production R helpers in the pinned container (2 CPUs, 4,000 MB). A host monitor sampled files and free disk space every 3 seconds.

- Extraction: 380.5 seconds; extraction plus compression: 424.4 seconds.
- Largest sampled output/scratch footprint: **9.46 GiB** (sampling may miss a brief higher peak).
- Final compressed FASTQs: approximately **1.45 GiB**.
- Host free space fell from **64.38 to 52.32 GiB**, a 12.05 GiB change. It had not recovered at the immediate post-container check despite only compressed reads remaining in the output directory.
- No residual fasterq scratch directory was present in the final output; the container reported only 16 KiB in `/tmp` and 8 KiB in `/root/.ncbi`.

The test passed, but it does **not** establish that the disk-exhaustion risk is fixed. Logical file sizes alone underestimate the observed host-space requirement on this setup. The cause of delayed/missing host-space reclamation remains unproven; these figures should not be generalized to other machines or datasets. Keep sequential accession preparation and substantial headroom; the earlier four-accession run succeeded after starting with roughly 66 GiB free. A complete low-disk recovery test and a storage-management change remain outstanding.

## Evidence

Local inputs, scripts, logs and machine-readable reports are retained under `real-runs/2026-11-13-validation-234/` outside the repository. No historical benchmark outputs were overwritten.

## Follow-up: isolated scratch and incomplete FASTQ pairs

The download helper now passes an explicit output directory and creates a separate scratch directory for each accession attempt. That scratch directory is removed after command success or failure, with an R exit handler as an additional safeguard. Unrelated files are preserved. This handles ordinary errors; it cannot guarantee cleanup after SIGKILL or a machine crash.

Completeness checks now accept exactly one single-end file or both paired-end mates. An isolated `_1` or `_2` mate is rejected instead of being mistaken for a complete single-end accession. These changes do not add cross-task caching: a failed Nextflow GSE task still needs its own recovery strategy.

A real fasterq safety-check test with explicit one-byte output and scratch limits failed as expected (`disk-limit exceeded`, exit 3), and left no files in the test directory. This exercises the tool's preflight failure path without filling the physical disk. Injected failures additionally verify scratch cleanup between retries and preservation of unrelated files. The complete regression suite passed.

The repeated full SRR7511294 test passed: extraction 376.7 s, extraction plus compression 421.9 s, largest sampled output/scratch footprint **8.19 GiB** versus 9.46 GiB previously. Final output contains only the two compressed mates. A single sequential comparison with 3-second sampling does not establish a stable peak reduction or speedup. The host returned to approximately 60 GiB free after completion on this run; the earlier delayed reclamation has not been conclusively explained. Keep adequate disk headroom rather than sizing storage only from final compressed FASTQs. Evidence is in `real-runs/2026-11-16-disk-check/disk/`.

The decompressed FASTQ SHA-256 hashes match the previous benchmark for both mates, confirming unchanged read content in this comparison.

## Complete Nextflow operon workflow

The complete scientific path passed with operons enabled and threshold 2, using GSE116952 (4 paired-end samples), GSE74554 (6 single-end samples), and Microarray GSE54577 (6 samples). RNA inputs were the previously prepared first 100,000 reads per FASTQ. The Microarray task succeeded in the first attempt; after NCBI returned HTTP 403 for supplementary files on a later attempt, its real full-table checkpoint was reused. The orthology map was verified byte-identical before reusing this checkpoint. This is not a fresh full-data download benchmark.

The successful attempt completed 12 Nextflow tasks in **480.95 seconds (8 min 1 s)**. Both RNA GSEs were grouped into one organism-level `RNA_ANALYZE` task. Final validation and an additional audit passed: **232 genes × 40 columns**, all **16 samples** present (6 + 4 + 6), and no `WAITING_FOR_OPERONS_TABLE` rows in either RNA result. The mode table was retained as evidence. These checks validate execution and annotation completeness, not biological accuracy of operon prediction.

| Task | Observed time |
| --- | --- |
| Orthology comparison | 3 min 7 s |
| RNA alignment, GSE116952 | 2 min 1 s |
| RNA alignment, GSE74554 | 2 min 40 s |
| Joint RNA analysis and operons | 2 min 8 s |
| Microarray checkpoint reuse | 19 s |
| Integration | 24 s |

The first attempt exposed a staging conflict in `my_index.files`: independently built Rsubread indexes contain different temporary filenames. RNA analysis no longer stages root-level `my_index.*` alignment artifacts, which it does not consume. Reference genomes remain subject to byte-level conflict checks. Regression tests verify both index exclusion and rejection of genuinely different reference sequences; the full regression suite passed.

Evidence, per-task timings, source hashes, input preparation details and preserved failed-attempt logs are in `real-runs/2026-11-16-operons-nextflow/`. Download stages consumed prepared read subsets; their timing does not measure downloading full FASTQs. Total time excludes failed attempts and preparation.
