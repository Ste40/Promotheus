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
