# Promotheus Nextflow compatibility and validation

> This document describes the former three-process workflow, now preserved as
> `main_legacy.nf`. The default `main.nf` uses isolated scientific stages. See
> [Modular Nextflow execution](MODULAR_NEXTFLOW.md). The historical real-data
> evidence here does not validate the new adapters; those have separate mock
> tests and a successful two-GSE real run documented in the modular guide.
> A formal real-data equivalence comparison is still required.

## Scientific behavior of the original pipeline

`Scripts/Introduction.sh` validates environment variables and invokes `Scripts/Microarray/main_pipeline.R`. That script is the true top-level orchestrator:

1. It builds an orthology map across the reference assemblies with RefSeq downloads and OrthoFinder.
2. The Microarray branch retrieves GEO metadata and raw files, handles Affymetrix, one-channel and two-channel platforms, performs within-experiment analysis, and maps probes to orthologous loci.
3. The RNA-Seq branch retrieves GEO/SRA metadata and reads, aligns them with Rsubread, parses and filters BAM files with samtools, optionally calls operons with Rockhopper, and identifies non-differentially expressed genes with DESeq2.
4. `interexp_analysis.R` intersects the mapped loci across experiments, applies ComBat batch correction, computes expression and percentile summaries, propagates differential-expression flags, and extracts 300 bp upstream promoter sequences.

The principal results are `/output/results/final_table.xlsx` and its equivalent TSV export; per-experiment tables are retained under `/output/results/fulltab/`.

## Compatibility strategy

The Nextflow workflow deliberately does not reimplement the complete R pipeline. Rewriting roughly 4,500 lines of stateful R/Bash would create a second scientific implementation and require a much larger numerical-validation study. Instead, `RUN_PROMOTHEUS` invokes the same published image (pinned to digest `sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7`) and passes the same ten scientific environment variables. It mounts ten repository scripts read-only over the image copies to normalize current RefSeq `ftp_path` values, handle empty filtered orthogroups, harden FASTQ downloads, correctly support disabled operon detection, reuse completed Microarray work, and accelerate independent Microarray and RNA-Seq work. Regression tests compare the optimized statistics, Kruskal-Wallis p-values and percentiles with the legacy implementation.

| Nextflow parameter | Legacy environment variable |
| --- | --- |
| `threads` | `THREADS` |
| `memory` | `MEMORY` |
| `correction_de` | `CORRECTION_DE` |
| `pvalue` | `PVALUE` |
| `operons_detection` | `OPERONS_DETECTION` |
| `amino_diff` | `AMINO_DIFF` |
| `similarity` | `SIMILARITY` |
| `asmcode` | `ASMcode` |
| `operons_threshold` | `OPERONS_THRESHOLD` |
| `multitaxa` | `MULTITAXA` |

The execution-only controls `max_gse_parallel`, `cache_dir`, `keep_cache`, `export_tsv`, and `container_platform` do not change statistical thresholds or biological algorithms. GSEs for the same organism remain ordered when operon-mode aggregation requires them; independent organisms and Microarray GSE calculations may run concurrently. Worker CPU and Rockhopper memory are divided from the totals declared by `threads` and `memory`.

Both interfaces mount the selected input directory at `/input` and result directory at `/output`. The Nextflow runner additionally mounts the staged, checksum-tracked compatibility scripts. Therefore, equivalence requires the same image digest, repository revision, inputs, parameters, network responses and tool nondeterminism.

## What Nextflow adds

- Structural validation of both Excel inputs before expensive downloads begin.
- SHA-256 manifests for inputs and key outputs.
- Content-addressed caching keyed by inputs, scientific parameters and overlay checksums.
- Per-accession full-tab checkpoints for auditable restart decisions.
- Explicit CPU and memory requests visible to an executor.
- Immediate propagation of Docker and validation failures.
- Timestamped trace, resource report, timeline and DAG artifacts under `pipeline_info/`; the trace includes per-process duration, CPU, peak memory, I/O and work directory.
- A deterministic no-network profile for CI.
- A post-run contract for the final table and promoter sequences.
- A TSV copy of the final Excel table with consistency validation.

These changes improve reproducibility, restart granularity and failure detection. Parallel workers reduce wall time when independent GSEs or SRA accessions are available, but they do not make OrthoFinder, alignment, Rockhopper, DESeq2 or ComBat intrinsically faster. Runtime remains sensitive to public-data servers, disk throughput and the AMD64 emulation required by the pinned image on ARM64 hosts.

## Automated coverage

`tests/shell/test_all.sh` performs:

1. Python syntax compilation and Bash syntax checks.
2. Positive validation of the two synthetic 29-column toy inputs.
3. Compatibility of all toy workbooks with the same R `openxlsx` reader used by Promotheus, when R is available.
4. Negative validation of malformed `ASMcode` and a missing RNA-Seq workbook.
5. A Docker-runner contract test for image, platform, scientific/execution variables, input/output/cache mounts and run record.
6. A complete DSL2 run using the mock engine.
7. Assertions for every propagated scientific parameter.
8. Byte-for-byte comparison of the mock published result with the formula-verified toy result.
9. Post-run schema validation, Excel/TSV consistency, accession checkpoints and key artifact checksums.

The GitHub Actions workflow runs the same suite on pushes and pull requests.

`tests/shell/test_container.sh` is an opt-in Docker test. It checks the pinned digest, architecture, entrypoint, bioinformatics executables, key R packages, byte identity of unchanged critical orchestration scripts, and execution of the compatibility overlays inside the container.

## Real-run evidence and remaining limitation

A real Docker-backed run of the supplied Lactococcus input snapshot completed `RUN_PROMOTHEUS` on 9 September 2026. The repaired workbook was structurally valid and contained 339 genes and 136 columns. The legacy Microarray code processed six of the 18 requested GSE accessions and skipped 12 whose current GEO supplementary files were absent or unreadable. The strict verifier detected the missing expression/percentile columns and rejected the partial result, as designed. See the README's validation-status section for the accession lists.

This run demonstrates real container execution, compatibility overlays, workbook repair and fail-closed validation. It does not establish full biological equivalence for every accession in the supplied snapshot.

## Evidence required for a formal biological-equivalence claim

The mock profile tests orchestration, not complete R/Bash numerical behavior. A formal equivalence study should run both entrypoints using the same immutable image digest and archived input/data snapshot, then compare:

- `orthologous_genes.txt` after stable sorting;
- every table in `results/fulltab/` with numeric tolerances appropriate to the R packages;
- the final gene set and column set;
- numeric expression, CV and percentile columns;
- DE flags and 300 bp promoter sequences.

External GEO/SRA/RefSeq content and unpinned packages can change over time. Exact equivalence claims require archiving the image digest and downloaded raw/reference inputs, not only the accession list.
