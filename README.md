# Promotheus: a bioinformatic pipeline for the analysis of bacterial transcriptomic data deposited on NCBI databases

### Description
Promotheus is a data-driven pipeline designed to retrieve, process, and integrate publicly available transcriptomic data from bacterial species and their strains. Its main goal is to support the selection of promoters that exhibit either stable (constitutive) or variable (inducible) expression profiles, facilitating biosensor design in synthetic biology applications.
Promotheus performs data retrieval from public repositories, raw data pre-processing, normalization, and batch effect correction. Gene IDs harmonization is performed through orthogroup identification, allowing a multi-strain and/or multi-species comparisons.

This project was carried out with the contribution of Stefano Quaranta, Debora Dallera and Lorenzo Pasotti.

Promotheus is distributed under the GNU GPL v3. Third-party components retain their own licenses; see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md), including the source and license information for the bundled Rockhopper JAR.

### Pipeline Overview

<img width="1066" alt="Tubemap_updated_github" src="https://github.com/user-attachments/assets/42159e26-9320-475d-b7b6-ec4db3d7f6ce" />


---

## Installation and execution modes

Promotheus supports two execution modes:

| Mode | Recommended for | Requirements |
| --- | --- | --- |
| Docker | Running the container entrypoint directly | Docker and this repository |
| Nextflow | Reproducible runs, input/output validation, provenance, reports and restart support | Nextflow, Java 17+ and Docker |

Both modes use the Promotheus scientific functions in the pinned Docker image. `main.nf` now schedules orthology preparation, per-genome comparisons and merge, per-GSE Microarray, RNA-Seq download, RNA-Seq alignment, RNA-Seq analysis, integration and validation processes. The former single-container workflow is preserved as `main_legacy.nf`.

See [Modular Nextflow execution](docs/MODULAR_NEXTFLOW.md) for the DAG, resource limits, restart behavior and validation limits. Validation includes synthetic tests, a completed two-GSE Docker run, a two-comparison real orthology timing benchmark, and two uncached orthology runs with identical results across all 127 comparisons in the archived genome set. Biological equivalence with the legacy orchestration still requires a dedicated comparison.

Install [Docker](https://docs.docker.com/get-docker/) and ensure that its daemon is running. To use Nextflow, also install [Nextflow 23.10 or newer](https://www.nextflow.io/docs/latest/install.html) and Java 17 or newer.

### Pulling the Docker Image

Once Docker is installed, pull the exact Promotheus image used by the Nextflow workflow:

```bash
docker pull ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7
```

## Parameters


| Environment Variable | Default | Description                                                                                                 | Value | **Required?**             |
|----------------------|---------|-------------------------------------------------------------------------------------------------------------|------------|----------------------------|
| `THREADS`           | 8       | Number of threads for parallel processing.                                                                  | Numeric    | No                         |
| `MEMORY`            | 4000    | Maximum memory allocation (MB) for Rockhopper execution.                                            | Numeric    | No                         |
| `CORRECTION_DE`     | none    | Correction method for multiple testing. See [p.adjust in R documentation](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust) for available options. | Character | No                         |
| `PVALUE`            | 0.05    | Significance threshold for identifying differentially expressed genes.                                      | Numeric    | No                         |
| `OPERONS_DETECTION` | TRUE   | Enable or disable operons detection.                                                                        | Logical    | No                         |
| `AMINO_DIFF`        | 3       | Defines the orthology matching in terms of protein length differences (amino acids) during orthogroup identification.                                       | Numeric    | No                         |
| `SIMILARITY`        | 30      | Minimum percentage of amino acid sequence similarity required to group genes into the same orthogroup for ID harmonization                                                          | Numeric    | No                         |
| `ASMcode`           | *NULL* | Assembly accession code for orthofinder mapping.                                                                 | Character | **Yes** (required)         |
| `OPERONS_THRESHOLD`      | 3 | Number of GSEs per microorganism used for operon identification. | Numeric | No |
| `MULTITAXA`      | FALSE | Enable or disable taxon-specific normalization with ComBat. | Logical | No |
| `MAX_GSE_PARALLEL` | 2 | Maximum independent GSE/organism workers. Total `THREADS` and `MEMORY` are divided between workers. | Positive integer | No |
| `PROMOTHEUS_KEEP_CACHE` | TRUE | Preserve reusable Microarray and orthology intermediates. | Logical | No |

### Important notes

- `ASMcode` **must** be specified. If not provided, the script will terminate. This genome is **not** used for reads mapping, but rather for orthology-based mapping of GeneIDs across species, reported in the final output table.
- If `OPERONS_THRESHOLD` is greater than the number of GSEs available for a given species, it will be adjusted to match the number of available GSEs.
- When `MULTITAXA` = TRUE, ComBat normalization is performed separately for each taxon, assigning the GSEs of each taxon to distinct batches.
- In modular Nextflow, RNA-Seq download and alignment run per GSE, including within the same organism. With operons disabled, subsequent analysis is also per GSE. With operons enabled, downstream analysis remains ordered within each organism to preserve operon-consensus behavior. Independent Microarray GSEs and the RNA-Seq branch can overlap.

## Input files

The pipeline needs two separate input tables for RNA-Seq and Microarray datasets. Each table must contain the following mandatory columns:

- **Accession**: Accession codes for various datasets (GSEs).
- **taxon**: Microorganism names, based on NCBI taxonomy name.
- **reference_genome**: Assembly accession codes for the reference genome used for reads alignment and gene counts quantification.

The examples of these two tables (RNA-Seq and Microarray) are provided in this repository, see Examples directory.

NOTE: It is recommended to use our GEO_Search.R script to generate these two tables.
See below for more information.

### Generating input files with GEO_Search.R R script

We provide GEO_Search.R script to help user to find **RNA-Seq** and **Microarray** experiments on NCBI GEO Datasets based on a provided **taxonomic ID (taxID)**.

Features:
- **Taxonomy Retrieval**: Retrieves the scientific name and subspecies for a given **taxID**.
- **Genome Assembly**: Fetches the reference genome assembly version for the organism.
- **GEO Datasets Search**: Searches for RNA-Seq and Microarray datasets from GEO, filtering by dataset type and minimum sample count (4).
- **Export**: Outputs filtered results to **Excel** files for RNA-Seq and Microarray experiments.

R (version 3.6 or higher recommended) must be installed in your system.
- Required R packages:
  - `devtools`
  - `rentrez`
  - `openxlsx`
  - `optparse`

To install the required packages, run the following in R:

```r
install.packages(c("devtools", "openxlsx", "optparse"))
devtools::install_github("ropensci/rentrez")
```

Then, the script can be run from command-line, specifying the taxID of interest:

```bash
Rscript GEO_Search.R --taxID <taxID>
```
NOTE: The script attempts to retrieve data from NCBI, so it may fail if the servers are down or experiencing issues.

## Usage

The standard input consists of one Microarray workbook and one RNA-Seq workbook. Their filenames must contain `Microarray` and `RNA_Seq`, respectively. Nextflow requires exactly one workbook of each type and rejects unrelated `.xlsx` files, so keep previous results such as `Output.xlsx` outside the input directory.

### Option 1: Docker

The pipeline runs within a Docker container. Set parameters with the `-e` flag and mount the input and output directories with `-v`. Run the command from the repository root: the ten additional read-only mounts apply the same compatibility and performance fixes used automatically by Nextflow.

WARNING: the scripts provided in this repository can't be executed without the appropriate environment and softwares/packages, which are included in Docker's image.

### Docker Run Example

Below is an example command to run the Promotheus pipeline. Adjust the paths and parameters as needed:

```bash
docker run --rm -it \
  --platform linux/amd64 \
  -e THREADS=4 \
  -e MEMORY=4000 \
  -e CORRECTION_DE=BH \
  -e PVALUE=0.01 \
  -e OPERONS_DETECTION=TRUE \
  -e AMINO_DIFF=5 \
  -e MULTITAXA=TRUE \
  -e SIMILARITY=30 \
  -e ASMcode=ASM317683v1 \
  -e OPERONS_THRESHOLD=3 \
  -e MAX_GSE_PARALLEL=2 \
  -e PROMOTHEUS_KEEP_CACHE=TRUE \
  -v /path/to/input/:/input \
  -v /path/to/output/:/output \
  -v "$PWD/Scripts/Microarray/ortologhi_new.R:/scripts/Microarray/ortologhi_new.R:ro" \
  -v "$PWD/Scripts/Microarray/dowload_file_assembly.R:/scripts/Microarray/dowload_file_assembly.R:ro" \
  -v "$PWD/Scripts/Microarray/metadata_retrieval.R:/scripts/Microarray/metadata_retrieval.R:ro" \
  -v "$PWD/Scripts/Microarray/intra_experiment_analysis.R:/scripts/Microarray/intra_experiment_analysis.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/Data_retrieval.R:/scripts/RNA-Seq/Data_retrieval.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/BAM_parsing.R:/scripts/RNA-Seq/BAM_parsing.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/NDE_detection.R:/scripts/RNA-Seq/NDE_detection.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:/scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/Operons_detection.R:/scripts/RNA-Seq/Operons_detection.R:ro" \
  -v "$PWD/Scripts/Microarray/interexp_analysis.R:/scripts/Microarray/interexp_analysis.R:ro" \
  ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7
```

### Option 2: Nextflow (recommended)

Clone this repository and run the command from its root directory, where `main.nf` and `nextflow.config` are located. Docker Desktop or the Docker daemon must already be running.

For example, prepare a clean input directory from the supplied workbooks:

```bash
mkdir -p promotheus-input
cp Examples/GEO_Datasets_Search_Microarrays_Lactococcus_lactis_2025-01-15.xlsx promotheus-input/
cp Examples/GEO_Datasets_Search_RNA_Seq_Lactococcus_lactis_2025-01-15.xlsx promotheus-input/
```

Do not pass the complete `Examples` directory because it also contains `Output.xlsx`, which is an output workbook rather than a pipeline input.

Run the standard Docker-backed workflow:

```bash
nextflow run main.nf \
  -profile standard \
  --input_dir "$PWD/promotheus-input" \
  --output_dir "$PWD/promotheus-output" \
  --asmcode ASM317683v1 \
  --threads 8 \
  --memory 4000 \
  --correction_de none \
  --pvalue 0.05 \
  --operons_detection true \
  --amino_diff 3 \
  --similarity 30 \
  --operons_threshold 3 \
  --multitaxa false \
  --max_gse_parallel 2 \
  --keep_cache true \
  --export_tsv true
```

`--asmcode` is required and must use the `ASM<digits>v<digits>` format. All other parameters have the defaults shown above. Large biological runs can take several hours and download substantial RefSeq, GEO and SRA data.

To continue an interrupted run, repeat the same command and append:

```bash
-resume
```

Nextflow reuses completed tasks whose code, inputs and parameters have not changed. Do not delete the `work/` directory before resuming.

The modular workflow uses Nextflow task caching for individual GSE stages. `--threads` and `--memory` set the total local resource budget; `--task_cpus` and `--task_memory` set resources per task. `--max_gse_parallel` limits concurrency. Keep `work/` to reuse downloads and alignments with `-resume`.

The former content-addressed `.promotheus-cache` and `--cache_dir` options remain available through `main_legacy.nf`. Modular runs keep their intermediates in the task work directories; `--keep_cache false` does not remove those directories. See the modular guide for staging disk requirements and curated experimental-group limitations.

The pinned image is `linux/amd64`. On an ARM64 system it therefore runs under Docker emulation. A compatible custom multi-architecture image can be selected with `--container_image owner/image:tag --container_platform auto`; `auto` omits the platform override and lets Docker select the native manifest. It does not make the pinned AMD64-only image native.

The main results are written under `--output_dir`:

```text
promotheus-output/
├── orthologous_genes.txt
├── RNA_Seq/
├── results/
│   ├── fulltab/
│   ├── final_table.xlsx
│   └── final_table.tsv
├── pipeline_info/
│   ├── execution_report_YYYY-MM-DD_HH-MM-SS.html
│   ├── execution_timeline_YYYY-MM-DD_HH-MM-SS.html
│   ├── execution_trace_YYYY-MM-DD_HH-MM-SS.txt
│   └── pipeline_dag_YYYY-MM-DD_HH-MM-SS.html
└── .promotheus-nextflow/
    ├── input_manifest.tsv
    ├── input_validation.json
    ├── run_record.tsv
    ├── checkpoint_manifest.tsv
    ├── output_manifest.tsv
    └── output_validation.json
```

`output_validation.json` must report `"status": "pass"` for a successfully validated run. `checkpoint_manifest.tsv` records a reusable full-tab status and checksum for every requested accession. The RNA-Seq downloader retries each failed `fasterq-dump` operation up to three times (with 30- and 60-second backoff), removes incomplete FASTQ files between attempts, accepts only one non-empty single-end file or two non-empty paired-end files per run accession, and retries within each GSE task. FASTQ compression uses gzip's fast level and compresses paired files concurrently; BAM filtering, sorting and indexing use the task's assigned CPUs; gene statistics and percentiles use allocation-efficient matrix operations that are regression-tested against the legacy calculations. Before output validation, Nextflow repairs the known malformed drawing relationships and worksheet dimension emitted by the legacy container's `openxlsx` version, then exports the same table as TSV for efficient downstream reuse. The validator checks Excel and TSV consistency, ZIP integrity, internal OOXML relationships, content-type references, worksheet bounds and the expected biological result columns. It also compares every GSE accession requested in the validated input manifest with the expression and percentile columns in `final_table.xlsx`; a partial analysis therefore fails validation instead of being reported as successful.

### Execution reports

Each run writes nf-core-style execution artifacts under `--output_dir/pipeline_info`, including an interactive HTML resource report, an HTML timeline, a tab-separated trace and the workflow DAG. Filenames contain the launch timestamp, so a subsequent run does not hide the previous report. The trace records status, exit code, elapsed time, CPU usage, peak resident/virtual memory, I/O counters and work directory for every Nextflow process. The modular workflow reports each GSE download, alignment and Microarray task separately, together with RNA-Seq analysis tasks, orthology, integration and validation. In `main_legacy.nf`, scientific calculations remain aggregated under `RUN_PROMOTHEUS`.

### Validation status and known data limitation

The corrected full workflow completed again on 2 November 2026 after resuming
from SRA transfer failures: 136 tasks accounted for, 23 reused from cache.
Validation passed with 352 genes, 28 columns and all ten samples; the orthology
map matched the independently validated full snapshot, and comparison with the
earlier final table found no differing values within the stated numeric
tolerances. See [complete-run evidence and timings](docs/MODULAR_NEXTFLOW.md#complete-corrected-run--29-october-to-2-november-2026).

The modular workflow completed a real Docker test on 28 October 2026 with `GSE8182` (six Microarray samples) and `GSE116952` (four RNA-Seq samples), reference `ASM317683v1` and operons disabled. Output validation passed with **352 genes and 28 columns**, including all ten samples. Six completed tasks were reused on resume after correcting GFF attribute lookup in integration. See [the modular validation evidence and limitations](docs/MODULAR_NEXTFLOW.md#real-execution-on-2728-october-2026).

A real Docker-backed run of the supplied Lactococcus input snapshot was completed on 9 September 2026. `RUN_PROMOTHEUS` produced a structurally valid `final_table.xlsx` containing 339 genes and 136 columns. Six of the 18 requested GSE accessions were processed (`GSE116952`, `GSE119647`, `GSE26536`, `GSE51494`, `GSE74554`, and `GSE8182`). The legacy Microarray implementation skipped the other 12 because their current GEO supplementary files were absent or unreadable: `GSE10254`, `GSE10256`, `GSE12962`, `GSE20373`, `GSE209756`, `GSE229414`, `GSE23987`, `GSE2823`, `GSE43875`, `GSE4872`, `GSE5761`, and `GSE9436`.

Consequently, strict output validation correctly reports the supplied 18-accession run as incomplete rather than returning a false success. This is an upstream-data/legacy-parser limitation, not a silent Nextflow success. Users should review `output_validation.json`; a production run is considered fully validated only when its status is `pass`.

## Output file

The output table contains gene expression values, percentiles, coefficient of variation (CV), and differential expression analysis results.

### Columns Description
- **Gene_name**: The name of the gene.
- **GeneID**: The unique identifier for the gene, according to the gene nomenclature of the user-provided ASMcode.
- **Expression_GSEXXXXX**: Gene expression values for each GSE.
- **Mean_Expression**: The average expression value of the gene across all samples in the datasets.
- **CV**: Coefficient of variation, indicating the variability of the gene expression across samples, computed as the standard deviation of log-transformed expression values. A higher CV suggests more variation in gene expression.
- **Percentile_GSEXXXXX**: Percentile values for the gene expression in different datasets.
- **Mean_Percentile**: The average percentile value across all datasets for the gene.
- **Percentile_SD**: Standard deviation of the percentile values across all datasets for the gene.
- **DE**: Differential expression status (YES|NO), in at least one GSE.
- **Sequence**: Contains the DNA sequence upstream of the gene.

An example of output table is provided in this repository (Output.xlsx)

---

## Legacy internal execution flow (Docker entrypoint and main_legacy.nf)

Promotheus has a single entrypoint (`Scripts/Introduction.sh`) that validates environment variables and starts `Scripts/Microarray/main_pipeline.R`. The R main script orchestrates the full workflow as follows:

1. **Parameter parsing/validation** from environment variables (`THREADS`, `MEMORY`, `CORRECTION_DE`, `PVALUE`, `OPERONS_DETECTION`, `OPERONS_THRESHOLD`, `AMINO_DIFF`, `SIMILARITY`, `MULTITAXA`, `ASMcode`).
2. **Directory setup** under `/output`:
   - `/output/tmp/microarray`
   - `/output/tmp/orthologous_genome`
   - `/output/results/fulltab`
   - `/output/RNA_Seq/results`
3. **Orthology step** (if `/output/orthologous_genes.txt` does not exist):
   - `Scripts/Microarray/ortologhi_new.R`
   - `Scripts/Microarray/dowload_file_assembly.R`
4. **Microarray branch** (if input contains files matching `Microarray`):
   - `Scripts/Microarray/pipeline_MA_bash.R`
   - metadata retrieval + per-experiment analyses (`intra_experiment_analysis.R`, one/two-channel functions, affymetrix).
5. **RNA-Seq branch** (if input contains files matching `RNA_Seq`):
   - `Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R`
   - metadata retrieval, download/alignment (`Data_retrieval.R`, `BAM_parsing.R`), optional operon detection, NDE detection, merge.
6. **Cross-experiment integration**:
   - `Scripts/Microarray/interexp_analysis.R`
   - intersects orthology-mapped genes across Microarray and RNA-Seq full tables, applies normalization/batch correction, computes summary metrics, writes final outputs.

### Data flow between folders

- **Input**: `/input/*Microarray*` and/or `/input/*RNA_Seq*` (xlsx tables).
- **Intermediate**:
  - `/output/tmp/microarray/*_full_tab.txt` (microarray per-GSE tables)
  - `/output/RNA_Seq/results/*full_tab.txt` (RNA-Seq per-GSE tables)
  - `/output/orthologous_genes.txt` (cross-genome gene mapping)
- **Final**:
  - `/output/results/fulltab/` (copied per-GSE full tables)
  - `/output/results/` final integrated result table(s), with expression/CV/percentiles/DE and promoter sequence annotation.

---

## Nextflow implementation and validation

The DSL2 workflow (`main.nf` + `nextflow.config`) runs the pinned `ste40/promotheus` image through isolated stage adapters. Scripts and inputs are mounted read-only; every task owns its output directory. The existing R scientific functions are reused.

### Workflow stages

- `VALIDATE_INPUTS`: validates both workbooks, records input checksums and creates a deterministic input shard per GSE.
- `ORTHOLOGY_PREPARE`: downloads/parses genomes and prepares the preferred reference and proteomes.
- `ORTHOLOGY_PAIR`: compares the preferred assembly with one other genome per task, with independent caching and restart.
- `ORTHOLOGY_MERGE`: restores the original genome order and publishes the cross-genome gene map.
- `MICROARRAY_GSE`: performs Microarray analysis in one task per GSE.
- `RNA_DOWNLOAD_GSE`: retrieves metadata, reads and reference files per GSE.
- `RNA_ALIGN_GSE`: aligns reads and prepares BAM files per GSE.
- `RNA_ANALYZE`: analyzes individual GSEs with operons disabled, or ordered groups of GSEs from the same organism with operons enabled.
- `INTEGRATE_RESULTS`: waits for all tables, integrates experiments, adds promoters, repairs Excel output and optionally exports TSV.
- `VERIFY_OUTPUTS`: checks workbook integrity, schema and coverage of every requested GSE.

Orthology and RNA-Seq downloads can overlap; Microarray and RNA-Seq tasks can also overlap. The [modular guide](docs/MODULAR_NEXTFLOW.md) explains the remaining operon barrier and resource allocation.

Orthology concurrency is controlled separately with `--max_orthology_parallel 4
--orthology_cpus 1 --orthology_memory 3000`; allocate enough total resources,
for example `--threads 8 --memory 16000`. With the default 4000 MB total budget,
only one 3000 MB pair fits. Add `-entry ORTHOLOGY_ONLY` to run just this branch.
On two real comparisons, two simultaneous tasks reduced elapsed time from
8 min 32 s to 4 min 59 s (1.71×). The investigation of variability on sequential
reruns found unstable matrix and candidate-gene ordering, now stabilized by the
modular adapter. See the [repeatability correction](docs/ORTHOLOGY_REPRODUCIBILITY.md)
and [timing benchmark](docs/MODULAR_NEXTFLOW.md#parallel-orthology-benchmark-on-28-october-2026).
These checks do not establish biological equivalence with legacy results.

For complete runs, opt into `-profile standard,balanced` to reserve CPU and
memory for one RNA-Seq task. For example, with `--threads 8 --memory 20000
--task_cpus 2 --task_memory 8000 --orthology_cpus 1 --orthology_memory 3000
--max_orthology_parallel 6`, the effective orthology limit is **4**. The limit is
the minimum of the requested concurrency and the CPU/memory slots remaining
after reserving one scientific task. An extra executor queue slot allows that
task to be submitted alongside the orthology jobs. A budget that cannot fit
both tasks fails before execution.

The reserve remains in effect after RNA-Seq finishes, so it can lengthen the
orthology branch; a shorter overall runtime is not guaranteed. Default runs
and `ORTHOLOGY_ONLY` should use `standard` without `balanced`. Changing this
scheduling option preserves reusable scientific results with `-resume`.

The previous workflow is retained as `main_legacy.nf`, including its cache and curated-group behavior. The September real-data evidence applies to that legacy orchestration. The modular adapter has separate successful two-GSE real-run evidence from October, described above; formal biological equivalence remains untested.

### Deterministic test profile

The `test` profile replaces Docker with a deterministic local runner. It tests input validation, parameter propagation, dataflow, publication and output verification without downloading NCBI data:

```bash
nextflow run main.nf \
  -profile test \
  --input_dir tests/fixtures/toy_input \
  --output_dir /tmp/promotheus-test-output \
  --asmcode ASM100001v1
```

The toy input workbooks mirror the 29-column Promotheus input contract and contain one synthetic Microarray GSE and one synthetic RNA-Seq GSE. The expected result contains four genes and formula-checked expression, variability, percentile, DE and 300 bp promoter fields. This profile validates orchestration only; synthetic accession codes are intentionally not submitted to GEO/SRA and therefore do not exercise the biological analysis inside the Docker image.

### Test suite

Run all local unit and integration tests with:

```bash
bash tests/shell/test_all.sh
```

With Docker running and the pinned image already pulled, verify the real container, the compatibility overlays, and the identity of unchanged critical scripts with:

```bash
bash tests/shell/test_container.sh
```

See `docs/NEXTFLOW_EQUIVALENCE.md` for the compatibility argument, test coverage, and limitations.
