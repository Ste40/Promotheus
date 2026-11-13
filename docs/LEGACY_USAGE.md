# Input preparation and legacy Docker usage

For the recommended Nextflow entrypoint, see the [README](../README.md).
The environment variables and Docker command below apply to the original entrypoint.

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


### Option 1: Docker

The pipeline runs within a Docker container. Set parameters with the `-e` flag and mount the input and output directories with `-v`. Run the command from the repository root: the additional read-only mounts apply the same compatibility and performance fixes used automatically by Nextflow.

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
  -v "$PWD/Scripts/RNA-Seq/MetaData_retrieval.R:/scripts/RNA-Seq/MetaData_retrieval.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/Groups_retrieval.R:/scripts/RNA-Seq/Groups_retrieval.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/BAM_parsing.R:/scripts/RNA-Seq/BAM_parsing.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/NDE_detection.R:/scripts/RNA-Seq/NDE_detection.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:/scripts/RNA-Seq/Pipeline_bash_RNA_Seq.R:ro" \
  -v "$PWD/Scripts/RNA-Seq/Operons_detection.R:/scripts/RNA-Seq/Operons_detection.R:ro" \
  -v "$PWD/Scripts/Microarray/interexp_analysis.R:/scripts/Microarray/interexp_analysis.R:ro" \
  ste40/promotheus@sha256:a7c5196b4ebf6f93ef8a4ff2698481d24e3f63306c9eaccd413bce3c48eebec7
```


## Output file

The output table contains gene expression values, percentiles, coefficient of variation (CV), and differential expression analysis results.

### Columns Description
- **Gene_name**: The name of the gene.
- **GeneID**: The unique identifier for the gene, according to the gene nomenclature of the user-provided ASMcode.
- **Expression_GSEXXXXX**: Gene expression values for each GSE.
- **Mean_Expression**: The average expression value of the gene across all samples in the datasets.
- **CV**: Sample standard deviation divided by the mean of the integrated expression values (`sd / mean`), reported as a ratio, not a percentage. At least two finite observations and a positive mean are required; otherwise analysis fails explicitly. This describes variability on the integrated expression scale (which may be log-transformed), not a CV reconstructed on the original abundance scale. Earlier versions incorrectly wrote the standard deviation into this column.
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

