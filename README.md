# Promotheus: a bioinformatic pipeline for the analysis of bacterial transcriptomic data deposited on NCBI databases

### Description
Promotheus is a data-driven pipeline designed to retrieve, process, and integrate publicly available transcriptomic data from bacterial species and their strains. Its main goal is to support the selection of promoters that exhibit either stable (constitutive) or variable (inducible) expression profiles, facilitating biosensor design in synthetic biology applications.
Promotheus performs data retrieval from public repositories, raw data pre-processing, normalization, and batch effect correction. Gene IDs harmonization is performed through orthogroup identification, allowing a multi-strain and/or multi-species comparisons.

To learn more about it, read the manuscript at XXX.

### Pipeline Overview

<img width="1066" alt="Tubemap_updated_github" src="https://github.com/user-attachments/assets/42159e26-9320-475d-b7b6-ec4db3d7f6ce" />


---

## Installation and Usage

Promotheus is deployed via Docker. If you do not already have Docker installed, follow the instructions for docker installation at https://docs.docker.com/

### Pulling the Docker Image

Once Docker is installed, you can pull the Promotheus image from Docker Hub by running:

```bash
docker pull ste40/promotheus
```
## Parameters


| Environment Variable | Default | Description                                                                                                 | Value | **Required?**             |
|----------------------|---------|-------------------------------------------------------------------------------------------------------------|------------|----------------------------|
| `THREADS`           | 8       | Number of threads for parallel processing.                                                                  | Numeric    | No                         |
| `MEMORY`            | 4000    | Maximum memory allocation (MB) for Rockhopper execution.                                            | Numeric    | No                         |
| `CORRECTION_DE`     | none    | Correction method for multiple testing. See [p.adjust in R documentation](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust) for available options. | Character | No                         |
| `PVALUE`            | 0.05    | Significance threshold for identifying differentially expressed genes.                                      | Numeric    | No                         |
| `OPERONS_DETECTION` | TRUE   | Enable or disable operons detection.                                                                        | Logical    | No                         |
| `AMINO_DIFF`        | 3       | Defines the orthology matching in terms of protein length differences (amino acids) during orthogoup identification.                                       | Numeric    | No                         |
| `SIMILARITY`        | 30      | Minimum percentage of amino acid sequence similarity required to group genes into the same orthogroup for ID harmonization                                                          | Numeric    | No                         |
| `ASMcode`           | *NULL* | Assembly accession code for orthofinder mapping.                                                                 | Character | **Yes** (required)         |
| `OPERONS_THRESHOLD`      | 3 | Number of GSEs per microorganism used for operon identification. | Character | No |
| `MULTITAXA`      | FALSE | Enable or disable taxon-specific normalization with ComBat. | Character | No |

### Important notes

- `ASMcode` **must** be specified. If not provided, the script will terminate. This genome is **not** used for reads mapping, but rather for orthology-based mapping of GeneIDs across species, reported in the final output table.
- If `OPERONS_THRESHOLD` is greater than the number of GSEs available for a given specie, it will be adjusted to match the number of available GSEs.
- When `MULTITAXA` = TRUE, ComBat normalization is performed separately for each taxon, assigning the GSEs of each taxon to distinct batches.

## Input files

The pipeline needs two separate input tables for RNA-Seq and Microarray datasets. Each table must contain the following mandatory columns:

- **Accession**: Accession codes for various datasets (GSEs).
- **taxon**: Microorganism names, based on NCBI taxonomy name.
- **reference_genome**: Assembly accession codes for the reference genome used for reads alignment and gene counts quantification.

The examples of these two tables (RNA-Seq and Microarray) are provided in this repository, see GEO_Datasets_Search_RNA_Seq_Lactococcus_lactis_2025-01-15.xlsx and GEO_Datasets_Search_Microarrays_Lactococcus_lactis_2025-01-15.xlsx.

NOTE: It is reccomended to use our GEO_Search.R script to generate these two tables.
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

The pipeline runs within a Docker container. You can set environment variables using the `-e` flag. Make sure to mount your input/output directories with the `-v` flag.

### Docker Run Example

Below is an example command to run the Promotheus pipeline. Adjust the paths and parameters as needed:

```bash
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
  -v /mnt/f/Umberto/BMS_Docker/input/:/input \
  -v /mnt/f/Umberto/BMS_Docker/output/:/output \
  ste40/promotheus
```
## Output file

The output table contains includes gene expression values, percentiles, coefficient of variation (CV), and Differential expression analysis results. 

### Columns Description
- **Gene_name**: The name of the gene.
- **GeneID**: The unique identifier for the gene.
- **Espressione_GSEXXXXX**: Gene expression values for each GSE.
- **Percentile_GSEXXXXX**: Percentile values for the gene expression in different datasets.
- **Espressione_media**: The average expression value of the gene across all samples in the datasets.
- **CV**: Coefficient of variation, indicating the variability of the gene expression across samples. A higher CV suggests more variation in gene expression.
- **Percentile_medio**: The average percentile value across all datasets for the gene.
- **Percentile_SD**: Standard deviation of the percentile values across all datasets for the gene.
- **DE**: Differential expression status (YES|NO)
- **sequence**: Contains the upstream gene nucleotide sequence.

An example of output table is provided in this repository (Output.xlsx)
