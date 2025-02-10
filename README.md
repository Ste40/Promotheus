# Promotheus: a bioinformatic pipeline for the analysis of bacterial transcriptomic data deposited on NCBI databases

### Description
Promotheus is a data-driven pipeline designed to retrieve, process, and integrate publicly available transcriptomic data from bacterial species and their strains. Its main goal is to support the selection of promoters that exhibit either stable (constitutive) or variable (inducible) expression profiles, facilitating biosensor design in synthetic biology applications.
Promotheus performs data retrieval from public repositories, raw data pre-processing, normalization, and batch effect correction. Gene IDs harmonization is performed through orthogroup identification, allowing a multi-strain and/or multi-species comparisons.

To learn more about it, read the manuscript at XXX.

### Pipeline Overview

![Tubemap_updated (2)](https://github.com/user-attachments/assets/5db7ef1b-80d8-44d8-9aa4-53ccadbff13f)


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
| `OPERONS_DETECTION` | FALSE   | Enable or disable operons detection.                                                                        | Logical    | No                         |
| `AMINO_DIFF`        | 3       | Defines the orthology matching in terms of protein length differences (amino acids) during orthogoup identification.                                       | Numeric    | No                         |
| `SIMILARITY`        | 30      | minimum percentage of amino acid sequence similarity required to group genes into the same orthogroup for ID harmonization                                                          | Numeric    | No                         |
| `ASMcode`           | *NULL* | The Assembly accession code.                                                                 | Character | **Yes** (required)         |
| `OPERONS_FILE`      | *NULL* | Path to the operons file; **only required if** `OPERONS_DETECTION` is set to `FALSE`.                       | Character | **Yes** if `OPERONS_DETECTION=FALSE` |

### Important Requirements

- `ASMcode` **must** be specified. If not provided, the script will terminate.
- If `OPERONS_DETECTION` is set to `FALSE`, you **must** provide an `OPERONS_FILE`; otherwise, the script will terminate.

---

## Usage

The pipeline runs within a Docker container. You can set environment variables using the `-e` flag. Make sure to mount your input/output directories with the `-v` flag.

### Docker Run Example

Below is an example command to run the PROMotheus Pipeline. Adjust the paths and parameters as needed:

```bash
docker run -it \
  -e THREADS=4 \
  -e MEMORY=4000 \
  -e CORRECTION_DE=BH \
  -e PVALUE=0.01 \
  -e OPERONS_DETECTION=TRUE \
  -e AMINO_DIFF=5 \
  -e SIMILARITY=30 \
  -e ASMcode=ASM317683v1 \
  -e OPERONS_FILE=/path/to/operons_file \
  -v /path/to/input:/input \
  -v /path/to/output:/output \
  ste40/promotheus
```

