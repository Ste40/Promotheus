#!/bin/bash

echo "
##################################
# Welcome to Promotheus pipeline #
##################################

Promotheus is a software for the analysis of gene expression data deposited on public databases.
Below are the available parameters:

- THREADS: Set the number of threads for parallel processing (default: 8, numeric).
- MEMORY: Set the maximum memory allocation for rockhopper execution (default: 4000, numeric).
- CORRECTION_DE: Specify the correction method for multiple testing. To see the available options, see p.adjust {stats} R documentation (default: none, character).
- PVALUE: Set the significance threshold for differentially expressed genes identification (default: 0.05, numeric).
- OPERONS_DETECTION: Enable operons detection (default: TRUE, logical).
- AMINO_DIFF: Set the number of amino acids by which proteins can differ in length (default: 3, numeric).
- SIMILARITY: Set the minimum similarity percentage between proteins (default: 30%, numeric).
- ASMcode: Specify the ASMcode for orthofinder mapping (**required**, character).
- OPERONS_THRESHOLD: Set the threshold for operon detection (default: 3, numeric).
- MULTITAXA: Enable indipendet specie batch effect correction (default: FALSE, logical).

Docker must be run in interactive mode, specified by -it.
Example Usage:
docker run -it \
  -e THREADS=4 \
  -e MEMORY=4000 \
  -e CORRECTION_DE=BH \
  -e PVALUE=0.01 \
  -e OPERONS_DETECTION=TRUE \
  -e AMINO_DIFF=5 \
  -e SIMILARITY=30 \
  -e ASMcode=ASM317683v1 \
  -e OPERONS_THRESHOLD=3 \
  -e MULTITAXA=TRUE \
  -v /path/to/input:/input \
  -v /path/to/output:/output \
  ste40/promotheus

NOTE: Adjust the volume paths for input and output according to your local file system.
"

# Imposta i valori di default se le variabili d'ambiente non sono già definite
: ${THREADS:=8}
: ${MEMORY:=4000}
: ${CORRECTION_DE:="none"}
: ${PVALUE:=0.05}
: ${OPERONS_DETECTION:=FALSE}
: ${AMINO_DIFF:=3}
: ${SIMILARITY:=30}
: ${ASMcode:=""}
: ${OPERONS_THRESHOLD:=3}
: ${MULTITAXA:=FALSE}

# Verifica che ASMcode sia stato specificato
if [ -z "$ASMcode" ]; then
  echo "Error: ASMcode is a required parameter. Please specify it via the ASMcode environment variable." >&2
  exit 1
fi

echo "Running the pipeline with the selected parameters:"
echo "THREADS: $THREADS"
echo "MEMORY: $MEMORY MB"
echo "CORRECTION_DE: $CORRECTION_DE"
echo "PVALUE: $PVALUE"
echo "OPERONS_DETECTION: $OPERONS_DETECTION"
echo "AMINO_DIFF: $AMINO_DIFF amino acids"
echo "SIMILARITY: $SIMILARITY%"
echo "ASMcode: $ASMcode"
echo "OPERONS_THRESHOLD: $OPERONS_THRESHOLD"
echo "MULTITAXA: $MULTITAXA"

# Esegue lo script R per l'analisi microarray
Rscript /scripts/Microarray/main_pipeline.R
