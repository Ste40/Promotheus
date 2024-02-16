# SteBora: a bioinformatic pipeline for the analysis of bacteria transcriptomics data deposited on NCBI databases

### Description
SteBora is a bioinformatic pipeline for the analysis of bacterial RNA-Seq and Microarray data deposited on Sequence Read Archive (SRA) and Gene Expression Omnibus (GEO) datasets.
The pipeline implements an automated workflow for differential gene expression analysis of transcriptomics datasets of a given bacteria strain.

This project was designed to support promoter selection for synthetic circuits engineering.
To learn more about it, read the manuscript at XXX.

### Pipeline overview

QUI CI VA LA MAPPA

## Installation and usage


Una volta installato docker sul proprio PC, e' necessario fare un pull dell'immagine, ad esempio:

'docker pull ste40/stebora_150124_debugging' oppure 'docker pull ste40/stebora_150124' 

Una volta fatto il pull della pipeline, e' possibile eseguirla:

'docker run -it -v /path/locale/tabelle:/input -v /path/locale/output:/output ste40/stebora_150124_debugging /bin/bash' 

**NOTA**: per navigare all'interno dell'immagine da terminale, va specificato '/bin/bash'. Questo non e' necessario per eseguire l'immagine ste40/stebora_150124: 

'docker run -it -v /path/locale/tabelle:/input -v /path/locale/output:/output ste40/stebora_150124' 

## Parametri 

La pipeline accetta diversi paramtri in input, i requisiti minimi sono i due volumi (-v), ovvero i path alle risorse locali sul tuo pc:

- -v /path/locale/tabelle:/input : **/path/locale/tabelle** e' il path sul tuo PC dove sono presenti le due tabelle .xlsx relative ai GSE microarray e RNA_Seq, con i nomi Experimental_Groups_MA.xlsx e Experimental_Groups_RNA_Seq.xlsx. mentre **/input** e' la cartella input presente nell'ambiente dell'immagine, di fatto con -v si monta un volume.
- -v /path/locale/output:/output : In questo caso si monta una cartella locale di output che e' dove saranno salvati tutti i file prodotti dalla pipeline.

docker acetta anche questi parametri in input:

- \`THREADS\`: Set the number of threads for parallel processing (default: 8, numeric).
- \`MEMORY\`: Set the maximum memory allocation for rockhopper execution (default: 4000, numeric).
- \`CORRECTION_DE\`: Specify the correction method for multiple testing. To see the availabe options, see p.adjust {stats} R documentation (default: none, character).
- \`PVALUE\`: Set the significance threshold for differentially expressed genes identification (default: 0.05, numeric). 
- \`OPERONS_DETECTION\`: Enable operons detection (default: FALSE, logical). 

La pipeline usera' i valori di default se non specificati. Possono essere specificati con '-e':

docker run -it -e THREADS=4 -e MEMORY=4000 -e CORRECTION_DE=BH -e PVALUE=0.01 -e OPERONS_DETECTION=TRUE -v /path/to/input:/input -v /path/to/output:/output ste40/stebora_150124_debugging /bin/bash

**NOTA**: CORRECTION_DE e PVALUE definiscono i rispettivi valori utilizzati dal KW per l'identificatione dei DE, non ho modificato altri pvalues o eventuali correzioni in altri punti della pipeline per l'analisi microarray.
THREADS setta i threads dove possibile, ad esempio orthofinder, allineamento ecc. Se ci sono funzioni specifiche nel quale sono settabili nella tua pipeline, possono essere aggiunti.

## Modificare gli script

Per modificare gli script, puoi montare direttamente alla cartella script, gli script in locale, ad esempio: 

'docker run -it -v /path/locale/tabelle:/input -v /path/locale/output:/output -v /path/locale/scripts:/scripts ste40/stebora_150124_debugging /bin/bash'

In questo modo, puoi modificare gli script in locale e lanciare la pipeline con gli script aggiornati. Nel momento in cui hai fatto le modifiche del caso puoi anche creare una nuova immagine, ti ho messo tra i collaboratori della repository. Basta anche solo che mi passi gli script e posso crearla io.  

