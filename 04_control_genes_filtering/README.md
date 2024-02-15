# Generate control genes group

## Download all human genes from NCBI + preprocessing
We used [NCBI Datasets command line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) (CLI v13.x API v1) to download a list of human genes in NCBI. The shell script also run Python script `get_short_info_from_summary.py` to preprocess `.json` file and generate `.tsv` file. 
```console
./download_ncbi_human.sh
```

## Download point mutation data from COSMIC
We downloaded COSMIC Mutation Data from [COSMIC Data Downloads](https://cancer.sanger.ac.uk/cosmic/archive-download) page. One can download directly from the page or use script to download the file.
```console
./download_cosmic_mutant_export.sh
```

## Preprocess COSMIC Mutation Data
We preprecoess COSMIC Mutation Data to extract only gene name and ID.
```console
./preprocessing_cosmic_mutant_export.sh
``` 

## Filter control genes using biomaRt
We used R script `get_control_genes.r` which utilizes [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) to retreive both NCBI RefSeq ID and Ensembl ID so that we can filter out genes in COSMIC Mutation Data from the list of human genes in NCBI.

## Finalize control genes list
Below shell script utilizes `awk` to extract only protein coding genes and remove mitocondrial genes and genes starting with LOC. The final file `ncbi_human_clean_no_mutation_protein_coding_no_LOC_no_MT.tsv` contains control genes.
```
./finalize_control_genes.sh
``` 