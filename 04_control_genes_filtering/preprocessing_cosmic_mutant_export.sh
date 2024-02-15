#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

awk -F"\t" 'BEGIN{OFS="\t"}{if ($1 ~ "_") {split($1,a,"_"); print a[1],$2} else {print $1,$2}}' ${BASE_PATH}CosmicMutantExport.tsv | awk -F"\t" 'BEGIN{OFS="\t"}NR>1{print $0}' | sort -u > ${BASE_PATH}CosmicMutantExport_genename_id_sorted_uniq.tsv