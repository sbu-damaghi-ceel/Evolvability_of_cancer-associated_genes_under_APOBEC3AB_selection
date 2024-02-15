#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

awk -F"\t" '{if (NR>1) {if ($12 ~ "PROTEIN_CODING") print $0} else {print $0}}' ${BASE_PATH}ncbi_human_clean_no_mutation.tsv > ${BASE_PATH}ncbi_human_clean_no_mutation_protein_coding.tsv 
awk -F"\t" '{if (NR>1) {if ($4 !~ "LOC") print $0} else {print $0}}' ${BASE_PATH}ncbi_human_clean_no_mutation_protein_coding.tsv > ${BASE_PATH}ncbi_human_clean_no_mutation_protein_coding_no_LOC.tsv
awk -F"\t" '{if (NR>1) {if ($4 !~ "MT-") print $0} else {print $0}}' ${BASE_PATH}ncbi_human_clean_no_mutation_protein_coding_no_LOC.tsv > ${BASE_PATH}ncbi_human_clean_no_mutation_protein_coding_no_LOC_no_MT.tsv