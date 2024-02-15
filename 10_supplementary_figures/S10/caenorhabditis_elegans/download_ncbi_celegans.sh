#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

datasets summary gene taxon 'caenorhabditis elegans' > ${BASE_PATH}all_celegans.json
python3 ./get_short_info_from_summary.py

cat ${BASE_PATH}ncbi_celegans.tsv | sed "s/', '/,/g" | sed "s/\['//g" | sed "s/']//g" > ${BASE_PATH}ncbi_celegans_clean.tsv
