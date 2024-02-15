#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

FILESDIR=/mnt/c/Users/CEEL-PC-005/Desktop/Joon/CDUR_trajectory_research/ortholog_analysis/files/

datasets summary gene taxon 'pteropus alecto' > ${BASE_PATH}all_palecto.json
python3 ./get_short_info_from_summary.py

cat ${BASE_PATH}ncbi_palecto.tsv | sed "s/', '/,/g" | sed "s/\['//g" | sed "s/']//g" > ${BASE_PATH}ncbi_palecto_clean.tsv
