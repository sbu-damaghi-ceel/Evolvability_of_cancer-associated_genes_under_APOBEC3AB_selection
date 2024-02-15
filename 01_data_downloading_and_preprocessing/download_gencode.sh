#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

release_num=40

# download and unzip data
wget -P ${BASE_PATH} http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${release_num}/gencode.v${release_num}.pc_transcripts.fa.gz
gunzip -c ${BASE_PATH}gencode.v${release_num}.pc_transcripts.fa.gz > ${BASE_PATH}gencode.v${release_num}.pc_transcripts.fa
