#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"

release_num=40

echo "Gene_name,Transcript_name,Ensembl_gene_id,Ensembl_transcript_id" > ${BASE_PATH}gencode.v${release_num}.pc_transcripts.nopary.cdsonly.name.id.csv
grep ">" ${BASE_PATH}gencode.v${release_num}.pc_transcripts.nopary.cdsonly.fa | cut -c 2- | awk -F"|" 'BEGIN{OFS="."}{print $6,$5,$2,$1}' | awk -F"." 'BEGIN{OFS=","}{print $1,$2,$3,$5}' >> ${BASE_PATH}gencode.v${release_num}.pc_transcripts.nopary.cdsonly.name.id.csv
