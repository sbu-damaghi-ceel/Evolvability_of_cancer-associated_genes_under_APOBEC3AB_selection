#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

# Download cosmic caner gene census 
AUTHSTRING=$(echo "YOUR_MAIL:YOUR_PASSWORD" | base64) # CHANGE HERE
downloadlink=$( curl -s -H "Authorization: Basic ${AUTHSTRING}" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/cancer_gene_census.csv | awk -F":" '{print $2":"$3}' | sed 's/}//g' | sed 's/\"//g' )
echo "Downloading cancer_gene_census.csv"
curl -o ${BASE_PATH}cancer_gene_census.csv ${downloadlink}