#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
INDIR=${BASE_PATH}cancer_genes_orthologs/

if [ ! -d ${INDIR} ]
then
    mkdir ${INDIR}
fi

GENE_ID_FILE=${BASE_PATH}"cancer_gene_census_id_for_orthologs.txt"

n=$(expr $(cat ${GENE_ID_FILE} | wc -l) - 2)
genes=($(awk -F"," 'NR>1{print $1}' ${GENE_ID_FILE}))
ids=($(awk -F"," 'NR>1{print $2}' ${GENE_ID_FILE}))

for ((i=0; i<=${n}; i++))
do 
echo ${genes[i]} start...

GENENAME=${genes[i]}
GENEID=${ids[i]}

datasets download ortholog gene-id ${GENEID} --filename ${INDIR}${GENENAME}.zip

unzip ${INDIR}${GENENAME}.zip -d ${INDIR}${GENENAME}
rm ${INDIR}${GENENAME}.zip

dataformat tsv gene \
--inputfile ${INDIR}${GENENAME}/ncbi_dataset/data/data_report.jsonl \
--fields tax-id,tax-name,transcript-accession,transcript-length,transcript-protein-accession,transcript-protein-length,transcript-protein-isoform > ${INDIR}${GENENAME}/ncbi_dataset/data/transcript_protein.tsv

dataformat tsv gene \
--inputfile ${INDIR}${GENENAME}/ncbi_dataset/data/data_report.jsonl \
--fields tax-id \
--elide-header > ${INDIR}${GENENAME}/ncbi_dataset/data/taxids.txt

echo ${genes[i]} done...
done
