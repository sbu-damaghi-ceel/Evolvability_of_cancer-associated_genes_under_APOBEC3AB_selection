#! /bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ceel/miniconda3/envs/cdur/lib # CHANGE HERE
export CFLAGS="-I/home/ceel/miniconda3/envs/cdur/include" # CHANGE HERE
export LDFLAGS="-L/home/ceel/miniconda3/envs/cdur/lib" # CHANGE HERE

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
CDUR_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/programs/CDUR" # CHANGE HERE

mkdir ${BASE_PATH}cdur_refseq_run1 ${BASE_PATH}cdur_refseq_run2 ${BASE_PATH}cdur_refseq_run3

nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}gencode.v40.pc_transcripts.nopary.cdsonly.refseq.fa -o ${BASE_PATH}cdur_refseq_run1/ -d 1 &> refseq_run1.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}gencode.v40.pc_transcripts.nopary.cdsonly.refseq.fa -o ${BASE_PATH}cdur_refseq_run2/ -d 1 &> refseq_run2.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}gencode.v40.pc_transcripts.nopary.cdsonly.refseq.fa -o ${BASE_PATH}cdur_refseq_run3/ -d 1 &> refseq_run3.log &

