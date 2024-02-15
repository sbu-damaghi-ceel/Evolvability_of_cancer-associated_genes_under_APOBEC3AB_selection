#! /bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ceel/miniconda3/envs/cdur/lib
export CFLAGS="-I/home/ceel/miniconda3/envs/cdur/include"
export LDFLAGS="-L/home/ceel/miniconda3/envs/cdur/lib"

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"
CDUR_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/programs/CDUR"

mkdir ${BASE_PATH}cdur_celegans_refseq_run1 ${BASE_PATH}cdur_celegans_refseq_run2 ${BASE_PATH}cdur_celegans_refseq_run3

for i in {1..3}
do
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part0.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part0.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part1.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part1.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part2.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part2.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part3.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part3.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part4.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part4.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part5.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part5.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part6.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part6.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}caenorhabditis_elegans_pc_transcripts_part7.fa -m ./motif.txt -o ${BASE_PATH}cdur_celegans_refseq_run${i}/ -d 1 &> refseq_run${i}_part7.log &
done

