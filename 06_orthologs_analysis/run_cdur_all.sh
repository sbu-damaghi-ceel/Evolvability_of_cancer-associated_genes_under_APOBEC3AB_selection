#! /bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ceel/miniconda3/envs/cdur/lib
export CFLAGS="-I/home/ceel/miniconda3/envs/cdur/include"
export LDFLAGS="-L/home/ceel/miniconda3/envs/cdur/lib"

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"
CDUR_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/programs/CDUR"

mkdir ${BASE_PATH}all_cancer_genes_orthologs_run2 ${BASE_PATH}all_trajs_all_cancer_genes_human_run2
mkdir ${BASE_PATH}all_control_genes_orthologs_run2 ${BASE_PATH}all_trajs_all_control_genes_human_run2

nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}all_cancer_genes_orthologs.fa -o ${BASE_PATH}all_cancer_genes_orthologs_run2/ -d 1 &> run2_log1.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}all_control_genes_orthologs.fa -o ${BASE_PATH}all_control_genes_orthologs_run2/ -d 1 &> run2_log2.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}all_trajs_all_cancer_genes_human.fa -o ${BASE_PATH}all_trajs_all_cancer_genes_human_run2/ -d 1 &> run2_log3.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}all_trajs_all_control_genes_human.fa -o ${BASE_PATH}all_trajs_all_control_genes_human_run2/ -d 1 &> run2_log4.log &


