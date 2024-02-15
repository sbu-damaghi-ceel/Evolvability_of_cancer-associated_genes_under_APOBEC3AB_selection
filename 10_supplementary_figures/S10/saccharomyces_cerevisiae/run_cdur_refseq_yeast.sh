#! /bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jhsong/miniconda3/envs/cdur/lib
export CFLAGS="-I/home/jhsong/miniconda3/envs/cdur/include"
export LDFLAGS="-L/home/jhsong/miniconda3/envs/cdur/lib"

BASE_PATH="/home/jhsong/projects/yeast/"
CDUR_PATH="/home/jhsong/programs/CDUR"

mkdir ${BASE_PATH}cdur_yeast_run1 ${BASE_PATH}cdur_yeast_run2 ${BASE_PATH}cdur_yeast_run3

nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}saccharomyces_cerevisiae_pc_transcripts.fa -m ./motif.txt -o ${BASE_PATH}cdur_yeast_run1/ -d 1 &> refseq_run1.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}saccharomyces_cerevisiae_pc_transcripts.fa -m ./motif.txt -o ${BASE_PATH}cdur_yeast_run2/ -d 1 &> refseq_run2.log &
nohup python3 $CDUR_PATH/CDUR.py -i ${BASE_PATH}saccharomyces_cerevisiae_pc_transcripts.fa -m ./motif.txt -o ${BASE_PATH}cdur_yeast_run3/ -d 1 &> refseq_run3.log &

