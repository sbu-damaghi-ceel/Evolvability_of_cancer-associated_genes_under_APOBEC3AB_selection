#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/10_additional_work_1/" # CHANGE HERE

FILESDIR=/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/10_additional_work_1/

datasets summary gene taxon 'pan troglodytes' > ${BASE_PATH}all_ptroglodytes.json
python3 ./get_short_info_from_summary.py

cat ${BASE_PATH}ncbi_ptroglodytes.tsv | sed "s/', '/,/g" | sed "s/\['//g" | sed "s/']//g" > ${BASE_PATH}ncbi_ptroglodytes_clean.tsv
