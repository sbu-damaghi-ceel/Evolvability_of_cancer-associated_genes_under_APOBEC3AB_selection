#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/12_additional_work_3/" # CHANGE HERE

FILESDIR=/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/12_additional_work_3/

#datasets summary gene taxon 'Saccharomyces cerevisiae S288C' > ${BASE_PATH}all_scerevisiae.json
#datasets summary gene taxon 'xenopus laevis' > ${BASE_PATH}all_xlaevis.json



python3 ./get_short_info_from_summary.py

cat ${BASE_PATH}ncbi_scerevisiae.tsv | sed "s/', '/,/g" | sed "s/\['//g" | sed "s/']//g" > ${BASE_PATH}ncbi_scerevisiae_clean.tsv
