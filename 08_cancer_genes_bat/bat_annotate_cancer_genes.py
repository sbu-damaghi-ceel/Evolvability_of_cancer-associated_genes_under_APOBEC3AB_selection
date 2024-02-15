import pandas as pd
from pathlib import Path, PurePath
import sys

import re

base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

cdur_file = PurePath(base_path, "cdur.pteropus_alecto_pc_transcripts.tc.csv")
cancer_genes_file = PurePath(base_path, "palecto_cancer_genes.csv")
control_genes_file = PurePath(base_path, "palecto_control_genes.csv")

cdur = pd.read_csv(cdur_file)
cancer_genes = pd.read_csv(cancer_genes_file)
control_genes = pd.read_csv(control_genes_file)

control_genes_list = [id.split(".")[0] for id in control_genes["Transcript Accession"].to_list()]
cancer_genes_list = [id.split(".")[0] for id in cancer_genes["Transcript Accession"].to_list()]

gene_class = []
for i in range(len(cdur["Transcript_name"].to_list())):
    id = cdur.iloc[i]["RefSeq_transcript_id"].split(".")[0]
    if id in control_genes_list:
        gene_class.append("control")
    elif id in cancer_genes_list:
        gene_class.append("cancer")
    else:
        gene_class.append("regular")

cdur["gene_class"] = gene_class

output_file = PurePath(base_path, "cdur.pteropus_alecto_pc_transcripts.tc.geneclass.csv")
cdur.to_csv(output_file, index=False, sep="\t")



