import pandas as pd

from pathlib import Path, PurePath
from math import isnan

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

control_genes_file_name = "ncbi_human_clean_no_mutation_2_protein_coding_no_LOC_no_MT_ngfc.tsv" # need to change this
control_genes_file = Path(PurePath(base_path, control_genes_file_name))

output_file_name = "control_genes_id_for_orthologs.txt"
output_file = Path(PurePath(base_path, output_file_name))

control_genes = pd.read_csv(control_genes_file, sep="\t")

with open(output_file, "w") as f:
    f.write("Gene_symbol,GeneID\n")
    for i, row in control_genes.iterrows():
        if not isnan(row["gene_id"]):
            print(row["symbol"], row["gene_id"])
            f.write(f"{row['symbol']},{int(row['gene_id'])}\n")

