import pandas as pd

from pathlib import Path, PurePath
from math import isnan

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

cancer_census_file_name = "cancer_gene_census.csv"
cancer_census_file = Path(PurePath(base_path, cancer_census_file_name))

output_file_name = "cancer_gene_census_id_for_orthologs.txt"
output_file = Path(PurePath(base_path, output_file_name))

cancer_census = pd.read_csv(cancer_census_file)

with open(output_file, "w") as f:
    f.write("Gene_symbol,GeneID\n")
    for i, row in cancer_census.iterrows():
        # we found 4 genes with no Entrez GeneId,
        # manual search via NCBI show
        # MDS2, MALAT1 as ncRNA
        # DUX4L1, HMGN2P46 as pseudogene
        # therefore concluded to ignore these 4 genes.

        if not isnan(row["Entrez GeneId"]):
            print(row["Gene Symbol"], row["Entrez GeneId"])
            f.write(f"{row['Gene Symbol']},{int(row['Entrez GeneId'])}\n")

