import pandas as pd
from pathlib import Path, PurePath

base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

cdur_file = PurePath(base_path, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.csv")
cancer_genes_file = PurePath(base_path, "cancer_gene_census.csv")
control_genes_file = PurePath(base_path, "ncbi_human_clean_no_mutation_protein_coding_no_LOC_no_MT.tsv")

cdur = pd.read_csv(cdur_file)
cancer_genes = pd.read_csv(cancer_genes_file)
control_genes = pd.read_csv(control_genes_file)

gencode_lst = [item[:len(item)-4] for item in cdur["Transcript_name"].to_list()]
cancer_gene_list = cancer_genes["Gene Symbol"].to_list()
control_gene_list = control_genes["symbol"].to_list()

index_lst = [key for key, val in enumerate(gencode_lst) if val in set(control_gene_list)]
index_lst2 = [key for key, val in enumerate(gencode_lst) if val in set(cancer_gene_list)]

gene_class = []
for i in range(len(cdur["Transcript_name"].to_list())):
    if i in index_lst:
        gene_class.append("control")
    elif i in index_lst2:
        gene_class.append("cancer")
    else:
        gene_class.append("regular")

cdur["gene_class"] = gene_class

output_file = PurePath(base_path, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.geneclass.csv")
cdur.to_csv(output_file, index=False, sep="\t")



