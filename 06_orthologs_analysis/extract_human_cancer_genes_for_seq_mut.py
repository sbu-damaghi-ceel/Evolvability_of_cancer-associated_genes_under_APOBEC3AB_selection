from pathlib import Path, PurePath
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pandas as pd
import re

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

cancer_genes_having_ten_orthologs_file = Path(PurePath(base_path, "cancer_genes_having_ten_orthologs.csv"))
cancer_genes_having_ten_orthologs = pd.read_csv(cancer_genes_having_ten_orthologs_file)

out_seq_records = []
for gene in cancer_genes_having_ten_orthologs["Gene_name"].to_list():
    fasta_file = Path(PurePath(base_path, "cancer_genes_orthologs", gene, f"ten_species_{gene.lower()}.fa"))

    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        id = seq_record.id

        if re.match("^Homo", id):
            new_record = SeqRecord(seq_record.seq, id=id, description="")
            out_seq_records.append(new_record)

output_file = PurePath(base_path, f"cancer_genes_human.fa")
SeqIO.write(out_seq_records, output_file, "fasta")