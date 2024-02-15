from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
cancer_genes_path = Path(PurePath(base_path, "cancer_genes_orthologs"))

cancer_genes_having_ten_orthologs_file = Path(PurePath(base_path, "cancer_genes_having_ten_orthologs.csv"))
cancer_genes_having_ten_orthologs = pd.read_csv(cancer_genes_having_ten_orthologs_file)

new_records = []
for path in cancer_genes_path.iterdir():
    if path.is_dir():
        genename = str(path).split("/")[-1]

        if genename in cancer_genes_having_ten_orthologs["Gene_name"].to_list(): 
            file_name = f"ten_species_{genename.lower()}_cdur.fa"
            input_file = PurePath(path, file_name)

            for seq_record in SeqIO.parse(input_file, "fasta"):
                seq = seq_record.seq

                new_record = SeqRecord(seq, id=seq_record.id, description="")
                new_records.append(new_record)

output_file = PurePath(base_path, "all_cancer_genes_orthologs.fa")
SeqIO.write(new_records, output_file, "fasta")