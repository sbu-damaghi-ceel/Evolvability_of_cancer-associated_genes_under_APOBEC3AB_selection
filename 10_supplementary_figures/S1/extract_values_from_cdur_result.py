import pandas as pd
import numpy as np

from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/15_additional_work_6/a3_stemloop" # CHANGE HERE

# fasta file input
refseq_file = Path(PurePath(base_path, f"gencode.v40.pc_transcripts.nopary.cdsonly.refseq.fa"))

# for every record in the file
transcript_names = []
ensembl_transcript_ids = []
ensembl_gene_ids = []
motif_underreps = []
mut_suscepts = []
for cnt, seq_record in enumerate(SeqIO.parse(refseq_file, "fasta")):    
    seq_id = seq_record.id
    seq_id_list = seq_id.split("|")
   
    filename = seq_id.replace("|","_")
    filename = filename.replace("-","_")
    filename = filename.replace(":","_")
    filename = f"{filename}_n3results.txt"

    run1_file = Path(PurePath(base_path, f"run1", filename))
    if run1_file.is_file():
        rlt1 = pd.read_table(run1_file, sep="\t", header=None, names=["Property", "Value"])
        
        motif_underrep1 = rlt1.loc[rlt1["Property"] == "belowV_C_"]["Value"].values[0]
        mut_suscept1 = 1-rlt1.loc[rlt1["Property"] == "repTrFrac_belowV_C_"]["Value"].values[0]

        transcript_id = seq_id_list[0]
        gene_id = seq_id_list[1]
        gene_symbol = seq_id_list[4]

        transcript_names.append(gene_symbol)
        ensembl_transcript_ids.append(transcript_id)
        ensembl_gene_ids.append(gene_id)

        motif_underreps.append(motif_underrep1)
        mut_suscepts.append(mut_suscept1)

    
df = pd.DataFrame(
    {
        "Transcript_name": transcript_names,
        "Ensembl_transcript_id": ensembl_transcript_ids,
        "Ensembl_gene_id": ensembl_gene_ids,
        "Motif_under.representation" : motif_underreps,
        "Mutational_susceptibility" : mut_suscepts
    }
)
    
df.to_csv(Path(PurePath(base_path, f"cdur.gencode.v40.pc_transcripts.vc.csv")), index=False)
