import pandas as pd
import numpy as np

from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/12_additional_work_3/" # CHANGE HERE

# fasta file input
refseq_file = Path(PurePath(base_path, f"saccharomyces_cerevisiae_pc_transcripts.fa"))

# for every record in the file
transcript_names = []
refseq_transcript_ids = []
motif_underrep1s = []
motif_underrep2s = []
motif_underrep3s = []
mut_suscept1s = []
mut_suscept2s = []
mut_suscept3s = []
motif_underrep_means = []
mut_suscept_means = []
motif_underrep_stds = [] 
mut_suscept_stds = []
for cnt, seq_record in enumerate(SeqIO.parse(refseq_file, "fasta")):    
    seq_id = seq_record.id
    seq_id_list = seq_id.split("_")
   
    filename = seq_id.replace("|","_")
    filename = filename.replace("-","_")
    filename = filename.replace(":","_")
    filename = f"{filename}_n3results.txt"

    run1_file = Path(PurePath(base_path, f"yeast", f"cdur_yeast_run1", filename))
    run2_file = Path(PurePath(base_path, f"yeast", f"cdur_yeast_run2", filename))
    run3_file = Path(PurePath(base_path, f"yeast", f"cdur_yeast_run3", filename))
    if run1_file.is_file() and run2_file.is_file() and run3_file.is_file():
        rlt1 = pd.read_table(run1_file, sep="\t", header=None, names=["Property", "Value"])
        rlt2 = pd.read_table(run2_file, sep="\t", header=None, names=["Property", "Value"])
        rlt3 = pd.read_table(run3_file, sep="\t", header=None, names=["Property", "Value"])
        
        motif_underrep1 = rlt1.loc[rlt1["Property"] == "belowT_C_"]["Value"].values[0]
        motif_underrep2 = rlt2.loc[rlt2["Property"] == "belowT_C_"]["Value"].values[0]
        motif_underrep3 = rlt3.loc[rlt3["Property"] == "belowT_C_"]["Value"].values[0]

        mut_suscept1 = 1-rlt1.loc[rlt1["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]
        mut_suscept2 = 1-rlt2.loc[rlt2["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]
        mut_suscept3 = 1-rlt3.loc[rlt3["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]

        motif_underrep_avg = np.mean([motif_underrep1, motif_underrep2, motif_underrep3])
        mut_suscept_avg = np.mean([mut_suscept1, mut_suscept2, mut_suscept3])

        motif_underrep_std = np.std([motif_underrep1, motif_underrep2, motif_underrep3])
        mut_suscept_std = np.std([mut_suscept1, mut_suscept2, mut_suscept3])

        refseq_id = "_".join(seq_id_list[:2])
        gene_symbol = "-".join(seq_id_list[5:-1])
        transcript_names.append(gene_symbol)
        refseq_transcript_ids.append(refseq_id)

        motif_underrep1s.append(motif_underrep1)
        motif_underrep2s.append(motif_underrep2)
        motif_underrep3s.append(motif_underrep3)
        mut_suscept1s.append(mut_suscept1)
        mut_suscept2s.append(mut_suscept2)
        mut_suscept3s.append(mut_suscept3)
        motif_underrep_means.append(motif_underrep_avg)
        mut_suscept_means.append(mut_suscept_avg)
        motif_underrep_stds.append(motif_underrep_std)
        mut_suscept_stds.append(mut_suscept_std)
    
df = pd.DataFrame(
    {
        "Transcript_name": transcript_names,
        "RefSeq_transcript_id": refseq_transcript_ids,
        "Motif_under-representation_run1" : motif_underrep1s,
        "Motif_under-representation_run2" : motif_underrep2s,
        "Motif_under-representation_run3" : motif_underrep3s,
        "Mutational_susceptibility_run1" : mut_suscept1s,
        "Mutational_susceptibility_run2" : mut_suscept2s,
        "Mutational_susceptibility_run3" : mut_suscept3s,
        "Motif_under-representation" : motif_underrep_means,
        "Mutational_susceptibility" : mut_suscept_means,
        "Motif_under-representation_std" : motif_underrep_stds,
        "Mutational_susceptibility_std" : mut_suscept_stds
    }
)
    
df.to_csv(Path(PurePath(base_path, f"cdur.saccharomyces_cerevisiae_pc_transcripts.tc.csv")), index=False)
