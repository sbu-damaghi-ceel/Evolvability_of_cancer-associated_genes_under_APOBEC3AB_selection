# import packages
from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

# GENCODE release number
release_num = 40

# make files' names and files' pathes
input_file_name = f"gencode.v{release_num}.pc_transcripts.nopary.cdsonly.fa"
input_file = Path(PurePath(base_path, input_file_name))

canonical_file_name = f"gencode.v{release_num}.pc_transcripts.nopary.cdsonly.canonical.id.csv"
canonical_file = Path(PurePath(base_path, canonical_file_name))
canonical = pd.read_csv(canonical_file)

refseq_file_name = f"gencode.v{release_num}.pc_transcripts.nopary.cdsonly.haverefseq.id.csv"
refseq_file = Path(PurePath(base_path, refseq_file_name))
refseq = pd.read_csv(refseq_file)

output_file_canonical_name = f"gencode.v{release_num}.pc_transcripts.nopary.cdsonly.canonical.fa"
output_file_canonical = Path(PurePath(base_path, output_file_canonical_name))

output_file_refseq_name = f"gencode.v{release_num}.pc_transcripts.nopary.cdsonly.refseq.fa"
output_file_refseq = Path(PurePath(base_path, output_file_refseq_name))

# make empty lists to store new record
records_canonical = []
records_refseq = []

# for every record in the file
for cnt, seq_record in enumerate(SeqIO.parse(input_file, "fasta")):    
    seq_id = seq_record.id
    seq_id_list = seq_id.split("|")

    ensembl_trancripts_id = seq_id_list[0].split(".")[0]

    # if record is canonical
    if ensembl_trancripts_id in canonical["x"].values:
        records_canonical.append(SeqRecord(Seq(seq_record.seq), id=seq_id ,description=""))

    # if recored has refseq id
    if ensembl_trancripts_id in refseq["x"].values:
        records_refseq.append(SeqRecord(Seq(seq_record.seq), id=seq_id ,description=""))

# save sequences
SeqIO.write(records_canonical, output_file_canonical, "fasta")
SeqIO.write(records_refseq, output_file_refseq, "fasta")

"""
    parity = 0
    
    if seq_id_list[7][0] == "C":
        cds_idx = 7
    elif seq_id_list[7][0] == "U":
        cds_idx = 8

    start = int(seq_id_list[cds_idx][4:].split("-")[0])
    end = int(seq_id_list[cds_idx][4:].split("-")[1])

    new_seq_cdsonly = seq_record.seq[start-1:end]
    record_cdsonly = SeqRecord(Seq(new_seq_cdsonly), id=seq_id ,description="")

    # check cds length is multiple of 3
    if len(new_seq_cdsonly)%3 != 0:
        parity = parity + 1
    
    # check cds starts with ATG
    if new_seq_cdsonly[0:3] != "ATG":
        parity = parity + 2
        
    # check cds ends with stop codon
    if new_seq_cdsonly[len(new_seq_cdsonly)-3:] not in ["TAA", "TAG", "TGA"]:
        parity = parity + 4
        
    # check duplicate PAR_Y
    if seq_id_list[0][len(seq_id_list[0])-5:] == "PAR_Y":
        parity = parity + 8
    
    # collect sequences that passes all check
    if parity == 0:              
        records_cdsonly.append(record_cdsonly)
    else:
        records_nopass_dict[f"parity{parity}"].append(seq_record)
    
# save passed sequences
SeqIO.write(records_cdsonly, output_file_cdsonly, "fasta")

# save sequences that are not passed by parity
for i in range(1,16,1):
    SeqIO.write(records_nopass_dict[f"parity{i}"], output_file_nopasses[i-1], "fasta")
"""





