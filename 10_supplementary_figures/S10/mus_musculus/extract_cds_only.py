# import packages
from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

# GENCODE release number
release_num = 33

# make files' names and files' pathes
input_file_name = f"gencode.vM{release_num}.pc_transcripts.fa"
input_file = Path(PurePath(base_path, input_file_name))

output_file_cdsonly_name = f"gencode.vM{release_num}.pc_transcripts.nopary.cdsonly.fa"
output_file_cdsonly = Path(PurePath(base_path, output_file_cdsonly_name))

nopass_names = []
for i in range(1,16,1):
    nopass_names.append(f"gencode.vM{release_num}.pc_transcripts.nopass.parity{i}.fa")

output_file_nopasses = []
for i in range(1,16,1):
    output_file_nopass = Path(
        PurePath(base_path, nopass_names[i-1])
    )
    output_file_nopasses.append(output_file_nopass)

# make empty lists to store new record
records_cdsonly = []
records_nopass_dict = {}
for i in range(1,16,1):
    records_nopass_dict[f"parity{i}"] = []

# for every record in the file
for cnt, seq_record in enumerate(SeqIO.parse(input_file, "fasta")):    
    seq_id = seq_record.id
    seq_id_list = seq_id.split("|")
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
