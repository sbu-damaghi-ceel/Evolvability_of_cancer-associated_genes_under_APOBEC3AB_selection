from pathlib import Path, PurePath
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio import Entrez

import pandas as pd
import time
import numpy as np
import sys

Entrez.email = "joon-hyun.song@stonybrook.edu"

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/12_additional_work_3/" # CHANGE HERE
ncbi_scerevisiae_file = Path(PurePath(base_path, "ncbi_scerevisiae_clean.tsv"))

ncbi_scerevisiae = pd.read_csv(ncbi_scerevisiae_file, sep="\t")
gene_ids = ncbi_scerevisiae.loc[ncbi_scerevisiae["type"] == "PROTEIN_CODING"]["gene_id"].to_list()

last_index = len(gene_ids)-1
max_cnt = 1000
is_done = True
starts = []
ends = []
ii = 0
while is_done:
    starts.append(ii)
    if ii+max_cnt >= last_index:
        ends.append(last_index)
        break
    else:
        ends.append(ii+max_cnt)
        ii = ii+max_cnt

print(starts)
print(ends)

seq_records = []

for batch_idx, start_idx in enumerate(starts):
    end_idx = ends[batch_idx]

    print("Batch efetch gene ids...")
    print(gene_ids[start_idx:end_idx])
    handle = Entrez.efetch(db="gene", id=gene_ids[start_idx:end_idx], rettype="gb", retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    print("Batch efetch gene ids complete...!!")

    print(f"Number of records: {len(record)}")
    print("Extracting RefSeq ids...")
    ref_ids = []
    for i in range(len(record)):
        if "Gene-commentary_products" in record[i]["Entrezgene_locus"][0].keys():
            for product in record[i]["Entrezgene_locus"][0]["Gene-commentary_products"]:  
                if 'Gene-commentary_label' in product.keys():
                    transcript_variant = product['Gene-commentary_label']
                else:
                    transcript_variant = ""

                ref_id = product['Gene-commentary_accession']

                if ref_id[1] != "R":
                    ref_ids.append(ref_id)
    print("Extracting RefSeq ids complete...!!")

    ref_last_index = last_index = len(ref_ids)-1
    ref_is_done = True
    ref_starts = []
    ref_ends = []
    ref_ii = 0
    while ref_is_done:
        ref_starts.append(ref_ii)
        if ref_ii+9000 >= ref_last_index:
            ref_ends.append(ref_last_index)
            break
        else:
            ref_ends.append(ref_ii+9000)
            ref_ii = ref_ii+9000

    print(ref_starts)
    print(ref_ends)
    for ref_batch_idx, ref_start_idx in enumerate(ref_starts):
        ref_end_idx = ref_ends[ref_batch_idx]
        print("Batch efetch RefSeq ids...")
        subhandle = Entrez.efetch(db="nucleotide", id=ref_ids[ref_start_idx:ref_end_idx], rettype="gb", retmode="xml")
        subrecord = Entrez.read(subhandle)
        subhandle.close()
        print("Batch efetch RefSeq ids complete...!!")

        print(f"Number of subrecords: {len(subrecord)}")
        print("Saving all transcripts into separate files...")
        for idx in range(len(subrecord)):
            cds_indices = []
            for i in range(len(subrecord[idx]["GBSeq_feature-table"])):
                for key, value in subrecord[idx]["GBSeq_feature-table"][i].items():           
                    if key == "GBFeature_key" and value == "CDS":
                        ref_id = subrecord[idx]["GBSeq_feature-table"][i]["GBFeature_intervals"][0]["GBInterval_accession"]
                        gene_symbol = subrecord[idx]["GBSeq_feature-table"][i]["GBFeature_quals"][0]["GBQualifier_value"]               
                        cds_indices.append(i)

            if cds_indices:
                froms = []
                tos = []
                for i in cds_indices:
                    for key, value in subrecord[idx]["GBSeq_feature-table"][i].items():
                        if key == "GBFeature_intervals":
                            froms.append(int(value[0]["GBInterval_from"]))
                            tos.append(int(value[0]["GBInterval_to"]))

                from_idx = froms[0]
                to_idx = tos[0]

                if transcript_variant:
                    if "/" in gene_symbol:
                        gene_symbol = gene_symbol.replace("/", ".")

                    new_id = f"{ref_id}_Saccharomyces_cerevisiae_{gene_symbol}_{'_'.join(transcript_variant.split())}"
                else:
                    if "/" in gene_symbol:
                        gene_symbol = gene_symbol.replace("/", ".")

                    new_id = f"{ref_id}_Saccharomyces_cerevisiae_{gene_symbol}"

                if "N" in subrecord[idx]["GBSeq_sequence"][from_idx-1:to_idx].upper():
                    pass
                else: 
                    seq_record = SeqRecord(Seq(subrecord[idx]["GBSeq_sequence"][from_idx-1:to_idx].upper()), id=new_id, description="")
                    SeqIO.write(seq_record, PurePath(base_path, "saccharomyces_cerevisiae_genes", f"{new_id}.fa"), "fasta")
                    seq_records.append(seq_record)  
                    print(seq_record)
                    print("")
    print("Saving all transcripts into separate files complete...!!")

print("Saving all transcripts into single file...")
output_file = PurePath(base_path, f"saccharomyces_cerevisiae_pc_transcripts.fa")
SeqIO.write(seq_records, output_file, "fasta")
print("Saving all transcripts into single file complete...!!")
