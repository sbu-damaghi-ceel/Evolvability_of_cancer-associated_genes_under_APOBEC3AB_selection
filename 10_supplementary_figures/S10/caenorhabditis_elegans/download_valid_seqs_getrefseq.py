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
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
ncbi_celegans_file = Path(PurePath(base_path, "ncbi_celegans_clean.tsv"))

ncbi_celegans = pd.read_csv(ncbi_celegans_file, sep="\t")
gene_ids = ncbi_celegans.loc[ncbi_celegans["type"] == "PROTEIN_CODING"]["gene_id"].to_list()

genes_files_path = Path(PurePath(base_path, "caenorhabditis_elegans_genes"))
genes_files_path.mkdir(parents=True, exist_ok=True)

last_index = len(gene_ids)-1
max_cnt = 9000
is_done = True
starts = []
ends = []
ii = 0
while is_done:
    starts.append(ii)
    if ii+9000 >= last_index:
        ends.append(last_index)
        break
    else:
        ends.append(ii+9000)
        ii = ii+9000

print(starts)
print(ends)

seq_records = []

for batch_idx, start_idx in enumerate(starts):
    end_idx = ends[batch_idx]

    print("Batch efetch gene ids...")
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

                new_id = f"{ref_id}_Caenorhabditis_elegans_{gene_symbol}_{'_'.join(transcript_variant.split())}"
            else:
                if "/" in gene_symbol:
                    gene_symbol = gene_symbol.replace("/", ".")

                new_id = f"{ref_id}_Caenorhabditis_elegans_{gene_symbol}"

            if "N" in subrecord[idx]["GBSeq_sequence"][from_idx-1:to_idx].upper():
                pass
            else: 
                seq_record = SeqRecord(Seq(subrecord[idx]["GBSeq_sequence"][from_idx-1:to_idx].upper()), id=new_id, description="")
                SeqIO.write(seq_record, PurePath(base_path, "caenorhabditis_elegans_genes", f"{new_id}.fa"), "fasta")
                seq_records.append(seq_record)  
                print(seq_record)
                print("")
    print("Saving all transcripts into separate files complete...!!")

print("Saving all transcripts into single file...")
output_file = PurePath(base_path, f"caenorhabditis_elegans_pc_transcripts.fa")
SeqIO.write(seq_records, output_file, "fasta")
print("Saving all transcripts into single file complete...!!")
