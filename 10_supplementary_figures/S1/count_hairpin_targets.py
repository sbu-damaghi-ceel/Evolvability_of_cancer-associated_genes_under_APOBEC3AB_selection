# import packages
import re

import multiprocessing as mp
import pandas as pd
import numpy as np 

from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
CODON_TABLE = CodonTable.standard_dna_table.forward_table
CODON_TABLE["TAA"] = "*"
CODON_TABLE["TAG"] = "*"
CODON_TABLE["TGA"] = "*"

# functions
def is_codon_synonymous(codon1, codon2): 
    if CODON_TABLE[codon1] == CODON_TABLE[codon2]:
        return True
    else:
        return False
    
def get_codon(seq, idx):
    if idx % 3 == 0:
        return seq[idx:idx+3]
    elif idx % 3 == 1:
        return seq[idx-1:idx+2]
    else:
        return seq[idx-2:idx+1]
    
def mutate_C_to_T(codon, idx):
    codon_list = [*codon]
    if codon_list[idx%3] == "C":
        codon_list[idx%3] = "T"
    else:
        codon_list[idx%3] = "A"

    return ''.join(codon_list)

def correct_index(idx, seq, strand, maxstem):
    if isinstance(idx, int):
        if strand == "+":
            return idx-maxstem
        else:
            seq_len = len(seq) - 2*maxstem
            return seq_len-1-(idx-maxstem)

    elif isinstance(idx, list):
        if strand == "+":
            return [i-maxstem for i in idx]
        else:
            seq_len = len(seq) - 2*maxstem
            return [seq_len-1-(i-maxstem) for i in idx]
    else:
        raise TypeError(
            "idx should be int or list"
        )

# change DNA sequence to integers
def convert_DNA_to_int(seq):
    mapping = {'N': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
    result = [mapping[nt] for nt in seq]
    return result

# find C of VpC and G of GpB indices
def find_C_indices_VpC(seq):
    patterns = ["AC", "CC", "GC"]
    patterns_regex = re.compile('|'.join('(?={0})'.format(re.escape(p)) for p in patterns))

    return [match.start()+1 for match in patterns_regex.finditer(seq)]

def find_G_indices_GpB(seq):
    patterns = ["GT", "GC", "GG"]
    patterns_regex = re.compile('|'.join('(?={0})'.format(re.escape(p)) for p in patterns))

    return [match.start() for match in patterns_regex.finditer(seq)]

# count APOBEC target VpC in hairpins
def score_stem(up_stem, down_stem):
    # How to compute pairing
        # is_pair:      check if the pair G:C or A:T pair
        # pair_type:    mark G:C and A:T pair differently
        #
        # if we map as 'A': 1, 'C': 2, 'G': 3, 'T': 4,
        # only G:C and A:T pairing will give sum of 5
        # 1 if the sum is 5, 0 otherwise 
        #                            --> is_pair
        #
        # in addition, if multiply mapped number,
        # G:C pair gives 6, A:T pair gives 4
        # substract 3 will result in 3 and 1, respectively 
        #                                              --> pair_type
        # elementwise multiplication of is_pair and pair_type gives np.array
        # where G:C pair marked as 3, A:T pair marked as 1, and 0 otherwise
        is_pair = np.array([1 if x==5 else 0 for x in np.array(convert_DNA_to_int(up_stem)) + np.array(convert_DNA_to_int(down_stem)[::-1])])
        pair_type = np.array(convert_DNA_to_int(up_stem))*np.array(convert_DNA_to_int(down_stem)[::-1]) - 3

        stem_score = np.sum(is_pair*pair_type)

        return stem_score


def search_hairpin_apobec_target(seq, looplen = 3, target_pos_in_loop = 3, maxstem = 20, threshold = 15): 
    # initialize out dict
    result = {
        'indices': [],
        'strands': [],
        'scores':  [],
        'threshold': threshold,
        'total_count': 0,
        'indices_target': [],
        'strands_target': [],
        'scores_target':  []
        }

    rev_seq = str(Seq(seq).reverse_complement())

    indices_plus  = find_C_indices_VpC(seq)
    indices_minus = find_C_indices_VpC(rev_seq)

    # (+) strand
    score_plus_strand = []
    for i in indices_plus:
        up_stem = seq[(i-target_pos_in_loop)+1-maxstem:(i-target_pos_in_loop)+1]
        down_stem = seq[i+(looplen-target_pos_in_loop)+1:i+(looplen-target_pos_in_loop)+1+maxstem]

        score_plus_strand.append(score_stem(up_stem, down_stem))

    # (-) strand
    score_minus_strand = []
    for i in indices_minus:
        up_stem = rev_seq[(i-target_pos_in_loop)+1-maxstem:(i-target_pos_in_loop)+1]
        down_stem = rev_seq[i+(looplen-target_pos_in_loop)+1:i+(looplen-target_pos_in_loop)+1+maxstem]

        score_minus_strand.append(score_stem(up_stem, down_stem))

    indices_plus_corrected = correct_index(indices_plus, seq, "+", maxstem)
    indices_minus_corrected = correct_index(indices_minus, seq, "-", maxstem)

    indices = np.array(indices_plus_corrected + indices_minus_corrected)
    strands = np.array(["+"]*len(indices_plus_corrected) + ["-"]*len(indices_minus_corrected))
    scores  = np.array(score_plus_strand + score_minus_strand)

    sort_idx = np.argsort(indices)

    result['indices'] = indices[sort_idx]
    result['strands'] = strands[sort_idx]
    result['scores']  = scores[sort_idx]
    result['total_count'] = np.sum(scores >= threshold)

    indices_over_threshold = np.where(result['scores'] >= threshold)[0]

    result['indices_target'] = result['indices'][indices_over_threshold]
    result['strands_target'] = result['strands'][indices_over_threshold]
    result['scores_target']  = result['scores'][indices_over_threshold]
    
    return result

def is_C_to_T_nonsynonymous(seq, i):
    ori_codon = get_codon(seq, i)
    mut_codon = mutate_C_to_T(ori_codon, i)

    return not is_codon_synonymous(ori_codon, mut_codon)

   
def get_CDUR_stats(seq, looplen = 3, target_pos_in_loop = 3, maxstem = 20, threshold = 15):
    left_plus = (maxstem+looplen+1)*"N"
    right_plus = (maxstem+looplen+1)*"N"
    new_seq = left_plus+seq+right_plus

    search_result = search_hairpin_apobec_target(new_seq, looplen, target_pos_in_loop, maxstem, threshold)
    indices = search_result['indices_target']

    # a. number of target
    V_C_loop = len(indices)
    # b. number of nonsynonymous
    nonsynonymity = list(map(is_C_to_T_nonsynonymous, V_C_loop*[seq], indices))
    repTrV_C_loop = nonsynonymity.count(False)
    # c. b/a
    if V_C_loop == 0:
        repTrFracV_C_loop = 0
    else:
        repTrFracV_C_loop = repTrV_C_loop/V_C_loop

    return (V_C_loop, repTrV_C_loop, repTrFracV_C_loop)

def get_CDUR_stats_from_fasta(file_path, *args, **kwargs):
    looplen = 3
    target_pos_in_loop = 3
    maxstem = 20 
    threshold = 15
    if "looplen" in kwargs.keys():
        looplen = kwargs["looplen"]
    if "target_pos_in_loop" in kwargs.keys():
        target_pos_in_loop = kwargs["target_pos_in_loop"]
    if "maxstem" in kwargs.keys():
        maxstem = kwargs[maxstem]
    if "threshold" in kwargs.keys():
        threshold = kwargs["threshold"]

    V_C_loops = []
    repTrV_C_loops = [] 
    repTrFracV_C_loops = []
    for seq_record in SeqIO.parse(file_path, "fasta"):
        V_C_loop, repTrV_C_loop, repTrFracV_C_loop = get_CDUR_stats(str(seq_record.seq), looplen, target_pos_in_loop, maxstem, threshold)

        V_C_loops.append(V_C_loop)
        repTrV_C_loops.append(repTrV_C_loop) 
        repTrFracV_C_loops.append(repTrFracV_C_loop)

    belowV_C_loop = np.sum(np.array(V_C_loops[1:]) < V_C_loops[0])/1000
    repTr_belowV_C_loop = np.sum(np.array(repTrV_C_loops[1:]) < repTrV_C_loops[0])/1000
    repTrFrac_belowV_C_loop = np.sum(np.array(repTrFracV_C_loops[1:]) < repTrFracV_C_loops[0])/1000

    if args:
        return args, belowV_C_loop, repTr_belowV_C_loop, repTrFrac_belowV_C_loop
    else:
        return belowV_C_loop, repTr_belowV_C_loop, repTrFrac_belowV_C_loop

####################################################
# parameter
looplen = 3
target_pos_in_loop = 3
maxstem = 20 
threshold = 15

print("Prepare to compute...")
# paths for files and directory
# FILE_PATH - location of fasta containing transcripts sequences, also output is saved here. 
# RUN_PATH  - location of CDUR output directory, where *_n3.fas files are in.
#             *_n3.fas contains a given sequence and 1,000 sequences shuffled using n3 method.
#
FILE_PATH = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/15_additional_work_6"
RUN_PATH  = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/15_additional_work_6/a3_stemloop/run1"
orig_fasta_file = PurePath(FILE_PATH, "gencode.v40.pc_transcripts.nopary.cdsonly.refseq.fa")
seq_records = list(SeqIO.parse(orig_fasta_file, "fasta"))

# prepare multiprocessing
kwargs = {
    "looplen": looplen,
    "target_pos_in_loop": target_pos_in_loop,
    "maxstem": maxstem, 
    "threshold": threshold
}

data_for_multiprocessing = []
for i in range(len(seq_records)):
    transcript_name = seq_records[i].id.split("|")[4]
    ensembl_transcript_id = seq_records[i].id.split("|")[0]
    ensembl_gene_id = seq_records[i].id.split("|")[1]

    file_name = seq_records[i].id.replace("|", "_").replace(":", "_").replace("-", "_") + "_n3.fasta"
    file_path = PurePath(RUN_PATH, file_name)

    data_for_multiprocessing.append((file_path, transcript_name, ensembl_transcript_id, ensembl_gene_id, kwargs))
print("Preparing done...")

print("Compute CDUR statistics... This takes long...")
# compute using multiprocessing
n_processor = 10
with mp.Pool(processes=n_processor) as pool:
    results = pool.starmap(get_CDUR_stats_from_fasta, data_for_multiprocessing)
print("Compute CDUR statistics done...")

print("Saving results...")
# make result into Pandas data frame and save into csv
output_col1 = [] # Transcript_name
output_col2 = [] # Ensembl_transcript_id
output_col3 = [] # Ensembl_gene_id
output_col4 = [] # Motif_under-representation
output_col5 = [] # Mutational_susceptibility
for rlt in results:
    transcript_name = rlt[0][0]
    ensembl_transcript_id = rlt[0][1]
    ensembl_gene_id = rlt[0][2]
    motif_underref = rlt[1]
    mut_susceptibility = 1-rlt[3]

    output_col1.append(transcript_name)
    output_col2.append(ensembl_transcript_id)
    output_col3.append(ensembl_gene_id)
    output_col4.append(motif_underref)
    output_col5.append(mut_susceptibility)

out_df = pd.DataFrame(
    {
        "Transcript_name": output_col1,
        "Ensembl_transcript_id": output_col2,
        "Ensembl_gene_id": output_col3,
        "Motif_under-representation": output_col4,
        "Mutational_susceptibility": output_col5
    }
)

outfile_path = PurePath(FILE_PATH, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.vcloop.csv")
out_df.to_csv(outfile_path, index=False)
print("Results saved...")
print("Completed!")