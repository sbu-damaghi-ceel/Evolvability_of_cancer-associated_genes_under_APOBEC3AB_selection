# import packages
from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np 
import re
#import scipy as 
# parameter
looplen = 3
maxstem = 20 
threshold = 15

left_plus = maxstem*"N"
right_plus = maxstem*"N"
test_seq = "\
ATGAAGAAGGTAACTGCAGAGGCTATTTCCTGGAATGAATCAACGAGTGAAACGAATAAC\
TCTATGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTCCAGACCTTC\
CTATTTATGTTGTTTTTTGTATTCTATGGAGGAATCGTGTTTGGAAACCTTCTTATTGTC\
ATAACAGTGGTATCTGACTCCCACCTTCACTCTCCCATGTACTTCCTGCTAGCCAACCTC\
TCACTCATTGATCTGTCTCTGTCTTCAGTCACAGCCCCCAAGATGATTACTGACTTTTTC\
AGCCAGCGCAAAGTCATCTCTTTCAAGGGCTGCCTTGTTCAGATATTTCTCCTTCACTTC\
TTTGGTGGGAGTGAGATGGTGATCCTCATAGCCATGGGCTTTGACAGATATATAGCAATA\
TGCAAGCCCCTACACTACACTACAATTATGTGTGGCAACGCATGTGTCGGCATTATGGCT\
GTCACATGGGGAATTGGCTTTCTCCATTCGGTGAGCCAGTTGGCGTTTGCCGTGCACTTA\
CTCTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAA\
CTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTGCTC\
ACTGTGTGTTCTTTTGTTCTTCTAATCATCTCATACACTATCATCCTAATGACCATCCAG\
CATCGCCCTTTAGATAAGTCGTCCAAAGCTCTGTCCACTTTGACTGCTCACATTACAGTA\
GTTCTTTTGTTCTTTGGACCATGTGTCTTTATTTATGCCTGGCCATTCCCCATCAAGTCA\
TTAGATAAATTCCTTGCTGTATTTTATTCTGTGATCACCCCTCTCTTGAACCCAATTATA\
TACACACTGAGGAACAAAGACATGAAGACGGCAATAAGACAGCTGAGAAAATGGGATGCA\
CATTCTAGTGTAAAGTTTTAG"

new_test_seq = left_plus+test_seq+right_plus
#rev_seq = str(Seq(new_test_seq).reverse_complement())

# function
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

# array1 = np.array(find_G_indices_GpB(new_test_seq))
# array2 = np.array(find_C_indices_VpC(new_test_seq))
# print(np.array(find_G_indices_GpB(new_test_seq)))
# new_way = len(new_test_seq) -1 - np.array(find_C_indices_VpC(rev_seq))
# print(new_way[::-1])
# print(np.array(find_C_indices_VpC(new_test_seq)))

# print(np.intersect1d(array1, array2))

# count APOBEC target VpC in hairpins
def count_hairpin_apobec_target(seq, looplen = 3, maxstem = 20, threshold = 15): 
    # (+) strand
    score_plus_strand = []
    for i in find_C_indices_VpC(seq):
        up_stem = new_test_seq[i-2-maxstem:i-2]
        down_stem = new_test_seq[i+1:i+1+maxstem]

        is_pair = np.array([1 if x==5 else 0 for x in np.array(convert_DNA_to_int(up_stem)) + np.array(convert_DNA_to_int(down_stem)[::-1])])
        pair_type = np.array(convert_DNA_to_int(up_stem))*np.array(convert_DNA_to_int(down_stem)[::-1]) - 3

        loop_score = np.sum(is_pair*pair_type)
        score_plus_strand.append(loop_score)

        #print(loop_score)
    #print(score_plus_strand)
    # (-) strand
    rev_seq = str(Seq(seq).reverse_complement())
    score_minus_strand = []
    for i in find_C_indices_VpC(rev_seq):
        up_stem = new_test_seq[i-maxstem:i]
        down_stem = new_test_seq[i+3:i+3+maxstem]

        is_pair = np.array([1 if x==5 else 0 for x in np.array(convert_DNA_to_int(up_stem)) + np.array(convert_DNA_to_int(down_stem)[::-1])])
        pair_type = np.array(convert_DNA_to_int(up_stem))*np.array(convert_DNA_to_int(down_stem)[::-1]) - 3

        loop_score = np.sum(is_pair*pair_type)
        score_minus_strand.append(loop_score)

        #print(loop_score)
    #print(score_minus_strand)
        
    #print(np.sum(np.array(score_plus_strand) >= threshold))
    #print(np.sum(np.array(score_minus_strand) >= threshold))
    
    return np.sum(np.array(score_plus_strand) >= threshold)+np.sum(np.array(score_minus_strand) >= threshold)
    
    

print(count_hairpin_apobec_target(new_test_seq))


