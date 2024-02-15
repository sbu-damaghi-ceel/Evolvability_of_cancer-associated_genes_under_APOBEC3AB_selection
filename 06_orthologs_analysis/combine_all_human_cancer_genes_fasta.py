from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
cancer_genes_path = Path(PurePath(base_path, "cancer_genes_human"))

n_traj_per_transcript = 10

new_records = []
for path in cancer_genes_path.iterdir():
    if path.is_dir():
        genename = str(path).split("/")[-1]    

        run_dir = PurePath(str(path), "run_1")

        for i in range(1,n_traj_per_transcript+1,1):
            file_name = f"fixed_{genename}_trajectory_{i}.fa"
            input_file = PurePath(run_dir, file_name)

            cnt = 0
            for seq_record in SeqIO.parse(input_file, "fasta"):
                new_records.append(seq_record)
                cnt = cnt + 1

output_file = PurePath(base_path, "all_trajs_all_cancer_genes_human.fa")
SeqIO.write(new_records, output_file, "fasta")