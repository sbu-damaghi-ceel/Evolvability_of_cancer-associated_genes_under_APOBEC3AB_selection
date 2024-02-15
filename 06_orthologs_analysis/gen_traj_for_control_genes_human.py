import sys

from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/sinabro") # CHANGE HERE
import sinabro.sinabro as snbr

import multiprocessing as mp

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
fasta_file = Path(PurePath(base_path, "control_genes_human.fa"))

def sim_traj(i, transcript_name, seq, transcript_path):
            traj = snbr.Trajectory(transcript_name, seq)
            traj.autofill(condition="nonsynonymous", method="mut_type", mut_type="T[C>T]")
            traj.save(transcript_path, prefix=transcript_name, suffix=str(i+1))

n_traj_per_transcript = 10
n_run = 1
n_processor = 10

for run in range(n_run):
    print(f"Run {run+1} starts...")
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        transcript_name = seq_record.id
        seq = seq_record.seq

        transcript_path = Path(PurePath(
            base_path, "control_genes_human", transcript_name,
            f"run_{str(run+1)}"
            ))

        transcript_path.mkdir(parents=True, exist_ok=True)
        trajs_args = []

        print(f"Gene {transcript_name} simulation starts...")
        for i in range(n_traj_per_transcript):
            trajs_args.append((i, transcript_name, seq, transcript_path))
        with mp.Pool(processes=n_processor) as pool:
            pool.starmap(sim_traj, trajs_args)
        print(f"Gene {transcript_name} simulation done...")
        print("")

    print(f"Run {run+1} done...")


            
