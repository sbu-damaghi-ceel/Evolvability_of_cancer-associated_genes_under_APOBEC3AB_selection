from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
control_genes_path = Path(PurePath(base_path, "control_genes_human"))

n_traj_per_transcript = 10

for path in control_genes_path.iterdir():
    if path.is_dir():
        genename = str(path).split("/")[-1]    

        run_dir = PurePath(str(path), "run_1")

        for i in range(1,n_traj_per_transcript+1,1):
            file_name = f"{genename}_trajectory_{i}.fa"
            out_name = f"fixed_{genename}_trajectory_{i}.fa"
            input_file = PurePath(run_dir, file_name)
            output_file = PurePath(run_dir, out_name)

            new_records = []
            cnt = 0
            for seq_record in SeqIO.parse(input_file, "fasta"):
                old_id = seq_record.id
                seq = seq_record.seq[1:-1]
                description = seq_record.description
                
                new_id = f"{old_id}_{i}_{cnt}"

                record = SeqRecord(seq, id=new_id, description="")
                new_records.append(record)
                cnt = cnt + 1

            SeqIO.write(new_records, output_file, "fasta")

