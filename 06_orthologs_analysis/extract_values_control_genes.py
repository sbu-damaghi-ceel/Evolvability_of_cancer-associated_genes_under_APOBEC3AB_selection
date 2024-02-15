import pandas as pd

from pathlib import Path, PurePath

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
orthologs_path = Path(PurePath(base_path, "all_control_genes_orthologs"))
human_trajs_path = Path(PurePath(base_path, "all_trajs_all_control_genes_human"))

genenames = []
species = []
motif_underreps = []
mut_suscepts = []
for path in orthologs_path.iterdir():
    if path.suffix == ".txt":
        results = pd.read_table(path, sep="\t", header=None, names=["Property", "Value"])

        motif_underrep = results.loc[results["Property"] == "belowT_C_"]["Value"].values[0]
        mut_suscept = 1-results.loc[results["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]

        genenames.append(path.stem.split("_")[2].upper())
        species.append("_".join([path.stem.split("_")[0], path.stem.split("_")[1]]))
        motif_underreps.append(motif_underrep)
        mut_suscepts.append(mut_suscept)

df = pd.DataFrame(
    {
        "Gene_name": genenames,
        "Species": species,
        "Motif_under-representation" : motif_underreps,
        "Mutational_susceptibility" : mut_suscepts
    }
)
    
df.to_csv(Path(PurePath(base_path, f"all_control_genes_orthologs_cdur.csv")), index=False)

n_traj_per_transcript = 10

genenames = []
trajs = []
traj_points = []
motif_underreps = []
mut_suscepts = []
for path in human_trajs_path.iterdir():
    if path.suffix == ".txt":
        results = pd.read_table(path, sep="\t", header=None, names=["Property", "Value"])

        motif_underrep = results.loc[results["Property"] == "belowT_C_"]["Value"].values[0]
        mut_suscept = 1-results.loc[results["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]

        genenames.append(path.stem.split("_")[2].upper())
        trajs.append(path.stem.split("_")[3])
        traj_points.append(path.stem.split("_")[4])
        motif_underreps.append(motif_underrep)
        mut_suscepts.append(mut_suscept)

df = pd.DataFrame(
    {
        "Gene_name": genenames,
        "Trajectory": trajs,
        "Traj_points": traj_points,
        "Motif_under-representation" : motif_underreps,
        "Mutational_susceptibility" : mut_suscepts
    }
)
    
df.to_csv(Path(PurePath(base_path, f"all_trajs_all_control_genes_human_cdur.csv")), index=False)

