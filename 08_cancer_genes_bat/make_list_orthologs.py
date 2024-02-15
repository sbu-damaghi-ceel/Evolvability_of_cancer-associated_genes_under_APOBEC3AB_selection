from pathlib import Path, PurePath
import pandas as pd


# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
cancer_genes_path = Path(PurePath(base_path, "cancer_genes_orthologs"))
control_genes_path = Path(PurePath(base_path, "control_genes_orthologs"))

species_of_interest = [
    "Pteropus alecto"
]

species_short = [
    "palecto"
]

for i, species in enumerate(species_of_interest):
    infos = pd.DataFrame()

    for path in cancer_genes_path.iterdir():
        if path.is_dir():
            genename = str(path).split("/")[-1]

            data_table_file = PurePath(path, "ncbi_dataset", "data", "transcript_protein.tsv")
            data_table = pd.read_csv(data_table_file, sep="\t")
               
            new_infos = data_table[data_table["Taxonomic Name"] == species]

            if not new_infos.empty:
                if infos.empty:
                    infos = new_infos
                else:
                    infos = pd.concat([infos, new_infos], ignore_index=True)

    out_file = Path(PurePath(base_path, f"{species_short[i]}_cancer_genes.csv"))
    infos.to_csv(out_file, index=False)

    print(infos)


for i, species in enumerate(species_of_interest):
    infos = pd.DataFrame()

    for path in control_genes_path.iterdir():
        if path.is_dir():
            genename = str(path).split("/")[-1]

            data_table_file = PurePath(path, "ncbi_dataset", "data", "transcript_protein.tsv")
            data_table = pd.read_csv(data_table_file, sep="\t")
               
            new_infos = data_table[data_table["Taxonomic Name"] == species]

            if not new_infos.empty:
                if infos.empty:
                    infos = new_infos
                else:
                    infos = pd.concat([infos, new_infos], ignore_index=True)

    out_file = Path(PurePath(base_path, f"{species_short[i]}_control_genes.csv"))
    infos.to_csv(out_file, index=False)

    print(infos)



