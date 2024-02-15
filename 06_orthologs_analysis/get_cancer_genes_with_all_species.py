from pathlib import Path, PurePath

import pandas as pd

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
cancer_genes_path = Path(PurePath(base_path, "cancer_genes_orthologs"))

species_of_interest = [
    "Homo sapiens",
    "Pan troglodytes",
    "Canis lupus familiaris",
    "Mus musculus",
    "Physeter catodon",
    "Loxodonta africana",
    "Myotis lucifugus",
    "Gallus gallus",
    "Xenopus tropicalis",
    "Petromyzon marinus"
]

gene_names = []
gene_ids = []

for path in cancer_genes_path.iterdir():
    if path.is_dir():
        genename = str(path).split("/")[-1]

        fasta_file = PurePath(path, "ncbi_dataset", "data", "gene.fna")
        infos_file = PurePath(path, "ncbi_dataset", "data", "transcript_protein.tsv")

        infos = pd.read_csv(infos_file, sep="\t")

        species_in_orthologs = infos["Taxonomic Name"].to_list()
        species_checklist = [False, False, False, False, False, False, False, False, False, False]
        for i, species in enumerate(species_of_interest):
            if species_of_interest[i] in species_in_orthologs:
                species_checklist[i] = True

        if sum(species_checklist) == 10:
            print(genename)
            gene_names.append(genename)

df = pd.DataFrame(
        {
        "Gene_name": gene_names
    }
)

df.to_csv(str(PurePath(base_path, "cancer_genes_having_ten_orthologs.csv")), index=False)