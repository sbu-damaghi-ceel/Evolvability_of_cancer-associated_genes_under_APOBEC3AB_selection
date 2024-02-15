import pandas as pd
import numpy as np

from pathlib import Path, PurePath

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
result_path = Path(PurePath(base_path, "pteropus_alecto", "n3_results"))

# for every record in the file
transcript_names = []
gene_symbols = []
gene_types = []
transcript_variants = []
refseq_ids = []
motif_underrep1s = []
motif_underrep2s = []
motif_underrep3s = []
mut_suscept1s = []
mut_suscept2s = []
mut_suscept3s = []
motif_underrep_means = []
mut_suscept_means = []
motif_underrep_stds = [] 
mut_suscept_stds = []


for path in result_path.iterdir():
    if path.is_file():
        filename = str(path).split("/")[-1]
        refseq_id = "_".join(filename.split("_")[0:2])

        fullname = [i for i in filename.split("_") if i]
        #print(fullname)
        delim = "_".join(filename.split("__")[-1].split("_")[0:-1]+["n3results"])
        name_type_removed = filename.split(delim)[0]
        fullname = [i for i in name_type_removed.split("_") if i]

        #print(fullname)

        if fullname[0] == "XM" or fullname[0] == "XR":
            cleaned_name = fullname[5:]
        else:
            cleaned_name = fullname[4:]
        #print(cleaned_name)
        if "transcript" in cleaned_name and "variant" in cleaned_name:
            transcript_variant = "_".join(cleaned_name[-3:])
            if len(cleaned_name[-4]) == 1:
                gene_symbol = "-".join(cleaned_name[-5:])
                gene_name = "-".join(cleaned_name[:-5])
            else:
                gene_symbol = cleaned_name[-4]
                gene_name = "-".join(cleaned_name[:-4])

        else:
            transcript_variant = "NA"
            if len(cleaned_name[-1]) == 1:
                gene_symbol = "-".join(cleaned_name[-2:])
                gene_name = "-".join(cleaned_name[:-2])
            else:
                gene_symbol = cleaned_name[-1]
                gene_name = "-".join(cleaned_name[:-1])

        typename = "_".join(filename.split("__")[-1].split("_")[0:-1])
        print(f"{refseq_id}\t{gene_symbol}\t{gene_name}\t{transcript_variant}\t{typename}")  

        rlt = pd.read_table(path, sep="\t", header=None, names=["Property", "Value"])

        if not rlt.empty:
            motif_underrep1 = rlt.loc[rlt["Property"] == "belowT_C_"]["Value"].values[0]
            motif_underrep2 = rlt.loc[rlt["Property"] == "belowT_C_"]["Value"].values[0]
            motif_underrep3 = rlt.loc[rlt["Property"] == "belowT_C_"]["Value"].values[0]

            mut_suscept1 = 1-rlt.loc[rlt["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]
            mut_suscept2 = 1-rlt.loc[rlt["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]
            mut_suscept3 = 1-rlt.loc[rlt["Property"] == "repTrFrac_belowT_C_"]["Value"].values[0]

            motif_underrep_avg = np.mean([motif_underrep1, motif_underrep2, motif_underrep3])
            mut_suscept_avg = np.mean([mut_suscept1, mut_suscept2, mut_suscept3])

            motif_underrep_std = np.std([motif_underrep1, motif_underrep2, motif_underrep3])
            mut_suscept_std = np.std([mut_suscept1, mut_suscept2, mut_suscept3])


            transcript_names.append(gene_name)
            gene_symbols.append(gene_symbol)
            transcript_variants.append(transcript_variant)
            gene_types.append(typename)

            refseq_ids.append(refseq_id)

            motif_underrep1s.append(motif_underrep1)
            motif_underrep2s.append(motif_underrep2)
            motif_underrep3s.append(motif_underrep3)
            mut_suscept1s.append(mut_suscept1)
            mut_suscept2s.append(mut_suscept2)
            mut_suscept3s.append(mut_suscept3)
            motif_underrep_means.append(motif_underrep_avg)
            mut_suscept_means.append(mut_suscept_avg)
            motif_underrep_stds.append(motif_underrep_std)
            mut_suscept_stds.append(mut_suscept_std)


    
df = pd.DataFrame(
    {
        "Gene_symbol": gene_symbols,
        "RefSeq_id": refseq_ids,
        "Motif_under-representation_run1" : motif_underrep1s,
        "Motif_under-representation_run2" : motif_underrep2s,
        "Motif_under-representation_run3" : motif_underrep3s,
        "Mutational_susceptibility_run1" : mut_suscept1s,
        "Mutational_susceptibility_run2" : mut_suscept2s,
        "Mutational_susceptibility_run3" : mut_suscept3s,
        "Motif_under-representation" : motif_underrep_means,
        "Mutational_susceptibility" : mut_suscept_means,
        "Motif_under-representation_std" : motif_underrep_stds,
        "Mutational_susceptibility_std" : mut_suscept_stds,
        #"Gene_name": transcript_names,
        #"Transcript_variant": transcript_variants,
        "Gene_type": gene_types
    }
)
    
df.to_csv(Path(PurePath(base_path, f"cdur.pteropus.alecto.bytom.csv")), index=False)

