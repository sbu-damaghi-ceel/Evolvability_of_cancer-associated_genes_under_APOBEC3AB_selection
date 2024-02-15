import json
from pathlib import Path, PurePath

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

json_file = Path(PurePath(base_path, "all_lafricana.json"))

with open(json_file, "r") as f:
    data = json.load(f)

#print(data["genes"][52]["gene"].keys())
#print(data["genes"][52]["gene"]['transcripts'])

#ensembl_transcripts = []
#for i in range(len(data["genes"][52]["gene"]['transcripts'])):
#    ensembl_transcripts.append(data["genes"][52]["gene"]['transcripts'][i]['ensembl_transcript'])

#print(ensembl_transcripts)
output_file = Path(PurePath(base_path, "ncbi_lafricana.tsv"))
with open(output_file, "w") as o:
    o.write(f"gene_id\tensembl_gene_ids\tensembl_transcripts_ids\tsymbol\tsynonyms\ttax_id\ttaxname\tchromosomes\tcommon_name\tdescription\torientation\ttype\n")
    for i in range(len(data["genes"])):
        gene_id = data["genes"][i]["gene"]["gene_id"]
        if "ensembl_gene_ids" in data["genes"][i]["gene"].keys():
            ensembl_gene_ids = data["genes"][i]["gene"]["ensembl_gene_ids"]
        else:
            ensembl_gene_ids = "na"

        symbol = data["genes"][i]["gene"]["symbol"]
        tax_id = data["genes"][i]["gene"]["tax_id"]
        taxname = data["genes"][i]["gene"]["taxname"]
        if "chromosomes" in data["genes"][i]["gene"].keys():
            chromosomes = data["genes"][i]["gene"]["chromosomes"]
        else:
            if "chromosome" in data["genes"][i]["gene"].keys():
                chromosomes = data["genes"][i]["gene"]["chromosome"]
            else:
                chromosomes = "na"

        common_name = data["genes"][i]["gene"]["common_name"]
        description = data["genes"][i]["gene"]["description"]
        if "orientation" in data["genes"][i]["gene"].keys():
            orientation = data["genes"][i]["gene"]["orientation"]
        else:
            orientation = "na"

        if "type" in data["genes"][i]["gene"].keys():
            gene_type = data["genes"][i]["gene"]["type"] 
        else:
            gene_type = "na"      

        if "synonyms" in data["genes"][i]["gene"].keys():
            synonyms = data["genes"][i]["gene"]["synonyms"] 
        else:
            synonyms = "na" 


        if "transcripts" in data["genes"][i]["gene"].keys():
            ensembl_transcripts = []
            for j in range(len(data["genes"][i]["gene"]["transcripts"])):
                if "ensembl_transcript" in data["genes"][i]["gene"]["transcripts"][j].keys(): 
                    ensembl_transcripts.append(data["genes"][i]["gene"]["transcripts"][j]["ensembl_transcript"])
                else:
                    ensembl_transcripts.append("na")
        else:
            ensembl_transcripts = "na"
        

        o.write(f"{gene_id}\t{ensembl_gene_ids}\t{ensembl_transcripts}\t{symbol}\t{synonyms}\t{tax_id}\t{taxname}\t{chromosomes}\t{common_name}\t{description}\t{orientation}\t{gene_type}\n")
