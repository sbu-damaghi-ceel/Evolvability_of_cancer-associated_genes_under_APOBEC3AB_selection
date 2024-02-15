library(biomaRt)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"

# get biomart from GRCh38 and 37
mart37 <- useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl", GRCh=37)
mart38 <- useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")

# import preprocessed COSMIC point mutation and NCBI genes
cosmic_mutation_file <- paste(base_path, "CosmicMutantExport_genename_id_sorted_uniq.tsv", sep = "")   
genename_id <- read.csv(cosmic_mutation_file, header = FALSE, sep = "\t")

ncbi_human_file <- paste(base_path, "ncbi_human_clean.tsv", sep = "")
ncbi_human <- read.csv(ncbi_human_file, header = TRUE, sep = "\t")

# extract gene names and ids
genenames <- as.vector(genename_id$V1)
id_versions <- as.vector(genename_id$V2)

ids <- c()
for (i in 1:length(id_versions)){
  ids <- c(ids, strsplit(id_versions[i], "\\.")[[1]][1])
}

# set attributes to check
attributes.of.interest = c(
  "external_gene_name",
  "ensembl_transcript_id_version",
  "ensembl_transcript_id",
  "ensembl_gene_id",
  "chromosome_name",
  "start_position",
  "end_position"
)

# get GRCh38 BM result
rlt <- getBM(
  attributes = attributes.of.interest,
  filters = "ensembl_transcript_id",
  values = ids,
  mart = mart38
)

idxs_in38 <- c()
idxs_notin38 <- c()
ensg <- c()
for (i in 1:length(ids)) {
  if (ids[i] %in% rlt$ensembl_transcript_id) {
    idxs_in38 <- c(idxs_in38, i)
    ensg <- c(ensg, rlt$ensembl_gene_id[rlt$ensembl_transcript_id == ids[i]])
  } else {
    idxs_notin38 <- c(idxs_notin38, i)
  }
}

out_table <- data.frame(
  gene_name = genenames[idxs_in38],
  ensembl_transcript_id_version = id_versions[idxs_in38],
  ensembl_transcript_id = ids[idxs_in38],
  ensembl_gene_id = ensg
)

write.table(
  out_table,
  file = paste(base_path, "CosmicMutantExport_detected_by_biomaRt_38.tsv", sep = ""),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ensg38 = ensg

# try GRCH37

rlt37 <- getBM(
  attributes = attributes.of.interest,
  filters = "ensembl_transcript_id",
  values = ids[idxs_notin38],
  mart = mart37
)

length(rlt37$external_gene_name)

ids37 = ids[idxs_notin38]

idxs_in37 <- c()
idxs_notin37 <- c()
enst_in37 <- c()
ensg_in37 <- c()
for (i in 1:length(ids37)) {
  if (ids37[i] %in% rlt37$ensembl_transcript_id) {
    idxs_in37 <- c(idxs_in37, i)
    enst_in37 <- c(enst_in37, rlt37$ensembl_transcript_id_version[rlt37$ensembl_transcript_id == ids37[i]])
    ensg_in37 <- c(ensg_in37, rlt37$ensembl_gene_id[rlt37$ensembl_transcript_id == ids37[i]])
  } else {
    idxs_notin37 <- c(idxs_notin37, i)
  }
}

out_table37 <- data.frame(
  gene_name = genenames[idxs_in37],
  ensembl_transcript_id_version = enst_in37,
  ensembl_transcript_id = ids37[idxs_in37],
  ensembl_gene_id = ensg_in37
)

write.table(
  out_table37,
  file = paste(base_path, "CosmicMutantExport_detected_by_biomaRt_37.tsv", sep = ""),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ids_left = ids37[idxs_notin37]

all_gene_names_detected = c(out_table$gene_name, out_table37$gene_name)

gene_names_not_detected <- c()
id_version_not_detected <- c()
id_not_detected <- c()
for (i in 1:length(ids)) {
  if (ids[i] %in% ids_left){
    gene_names_not_detected <- c(gene_names_not_detected, genenames[i])
    id_version_not_detected <- c(id_version_not_detected, id_versions[i])
    id_not_detected <- c(id_not_detected, ids[i])
  }
}
length(gene_names_not_detected)

gene_only <- c()
idx_gene_only <- c()
for (i in 1:length(gene_names_not_detected)){
  if (!(gene_names_not_detected[i] %in% all_gene_names_detected)) {
    gene_only <- c(gene_only, gene_names_not_detected[i])
    idx_gene_only <- c(idx_gene_only, i)
  }
}

length(gene_only)

out_table_really_not_detected <- data.frame(
  gene_name = gene_names_not_detected[idx_gene_only],
  ensembl_transcript_id_version = id_version_not_detected[idx_gene_only],
  ensembl_transcript_id = id_not_detected[idx_gene_only]
)

ids[idxs]

write.table(
  out_table_really_not_detected,
  file = paste(base_path, "CosmicMutantExport_not_detected_in_anyway.tsv", sep = ""),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ensembl_gene_ids <- c(out_table$ensembl_gene_id, out_table37$ensembl_gene_id)
ensembl_gene_ids_unique <- unique(ensembl_gene_ids)
idx_detected_in_ncbi <- c()
idx_not_detected_in_ncbi <- c()
ids_got <- c()
for (i in 1:length(ncbi_human$gene_id)) {
  if (grepl(",",ncbi_human$ensembl_gene_ids[i])) {
    to_check = unlist(strsplit(ncbi_human$ensembl_gene_ids[i],","))
    
    flag = 0
    for (id in to_check){
      if (id %in% ensembl_gene_ids_unique){
        flag = flag + 1
      }
    }
    
    if (flag != 0) {
      idx_detected_in_ncbi <- c(idx_detected_in_ncbi, i)
      ids_got <- c(ids_got, ncbi_human$ensembl_gene_ids[i])
    } else {
      idx_not_detected_in_ncbi <- c(idx_not_detected_in_ncbi, i)
    }
  } else {
    if (ncbi_human$ensembl_gene_ids[i] %in% ensembl_gene_ids_unique) {
      idx_detected_in_ncbi <- c(idx_detected_in_ncbi, i)
      ids_got <- c(ids_got, ncbi_human$ensembl_gene_ids[i])
    } else {
      idx_not_detected_in_ncbi <- c(idx_not_detected_in_ncbi, i)
    }
  }
}

length(ensembl_gene_ids_unique)
length(idx_detected_in_ncbi)
length(idx_not_detected_in_ncbi)
length(ids_got)


write.table(
  ncbi_human[idx_not_detected_in_ncbi,],
  file = paste(base_path, "ncbi_human_clean_no_mutation.tsv", sep = ""),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
