library(biomaRt)
library(dplyr)

# use GRCh38 to get information
mart38 <- useEnsembl(biomart="ensembl", dataset = "mmusculus_gene_ensembl")

# set attributes to view
# you can see possible attributes using: listAttributes(mart38)
attributes.of.interest = c(
  "external_gene_name",
  "ensembl_transcript_id_version",
  "ensembl_transcript_id",
  "ensembl_gene_id",
  "description",
  "gene_biotype",
  "transcript_gencode_basic",
  "transcript_appris",
  "transcript_is_canonical",
  "refseq_mrna"
)

# import name and id of genes from gencode (cleaned for valid sequence)
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"
name_id_file <- paste(base_path, "gencode.vM33.pc_transcripts.nopary.cdsonly.name.id.csv", sep="")
genename_id <- read.csv(name_id_file, sep = ",")

# retrieve information from biomart
rlt <- getBM(
  attributes = attributes.of.interest,
  filters = "ensembl_transcript_id",
  values = genename_id$Ensembl_transcript_id,
  mart = mart38
)

# get only distinct gene id
rlt_distinct <- rlt %>% distinct(ensembl_transcript_id, .keep_all = TRUE)

# get transcripts marked as canonical
rlt_canonical <- rlt_distinct %>% filter(transcript_is_canonical == 1)
write.csv2(rlt_canonical$ensembl_transcript_id, paste(base_path, "gencode.vM33.pc_transcripts.nopary.cdsonly.canonical.id.csv", sep=""), row.names = FALSE, quote=FALSE)

# get transcript having refseq id
rlt_refseq <- rlt_distinct %>% filter(refseq_mrna != "")

# get ids not found from GRCh38
id_not_found <- setdiff(genename_id$Ensembl_transcript_id, rlt_distinct$ensembl_transcript_id)

# This gave 8 records
# "ENSMUSG00000039337" "ENSMUSG00000039329" "ENSMUSG00000023083" "ENSMUSG00000024448" 
# Tex19.2              Tex19.1              H2-M10.2             H2-M10.1                      
# NM_027622.3          NM_028602.3          NM_177923.1          NM_001360498.1,NM_013544.3
#
# "ENSMUSG00000058124" "ENSMUSG00000048231" "ENSMUSG00000037246" "ENSMUSG00000037130"
# H2-M10.3             H2-M10.4             H2-M10.5             H2-M10.6        
# NM_201608.2          NM_177634.2          NM_177637.3          NM_201611.2

manual_add <- c("ENSMUSG00000039337", "ENSMUSG00000039329", "ENSMUSG00000023083", "ENSMUSG00000024448", "ENSMUSG00000058124", "ENSMUSG00000048231", "ENSMUSG00000037246", "ENSMUSG00000037130")

# combine transcripts having refseq id
rlt_refseq_all <- c(rlt_refseq$ensembl_transcript_id, manual_add)
write.csv(rlt_refseq_all, paste(base_path, "gencode.vM33.pc_transcripts.nopary.cdsonly.haverefseq.id.csv", sep=""), row.names = FALSE, quote=FALSE)

