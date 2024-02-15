# import package
library(dplyr)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/14_additional_work_5/" # CHANGE HERE
cdur_file <- paste(base_path, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.csv", sep = "")
cdur_data <- read.table(file = cdur_file,  sep = ',', header = TRUE)

# for gene ontology analysis
get_genename <- function(str) { 
  if (substring(str, 1, 4) == "ENST") {
    return(str)
  } else {
    return (substring(str, 1, nchar(str)-4))  
  }
  
}

x <- runif(10)
y <- runif(10)

for (i in 1:10) {
  rand_box <- cdur_data %>% 
    filter(Motif_under.representation >= x[i]-0.025 & Motif_under.representation <= x[i]+0.025) %>% 
    filter(Mutational_susceptibility >= y[i]-0.025 & Mutational_susceptibility <= y[i]+0.025)
  
  rand_box_file <- paste(base_path, "rand_box_", as.character(x[i]), "_", as.character(y[i]), ".txt", sep="")
  write(unique(unlist(lapply(rand_box$Transcript_name, get_genename))), rand_box_file)
  
}

# filter left-top
left_top <- cdur_data %>% 
  filter(Motif_under.representation <= 0.05) %>% 
  filter(Mutational_susceptibility <= 0.05)

left_top_file <- paste(base_path, "top_left.txt", sep="")
write(unique(unlist(lapply(left_top$Transcript_name, get_genename))), left_top_file)

length(left_top$Transcript_name)

outside_box <- cdur_data %>% 
  filter(! Transcript_name %in% left_top$Transcript_name) 

gene_outside <- unique(unlist(lapply(outside_box$Transcript_name, get_genename)))

sample(gene_outside, length(left_top$Transcript_name))

for (i in 11:100) {
  random_n_file <- paste(base_path, "random_n_", as.character(i), ".txt", sep="")
  write(sample(gene_outside, length(left_top$Transcript_name)), random_n_file)
  
}

