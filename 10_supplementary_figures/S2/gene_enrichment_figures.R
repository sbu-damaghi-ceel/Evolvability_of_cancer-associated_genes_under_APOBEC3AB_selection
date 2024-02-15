library(ggplot2)
library(dplyr)
library(extrafont)
library(egg)
library(stringr) 

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/14_additional_work_5/"

# https://www.genome.jp/kegg/pathway.html
cancer_pathways <- c("Pathways in cancer", 
                     "Transcriptional misregulation in cancer", 
                     "MicroRNAs in cancer", 
                     "Proteoglycans in cancer",
                     "Chemical carcinogenesis",
                     "Chemical carcinogenesis - DNA adducts",
                     "Chemical carcinogenesis - receptor activation",
                     "Chemical carcinogenesis - reactive oxygen species",
                     "Viral carcinogenesis",
                     "Central carbon metabolism in cancer",
                     "Choline metabolism in cancer",
                     "PD-L1 expression and PD-1 checkpoint pathway in cancer",
                     "Colorectal cancer",
                     "Pancreatic cancer",
                     "Hepatocellular carcinoma",
                     "Gastric cancer",
                     "Glioma",
                     "Thyroid cancer",
                     "Acute myeloid leukemia",
                     "Chronic myeloid leukemia",
                     "Basal cell carcinoma",
                     "Melanoma",
                     "Renal cell carcinoma",
                     "Bladder cancer",
                     "Prostate cancer",
                     "Endometrial cancer",
                     "Breast cancer",
                     "Small cell lung cancer",
                     "Non-small cell lung cancer"
)

cnt <- c()
avgs <- c()
samples <- c()
fold_enrichs <- c()
for (i in 1:20){
  kegg_file <- paste(base_path, "enrichment_", as.character(i), ".csv", sep="")
  kegg_data <- read.table(file = kegg_file,  sep = ',', header = TRUE)
  
  tmp <- kegg_data %>% filter(Pathway %in% cancer_pathways)
  
  samples <- c(samples, rep(paste("rand_", as.character(i), sep = ""), length(tmp$Pathway)))
  fold_enrichs <- c(fold_enrichs, tmp$Fold.Enrichment)
  cnt <- c(cnt, length(tmp$Pathway))
  avgs <- c(avgs, mean(tmp$Fold.Enrichment))
  
}
mean(cnt)
mean(na.omit(avgs))

topleft_file <- paste(base_path, "enrichment_top_left.csv", sep="")
topleft_data <- read.table(file = topleft_file,  sep = ',', header = TRUE)

topleft <- topleft_data %>% filter(Pathway %in% cancer_pathways)
length(topleft$Pathway)
mean(topleft$Fold.Enrichment)
min(topleft$Fold.Enrichment)

samples <- c(samples, rep(paste("Top-left"), length(topleft$Pathway)))
fold_enrichs <- c(fold_enrichs, topleft$Fold.Enrichment)

samples <- c(samples, "rand_2", "rand_3", "rand_4", "rand_5", "rand_8", "rand_10", 
                      "rand_13","rand_16","rand_17","rand_19", "rand_20")
fold_enrichs <- c(fold_enrichs, rep(1, 11))

plot.df <- data.frame(
  sample <- samples,
  Fold.Enrichment <- fold_enrichs
)

plot1 <- ggplot(plot.df, aes(x = factor(sample, levels = c("Top-left", "rand_1", "rand_2", "rand_3",
                                                 "rand_4",   "rand_5", "rand_6", "rand_7",
                                                 "rand_8",   "rand_9", "rand_10","rand_11",
                                                 "rand_12",  "rand_13","rand_14","rand_15",
                                                 "rand_16",  "rand_17","rand_18","rand_19",
                                                 "rand_20")),
                    y = Fold.Enrichment,
                    color = factor(sample, levels = c("Top-left", "rand_1", "rand_2", "rand_3",
                                                      "rand_4",   "rand_5", "rand_6", "rand_7",
                                                      "rand_8",   "rand_9", "rand_10","rand_11",
                                                      "rand_12",  "rand_13","rand_14","rand_15",
                                                      "rand_16",  "rand_17","rand_18","rand_19",
                                                      "rand_20")))) +
  geom_boxplot() +
  geom_point() +
  ylim(1, 3) +
  xlab("Sampling") +
  scale_x_discrete(drop = FALSE) +
  scale_color_discrete(name = "Sampling") +
  guides(color=guide_legend(ncol=2)) +
  theme_light() +
  theme(legend.position = 'right',
        #plot.title = element_text(size = 16, family = "Arial", hjust=0.5),
        axis.title.x = element_text(size = 8, family = "Arial"), #24 7
        axis.title.y = element_text(size = 8, family = "Arial"), #24 7
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6, family = "Arial"), # 18 5
        axis.text.y = element_text(size = 6, family = "Arial"), # 18 5
        legend.text = element_text(size = 6, family = "Arial"),
        legend.title = element_text(size = 7, family = "Arial"))
plot1
out_file <- paste(base_path, "compare_random_selection.png", sep="")
ggsave(out_file, plot = plot1, dpi = 1200, width = 150, height = 100, units = "mm")
