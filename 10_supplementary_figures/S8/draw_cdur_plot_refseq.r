# import package
library(ggplot2)
library(ggExtra)
library(extrafont)
library(dplyr)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/13_additional_work_4/" # CHANGE HERE
cdur_file <- paste(base_path, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.csv", sep = "")
cdur_data <- read.table(file = cdur_file,  sep = ',', header = TRUE)

# grep "chrX" gencode.v40.annotation.gtf | awk -F"\t" '{if ($3 == "transcript") print $0}' | awk -F"\t" '{print $9}' | awk -F'\"' 'BEGIN{OFS=","}{print $2,$4}' > gencode.v40.chrX.ensembl_id.csv
chrX_file <- paste(base_path, "gencode.v40.chrX.ensembl_id.csv", sep = "")
chrX <- read.table(file = chrX_file,  sep = ',', header = FALSE, col.names = c("Ensembl_gene_id", "Ensembl_transcript_id"))
chrY_file <- paste(base_path, "gencode.v40.chrY.ensembl_id.csv", sep = "")
chrY <- read.table(file = chrY_file,  sep = ',', header = FALSE, col.names = c("Ensembl_gene_id", "Ensembl_transcript_id"))

mark_XY <- function(gene){
  if (gene %in% chrX$Ensembl_transcript_id) { return("chrX")}
  else if (gene %in% chrY$Ensembl_transcript_id) { return("chrY")}
  else { return("auto") }
}

cdur_data$chr_type <- sapply(cdur_data$Ensembl_transcript_id, mark_XY)

# statistical analysis
repx <- cdur_data %>% filter(chr_type == "chrX") %>% select(Motif_under.representation)
repx <- repx$Motif_under.representation
repy <- cdur_data %>% filter(chr_type == "chrY") %>% select(Motif_under.representation)
repy <- repy$Motif_under.representation

ks.test(repx, repy)
wilcox.test(repx, repy)

susx <- cdur_data %>% filter(chr_type == "chrX") %>% select(Mutational_susceptibility)
susx <- susx$Mutational_susceptibility
susy <- cdur_data %>% filter(chr_type == "chrY") %>% select(Mutational_susceptibility)
susy <- susy$Mutational_susceptibility

ks.test(susx, susy)
wilcox.test(susx, susy)





cdur_data %>% group_by(chr_type) %>% summarise(n=n())

cdur_plot <- ggplot(cdur_data, aes(x=Motif_under.representation, y=Mutational_susceptibility)) + 
  geom_point(aes(color=factor(chr_type, levels = c("chrX", "chrY", "auto"))), size=0.4, alpha=0.6, stroke = 0.01) +
  geom_hline(yintercept = 0.05, linewidth=0.2) +
  geom_hline(yintercept = 0.95, linewidth=0.2) +
  geom_vline(xintercept = 0.05, linewidth=0.2) +
  geom_vline(xintercept = 0.95, linewidth=0.2)

cdur_plot <- cdur_plot + 
  scale_y_reverse() +
  labs(
    x= "Motif under-representation (p-value)",
    y= "Mutational susceptibility (p-value)"
  ) +
  theme_light() +
  scale_color_manual(
    name = "Chromosome",
    labels = c("X", "Y", "Autosomal"),
    values = c("#001f54", "#ea3546","grey")
    ) +
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(
    legend.position="bottom",
    axis.title.x = element_text(size = 8, family = "Arial"), 
    axis.title.y = element_text(size = 8, family = "Arial"),
    axis.text.x = element_text(size = 6, family = "Arial"),
    axis.text.y = element_text(size = 6, family = "Arial"),
    panel.grid.minor = element_blank()
  )

cdur_plot <- ggMarginal(
  cdur_plot,                  
  type="density", 
  size = 5,
  alpha = 0.1,
  groupColour = TRUE, groupFill = TRUE,
  linewidth = 0.1
)
cdur_plot

out_file <- paste(base_path, "cdur_plot_refseq_transcripts_tc_XY.png", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 1200, width = 85, height = 95, units = "mm")

out_file <- paste(base_path, "cdur_plot_refseq_transcripts_small_tc_XY.png", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 1200, width = 65, height = 75, units = "mm")

