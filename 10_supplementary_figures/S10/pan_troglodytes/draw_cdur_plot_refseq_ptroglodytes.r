# import package
library(ggplot2)
library(ggExtra)
library(extrafont)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"
cdur_file <- paste(base_path, "cdur.pteropus_alecto_pc_transcripts.tyc.csv", sep = "")
cdur_data <- read.table(file = cdur_file,  sep = ',', header = TRUE)

cdur_plot <- ggplot(cdur_data, aes(x=Motif_under.representation, y=Mutational_susceptibility)) + 
  geom_point(size=0.4, alpha=0.25, fill = "black", stroke = 0.01) +
  geom_hline(yintercept = 0.05, linewidth=0.2) +
  geom_hline(yintercept = 0.95, linewidth=0.2) +
  geom_vline(xintercept = 0.05, linewidth=0.2) +
  geom_vline(xintercept = 0.95, linewidth=0.2)

cdur_plot <- cdur_plot + 
  scale_y_reverse() +
  labs(
    x = "Motif under-representation (p-value)",
    y = "Mutational susceptibility (p-value)"
  ) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    panel.grid.minor = element_blank()
  )

cdur_plot <- ggMarginal(
  cdur_plot,                  
  type="histogram", 
  size = 10,
  xparams = list(binwidth = 0.01, linewidth = 0.1), 
  yparams = list(binwidth = 0.01, linewidth = 0.1)
)

out_file <- paste(base_path, "cdur_plot_palecto_refseq_transcripts_tyc.tiff", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 600, width = 134, height = 134, device='tiff', units = "mm")

out_file <- paste(base_path, "cdur_plot_palecto_refseq_transcripts_tyc.png", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 1200, width = 85, height = 85, units = "mm")

out_file <- paste(base_path, "cdur_plot_palecto_refseq_transcripts_small_tyc.png", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 1200, width = 65, height = 65, units = "mm")

#
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"
bat_file <- paste(base_path, "cdur.pteropus_alecto_pc_transcripts.tc.csv", sep = "")
bat_data <- read.table(file = bat_file,  sep = ',', header = TRUE)
human_file <- paste(base_path, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.csv", sep = "")
human_data <- read.table(file = human_file,  sep = ',', header = TRUE)
bat_data

rep_bat <- bat_data$Motif_under.representation
rep_human <- human_data$Motif_under.representation

rep_data <- data.frame(
  species = c(rep("human", length(rep_human)),rep("bat", length(rep_bat))),
  Motif_under.representation = c(rep_human, rep_bat)
)

mur_plot <- ggplot(rep_data, aes(x=Motif_under.representation, fill=species)) + 
  geom_density(alpha=0.25, linewidth = 0.2) +
  guides(fill = guide_legend(title = "Species", order = 1)) +
  xlab("Motif under-representation (p-value)") +
  ylab("Density") +
  theme_light() +
  theme(
    legend.position = 'right',
    #plot.title = element_text(size = 16, family = "Arial", hjust=0.5),
    axis.title.x = element_text(size = 8, family = "Arial"), #24 7
    axis.title.y = element_text(size = 8, family = "Arial"), #24 7
    axis.text.x = element_text(size = 6, family = "Arial"), # 18 5
    axis.text.y = element_text(size = 6, family = "Arial"), # 18 5
    legend.text = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 7, family = "Arial"),
    legend.key.size = unit(2,"mm")
  )

out_file <- paste(base_path, "human_vs_bat_motif_under-representation.png", sep="")
ggsave(out_file, plot = mur_plot, dpi = 1200, width = 115, height = 65, units = "mm")



ks.test(rep_bat, rep_human, exact = TRUE)
wilcox.test(rep_bat, rep_human)

sus_bat <- bat_data$Mutational_susceptibility
sus_human <- human_data$Mutational_susceptibility

ks.test(sus_bat, sus_human, exact = TRUE)
wilcox.test(sus_bat, sus_human)

sus_data <- data.frame(
  species = c(rep("human", length(sus_human)),rep("bat", length(sus_bat))),
  Mutational_susceptibility = c(sus_human, sus_bat)
)

sus_plot <- ggplot(sus_data, aes(y=Mutational_susceptibility, fill=species)) + 
  geom_density(alpha=0.25, linewidth = 0.2) +
  guides(fill = guide_legend(title = "Species", order = 1)) +
  ylab("Mutational susceptibility (p-value)") +
  xlab("Density") +
  theme_light() +
  scale_y_reverse() +
  theme(
    legend.position = 'right',
    #plot.title = element_text(size = 16, family = "Arial", hjust=0.5),
    axis.title.x = element_text(size = 8, family = "Arial"), #24 7
    axis.title.y = element_text(size = 8, family = "Arial"), #24 7
    axis.text.x = element_text(size = 6, family = "Arial"), # 18 5
    axis.text.y = element_text(size = 6, family = "Arial"), # 18 5
    legend.text = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 7, family = "Arial"),
    legend.box.margin = margin(0,0,0,-10),
    legend.key.size = unit(2,"mm")
  )

out_file <- paste(base_path, "human_vs_bat_mutational_susceptibility.png", sep="")
ggsave(out_file, plot = sus_plot, dpi = 1200, width = 60, height = 65, units = "mm")

ggplot(sus_data, aes(y=Mutational_susceptibility, fill=species)) + geom_density(alpha=0.25) + scale_y_reverse() 

mean(rep_bat)
mean(rep_human)
mean(sus_bat)
mean(sus_human)

library("sm")
sm.density.compare(rep_bat, rep_human, model="equal", nboot= 500, ngrid= 100)
sm.density.compare(sus_bat, sus_human, model="equal", nboot= 500, ngrid= 100)

