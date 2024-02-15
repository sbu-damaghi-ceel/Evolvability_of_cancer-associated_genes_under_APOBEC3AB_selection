# import package
library(ggplot2)
library(ggExtra)
library(extrafont)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"
cdur_file <- paste(base_path, "cdur.loxodonta_africana_pc_transcripts.tc.csv", sep = "")
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

out_file <- paste(base_path, "cdur_plot_elephant_refseq_transcripts_tc.tiff", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 600, width = 134, height = 134, device='tiff', units = "mm")

out_file <- paste(base_path, "cdur_plot_elephant_refseq_transcripts_tc.png", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 1200, width = 85, height = 85, units = "mm")

out_file <- paste(base_path, "cdur_plot_elephant_refseq_transcripts_small_tc.png", sep="")
ggsave(out_file, plot = cdur_plot, dpi = 1200, width = 65, height = 65, units = "mm")

