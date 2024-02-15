# import package
library(dplyr)
library(ggplot2)
library(ggExtra)
library(lubridate)
library(extrafont)
library(gridExtra)
library(egg)
library(cowplot)
library(ggpubr)

# file directories
base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"

data_file <- paste(base_path, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.tc.geneclass.csv", sep = "")
df.data <- read.table(file = data_file,  sep = '\t', header = TRUE)

df.data <- df.data %>% filter(gene_class != "regular")

rep_cancer <- df.data %>% filter(gene_class == "cancer") %>% select(Motif_under.representation)
rep_cancer <- rep_cancer$Motif_under.representation
rep_control <- df.data %>% filter(gene_class == "control") %>% select(Motif_under.representation)
rep_control <- rep_control$Motif_under.representation

ks.test(rep_cancer, rep_control)
wilcox.test(rep_cancer, rep_control)

sus_cancer <- df.data %>% filter(gene_class == "cancer") %>% select(Mutational_susceptibility)
sus_cancer <- sus_cancer$Mutational_susceptibility
sus_control <- df.data %>% filter(gene_class == "control") %>% select(Mutational_susceptibility)
sus_control <- sus_control$Mutational_susceptibility

ks.test(sus_cancer, sus_control)
wilcox.test(sus_cancer, sus_control)

comparison.plot <- ggplot(df.data, aes(x=Motif_under.representation, y=Mutational_susceptibility)) + 
  geom_point(data=df.data,aes(color=factor(gene_class), size=factor(gene_class), alpha=factor(gene_class)), stroke = 0.01) +
  geom_hline(yintercept = 0.05, size=0.2) +
  geom_hline(yintercept = 0.95, size=0.2) +
  geom_vline(xintercept = 0.05, size=0.2) +
  geom_vline(xintercept = 0.95, size=0.2)


comparison.plot <- comparison.plot + 
  scale_y_reverse() +
  labs(
    x = "Motif under-representation (p-value)",
    y = "Mutational susceptibility (p-value)"
  ) +
  theme_light() +
  scale_color_manual(name = "Type",
                     labels = c("Cancer", "Control"),
                     values = c(
                       "red","blue"
                     )) +
  scale_size_manual(name = "Type",
                    labels = c("Cancer", "Control"),
    values = c(0.5,0.5)
  ) +
  scale_alpha_manual(name = "Type",
                     labels = c("Cancer", "Control"),
    values = c(1,1)
  ) +
  theme(
    legend.position = 'bottom',
    axis.title.x = element_text(size = 8, family = "Arial"), 
    axis.title.y = element_text(size = 8, family = "Arial"), 
    axis.text.x = element_text(size = 6, family = "Arial"), 
    axis.text.y = element_text(size = 6, family = "Arial"), 
    legend.text = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    panel.grid.minor = element_blank()
    ) +
  guides(color = guide_legend(override.aes = list(size = 3)))


comparison.plot <- ggMarginal(
  comparison.plot,                  
  type="density", 
  size = 5,
  alpha = 0.1,
  groupColour = TRUE, groupFill = TRUE,
  linewidth = 0.1
)

out_file <- paste(base_path, "cancer_vs_control.png", sep="")
ggsave(out_file, plot = comparison.plot, dpi = 1200, width = 85, height = 97, units = "mm")

