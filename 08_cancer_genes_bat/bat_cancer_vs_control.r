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

data_file <- paste(base_path, "cdur.pteropus_alecto_pc_transcripts.tc.geneclass.csv", sep = "")
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

traj.plot <- ggplot(df.data, aes(x=Motif_under.representation, y=Mutational_susceptibility)) + 
  geom_point(data=df.data,aes(color=factor(gene_class), size=factor(gene_class), alpha=factor(gene_class)), stroke = 0.01) +
  geom_hline(yintercept = 0.05, size=0.2) +
  geom_hline(yintercept = 0.95, size=0.2) +
  geom_vline(xintercept = 0.05, size=0.2) +
  geom_vline(xintercept = 0.95, size=0.2)


traj.plot <- traj.plot + 
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
    #plot.title = element_text(size = 16, family = "Arial", hjust=0.5),
    axis.title.x = element_text(size = 8, family = "Arial"), #24 7
    axis.title.y = element_text(size = 8, family = "Arial"), #24 7
    axis.text.x = element_text(size = 6, family = "Arial"), # 18 5
    axis.text.y = element_text(size = 6, family = "Arial"), # 18 5
    legend.text = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    ) +
  guides(color = guide_legend(override.aes = list(size = 3)))


traj.plot <- ggMarginal(
  traj.plot,                  
  type="density", 
  size = 5,
  alpha = 0.1,
  groupColour = TRUE, groupFill = TRUE,
  linewidth = 0.1
)


traj.plot

out_file <- paste(base_path, "bat_cancer_vs_control.png", sep="")
ggsave(out_file, plot = traj.plot, dpi = 1200, width = 85, height = 97, units = "mm")







################################################################################
################################################################################

sorted_genes_x <- df.data %>% filter(gene_class == "cancer") %>% arrange(belowT_C_) %>% mutate(Name = substring(Protein, 1, nchar(Protein)-4)) %>% select(Name) %>% distinct()
sorted_genes_x <- sorted_genes_x$Name

top10 <- sorted_genes_x[1:10]
fileConn <- file(paste(out_dir, "cgc_left10.txt", sep = "/"))
writeLines(top10, fileConn)
close(fileConn)

top50 <- sorted_genes_x[1:50]
fileConn <- file(paste(out_dir, "cgc_left50.txt", sep = "/"))
writeLines(top50, fileConn)
close(fileConn)

top100 <- sorted_genes_x[1:100]
fileConn <- file(paste(out_dir, "cgc_left100.txt", sep = "/"))
writeLines(top100, fileConn)
close(fileConn)

top200 <- sorted_genes_x[1:200]
fileConn <- file(paste(out_dir, "cgc_left200.txt", sep = "/"))
writeLines(top200, fileConn)
close(fileConn)

bottom10 <- sorted_genes_x[704:713]
fileConn <- file(paste(out_dir, "cgc_right10.txt", sep = "/"))
writeLines(bottom10, fileConn)
close(fileConn)

bottom50 <- sorted_genes_x[664:713]
fileConn <- file(paste(out_dir, "cgc_right50.txt", sep = "/"))
writeLines(bottom50, fileConn)
close(fileConn)

bottom100 <- sorted_genes_x[614:713]
fileConn <- file(paste(out_dir, "cgc_right100.txt", sep = "/"))
writeLines(bottom100, fileConn)
close(fileConn)

bottom200 <- sorted_genes_x[514:713]
fileConn <- file(paste(out_dir, "cgc_right200.txt", sep = "/"))
writeLines(bottom200, fileConn)
close(fileConn)

cgc_file <- paste(in_dir, "cancer_gene_census.csv", sep = "/")
cgc_data <- read.table(file = cgc_file,  sep = ',', header = TRUE)
cgc_data %>% filter(Gene.Symbol %in% top10) %>% select(c(Gene.Symbol, Tumour.Types.Somatic.))
cgc_data %>% filter(Gene.Symbol %in% top50) %>% select(c(Gene.Symbol, Tumour.Types.Somatic.))
cgc_data %>% filter(Gene.Symbol %in% top100) %>% select(c(Gene.Symbol, Tumour.Types.Somatic.))
cgc_data %>% filter(Gene.Symbol %in% top200) %>% select(c(Gene.Symbol, Tumour.Types.Somatic.))
