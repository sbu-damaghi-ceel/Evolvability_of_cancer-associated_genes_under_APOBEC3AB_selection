library(ggplot2)
library(dplyr)
library(extrafont)
library(egg)
library(stringr)
library(ggbeeswarm)
library(ggpubr)
library(ggsignif)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
orthologs_file <- paste(base_path, "all_control_genes_orthologs_cdur.csv", sep = "")
orthologs_data <- read.table(file = orthologs_file,  sep = ',', header = TRUE)
trajs_human_file <- paste(base_path, "all_trajs_all_control_genes_human_cdur.csv", sep = "")
trajs_human_data <- read.table(file = trajs_human_file,  sep = ',', header = TRUE)

orthologs_genenames <- orthologs_data$Gene_name %>% unique()
human_trajs_genenames <- trajs_human_data$Gene_name %>% unique()

stdmmrs <- c()
stdmss <- c()
for (gene in human_trajs_genenames) {
  stds <- trajs_human_data %>% 
    filter(Gene_name == gene) %>% 
    group_by(Trajectory) %>%
    filter(Traj_points == max(Traj_points)) %>%
    arrange(Trajectory) %>%
    ungroup() %>%
    summarise(stdmmr = sd(Motif_under.representation), stdms = sd(Mutational_susceptibility))
  
  stdmmrs <- c(stdmmrs, stds$stdmmr)
  stdmss <- c(stdmss, stds$stdms)
}

std_of_human_trajs <- data.frame(
  Gene_name = human_trajs_genenames,
  Std_motif_underrep = stdmmrs,
  Std_mut_suscept = stdmss
)


std_of_ortholgos <- orthologs_data %>%
  group_by(Gene_name) %>%
  summarise(stdmmr = sd(Motif_under.representation), stdms = sd(Mutational_susceptibility))

df_plot <- data.frame(
  Gene_name = human_trajs_genenames,
  Std_of = c(rep("human_trajs", length(human_trajs_genenames)), rep("orthologs", length(human_trajs_genenames))),
  Std_motif_underrep = c(stdmmrs, std_of_ortholgos$stdmmr),
  Std_mut_suscept = c(stdmss, std_of_ortholgos$stdms)

)

df_plot

p_motif_underrep <- t.test(
  df_plot[df_plot$Std_of == "human_trajs", "Std_motif_underrep"],
  df_plot[df_plot$Std_of == "orthologs", "Std_motif_underrep"],
  paired = TRUE
)$p.value

motifrep_violin <- ggplot(df_plot, aes(x=factor(Std_of), y=Std_motif_underrep, fill=Std_of)) +
  geom_violin(scale = "count", trim = FALSE, alpha = 0.4, linewidth=0.1) +
  geom_boxplot(width=0.1, color="#d62828", size= 1, alpha=0.4, linewidth=0.3,outlier.stroke = 0, outlier.size = 0.1) +
  geom_beeswarm(cex = 1, size = 0.5, stroke = 0.01) +
  geom_signif(
    annotation = formatC(p_motif_underrep),
    y_position = 0.55,
    xmin = 1, xmax = 2,
    textsize = 8/.pt
  )+
  labs(
    x="",
    y="Motif under-representation standard deviation",
    fill = "Type"
  ) +
  scale_fill_manual(values = c("#ff006e","#3a86ff"))+
  scale_x_discrete(labels = c("Sequential mutations", "Orthologs")) +
  theme_light() +
  theme(
    axis.title.x = element_blank(), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"),
    legend.title =  element_text(size = 7, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom",
    legend.key.size = unit(4,"mm"),
    legend.box.margin=margin(-10,0,0,0),
    plot.margin = unit(c(1,1,1,1), "mm")
  ) + guides(fill = "none") + ylim(-0.1, 0.56)
  
motifrep_violin
out_file <- paste(base_path, "control_genes_std_motifrep.png", sep = "")
ggsave(out_file, plot = motifrep_violin, dpi = 1200, width = 45, height = 80, units = "mm")


p_mut_suscept <- t.test(
  df_plot[df_plot$Std_of == "human_trajs", "Std_mut_suscept"],
  df_plot[df_plot$Std_of == "orthologs", "Std_mut_suscept"],
  paired = TRUE
)$p.value

mutsuscept_violin <- ggplot(df_plot, aes(x=factor(Std_of), y=Std_mut_suscept, fill=Std_of)) +
  geom_violin(scale = "count", trim = FALSE, alpha = 0.4, linewidth=0.1) +
  geom_boxplot(width=0.1, color="#d62828", size= 1, alpha=0.4, linewidth=0.3,outlier.stroke = 0, outlier.size = 0.1) +
  geom_beeswarm(cex = 1, size = 0.5, stroke = 0.01) +
  geom_signif(
    annotation = formatC(p_mut_suscept),
    y_position = 0.55,
    xmin = 1, xmax = 2,
    textsize = 8/.pt
  )+
  labs(
    x="",
    y="Mutational susceptibility standard deviation",
    fill = "Type"
  ) +
  scale_fill_manual(values = c("#ff006e","#3a86ff"))+
  scale_x_discrete(labels = c("Sequential mutations", "Orthologs")) +
  theme_light() +
  theme(
    axis.title.x = element_blank(), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"),
    legend.title =  element_text(size = 7, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom",
    legend.key.size = unit(4,"mm"),
    legend.box.margin=margin(-10,0,0,0),
    plot.margin = unit(c(1,1,1,1), "mm")
  ) + guides(fill = "none") + ylim(-0.1, 0.56)

mutsuscept_violin
out_file <- paste(base_path, "control_genes_std_mutsuscept.png", sep = "")
ggsave(out_file, plot = mutsuscept_violin, dpi = 1200, width = 45, height = 80, units = "mm")



