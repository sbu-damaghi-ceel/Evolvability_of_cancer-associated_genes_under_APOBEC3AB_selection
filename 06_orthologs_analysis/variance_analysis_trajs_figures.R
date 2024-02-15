library(ggplot2)
library(dplyr)
library(extrafont)
library(egg)
library(stringr)
library(ggbeeswarm)
library(ggpubr)
library(ggsignif)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE
orthologs_file <- paste(base_path, "all_cancer_genes_orthologs_cdur.csv", sep = "")
orthologs_data <- read.table(file = orthologs_file,  sep = ',', header = TRUE)
trajs_human_file <- paste(base_path, "all_trajs_all_cancer_genes_human_cdur.csv", sep = "")
trajs_human_data <- read.table(file = trajs_human_file,  sep = ',', header = TRUE)

orthologs_genenames <- orthologs_data$Gene_name %>% unique()
human_trajs_genenames <- trajs_human_data$Gene_name %>% unique()

gene_plot_df <- orthologs_data %>% filter(Gene_name == "BRCA1")

species_order = c(
  "Homo_sapiens", "Pan_troglodytes",
  "Canis_lupus", "Mus_musculus", 
  "Physeter_catodon", "Loxodonta_africana", 
  "Myotis_lucifugus", "Gallus_gallus", 
  "Xenopus_tropicalis", "Petromyzon_marinus"
)

num_trajs = 10

gene_plot <- ggplot(gene_plot_df, aes(x=Motif_under.representation, y=Mutational_susceptibility, color = factor(Species, levels = species_order))) +
  scale_y_reverse() +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.05, size=0.2) +
  geom_hline(yintercept = 0.95, size=0.2) +
  geom_vline(xintercept = 0.05, size=0.2) +
  geom_vline(xintercept = 0.95, size=0.2)

for (i in 1:num_trajs) {
  gene_traj <- trajs_human_data %>% 
    filter(Gene_name == "BRCA1") %>%
    filter(Trajectory == i) %>%
    filter(Traj_points > 0)
  
  ori_point <- gene_plot_df %>% filter(Species == "Homo_sapiens")
  
  traj_xs <-c(ori_point$Motif_under.representation, gene_traj$Motif_under.representation)
  traj_ys <-c(ori_point$Mutational_susceptibility, gene_traj$Mutational_susceptibility)
  
  gene_traj_df <- data.frame(
    Motif_under.representation = traj_xs,
    Mutational_susceptibility = traj_ys,
    Species = rep("Homo_sapiens", length(traj_xs))
  )
  
  last_point <- gene_traj %>% slice(n())
  last_point$Species <- "Homo_sapiens"
  
  gene_plot <- gene_plot +
    geom_line(data=gene_traj_df, 
              aes(x=Motif_under.representation, y=Mutational_susceptibility), 
              size = 0.5,
              alpha = 0.8) +
    geom_point(data=last_point,
               aes(x=Motif_under.representation, y=Mutational_susceptibility),
               size = 0.5,
               shape = 4,
               alpha = 0.8)
  
}

gene_plot <- gene_plot +
  scale_color_discrete(name = "Species",
                    labels =c(
                      "Homo sapiens", "Pan troglodytes",
                      "Canis lupus", "Mus musculus", 
                      "Physeter catodon", "Loxodonta africana", 
                      "Myotis lucifugus", "Gallus gallus", 
                      "Xenopus tropicalis", "Petromyzon marinus"
                    ))+
  ggtitle("BRCA1 Orthologs") +
  labs(
    x = "Motif under-representation (p-value)",
    y = "Mutational susceptibility (p-value)"
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 8, family = "Arial", hjust = 0.5),
    legend.title = element_text(size = 8, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(1,"mm"),
    legend.box.margin=margin(-10,0,0,-15),
    plot.margin = unit(c(1,1,1,1), "mm")
  )  +
  guides(color=guide_legend(nrow=4, byrow=TRUE))
gene_plot

plot.fixed <- set_panel_size(gene_plot,
                             width = unit(75, "mm"),
                             height = unit(75, "mm"))
  
out_file <- paste(base_path, "BRCA1.png", sep = "")
ggsave(out_file, plot = plot.fixed, dpi = 1200, width = 88, height = 110, units = "mm")



