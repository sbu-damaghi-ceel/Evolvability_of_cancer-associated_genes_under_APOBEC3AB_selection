library(ggplot2)
library(dplyr)
library(extrafont)
library(egg)
library(stringr) 

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/"
kegg_file <- paste(base_path, "left-top_kegg_enrichment.csv", sep = "")
kegg_data <- read.table(file = kegg_file,  sep = ',', header = TRUE)
kegg_data

mycolors <- c("#2d00f7", "#6a00f4", "#8900f2", "#a100f2", "#b100e8", "#bc00dd", "#d100d1", "#db00b6", "#e500a4", "#f20089")

enrichment_plot <- kegg_data %>%
  mutate(nlog10FDR = floor(-log10(Enrichment.FDR))) %>%
  mutate(nlog10FDR=factor(nlog10FDR, levels=c(4,5,6,7,8,9))) %>%
  arrange(Fold.Enrichment) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Pathway=factor(Pathway, levels=Pathway)) %>%   # This trick update the factor levels
  ggplot( aes(x=Pathway, y=Fold.Enrichment, color=nlog10FDR)) +
  geom_segment( aes(xend=Pathway, yend=0), size = 0.8) +
  geom_point( aes(size=nGenes, color=nlog10FDR)) +
  coord_flip() +
  labs(x="", y= "Fold Enrichment") +
  theme(
    axis.title.x = element_text(size = 7, family = "Arial"), #24
    axis.title.y = element_text(size = 7, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    legend.title =  element_text(size = 7, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.key.size = unit(2,"mm"),
    legend.box.margin=margin(0,0,0,-10),
    plot.margin = unit(c(1,1,1,1), "mm")
  ) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  scale_colour_manual(
    name = "-log10(FDR)",
    values = mycolors[c(1,2,4,6,8,10)],
    drop = FALSE
    
  ) + 
  guides(size=guide_legend(title="N. of Genes", order = 1), color= guide_legend(order=2)) + scale_size(range = c(0.3,2))

plot.fixed <- set_panel_size(enrichment_plot,
                             width = unit(40, "mm"),
                             height = unit(90, "mm"))

#enrichment_plot
#out_file <- paste(base_path, "left-top_kegg_enrichment.tiff", sep="")
#ggsave(out_file, plot = enrichment_plot, dpi = 600, width = 170, height = 100, device='tiff', units = "mm")

#out_file <- paste(base_path, "left-top_kegg_enrichment.png", sep="")
#ggsave(out_file, plot = enrichment_plot, dpi = 600, width = 100, height = 65, units = "mm")

out_file <- paste(base_path, "left-top_kegg_enrichment2.png", sep="")
ggsave(out_file, plot = plot.fixed, dpi = 1200, width = 90, height = 100, units = "mm")

################################################################################
gobp_file <- paste(base_path, "left-top_gobp_enrichment.csv", sep = "")
gobp_data <- read.table(file = gobp_file,  sep = ',', header = TRUE)
gobp_data

mycolors <- c("#2d00f7", "#6a00f4", "#8900f2", "#a100f2", "#b100e8", "#bc00dd", "#d100d1", "#db00b6", "#e500a4", "#f20089")

enrichment_plot <- gobp_data %>%
  mutate(nlog10FDR = 2*(floor(-log10(Enrichment.FDR)/2))) %>%
  mutate(nlog10FDR=factor(nlog10FDR, levels=c(16,18,20,22,24,26,28))) %>%
  arrange(Fold.Enrichment) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Pathway=factor(Pathway, levels=Pathway)) %>%   # This trick update the factor levels
  ggplot( aes(x=Pathway, y=Fold.Enrichment, color=nlog10FDR)) +
  geom_segment( aes(xend=Pathway, yend=0), size = 0.8) +
  geom_point( aes(size=nGenes, color=nlog10FDR)) +
  coord_flip() +
  labs(x="", y= "Fold Enrichment") +
  theme(
    axis.title.x = element_text(size = 7, family = "Arial"), #24
    axis.title.y = element_text(size = 7, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    legend.title =  element_text(size = 7, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.key.size = unit(2,"mm"),
    legend.box.margin=margin(0,0,0,-10),
    plot.margin = unit(c(1,1,1,1), "mm")

  ) + 
  scale_colour_manual(
    name = "-log10(FDR)",
    values = mycolors[c(1,2,4,5,6,8,10)],
    drop = FALSE
    
  ) + 
  guides(size=guide_legend(title="N. of Genes", order = 1), color= guide_legend(order=2)) + scale_size(range = c(0.3,2))



plot.fixed <- set_panel_size(enrichment_plot,
                                  width = unit(40, "mm"),
                                  height = unit(65, "mm"))

#enrichment_plot
#out_file <- paste(base_path, "left-top_gobp_enrichment.tiff", sep="")
#ggsave(out_file, plot = enrichment_plot, dpi = 600, width = 170, height = 100, device='tiff', units = "mm")

#out_file <- paste(base_path, "left-top_gobp_enrichment.png", sep="")
#ggsave(out_file, plot = enrichment_plot, dpi = 600, width = 100, height = 65, units = "mm")

out_file <- paste(base_path, "left-top_gobp_enrichment.png", sep="")
ggsave(out_file, plot = plot.fixed, dpi = 1200, width = 125, height = 75, units = "mm")


