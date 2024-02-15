# import package
library(ggplot2)
library(ggExtra)
library(extrafont)

base_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Song2023_additional_works_1/15_additional_work_6/"
stemloop_file <- paste(base_path, "cdur.gencode.v40.pc_transcripts.nopary.cdsonly.refseq.vcloop.csv", sep = "")
vc_file <- paste(base_path, "cdur.gencode.v40.pc_transcripts.vc.csv", sep = "")
stemloop_data <- read.table(file = stemloop_file,  sep = ',', header = TRUE)
vc_data <- read.table(file = vc_file,  sep = ',', header = TRUE)


cdur_stemloop_plot <- ggplot(stemloop_data, aes(x=Motif_under.representation, y=Mutational_susceptibility)) + 
  geom_point(size=0.4, alpha=0.25, fill = "black", stroke = 0.01) +
  geom_hline(yintercept = 0.05, linewidth=0.2) +
  geom_hline(yintercept = 0.95, linewidth=0.2) +
  geom_vline(xintercept = 0.05, linewidth=0.2) +
  geom_vline(xintercept = 0.95, linewidth=0.2)

cdur_stemloop_plot <- cdur_stemloop_plot + 
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

cdur_stemloop_plot <- ggMarginal(
  cdur_stemloop_plot,                  
  type="histogram", 
  size = 10,
  xparams = list(binwidth = 0.01, linewidth = 0.1), 
  yparams = list(binwidth = 0.01, linewidth = 0.1)
)
cdur_stemloop_plot


cdur_vc_plot <- ggplot(vc_data, aes(x=Motif_under.representation, y=Mutational_susceptibility)) + 
  geom_point(size=0.4, alpha=0.25, fill = "black", stroke = 0.01) +
  geom_hline(yintercept = 0.05, linewidth=0.2) +
  geom_hline(yintercept = 0.95, linewidth=0.2) +
  geom_vline(xintercept = 0.05, linewidth=0.2) +
  geom_vline(xintercept = 0.95, linewidth=0.2)

cdur_vc_plot <- cdur_vc_plot + 
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

cdur_vc_plot <- ggMarginal(
  cdur_vc_plot,                  
  type="histogram", 
  size = 10,
  xparams = list(binwidth = 0.01, linewidth = 0.1), 
  yparams = list(binwidth = 0.01, linewidth = 0.1)
)
cdur_vc_plot




out_file <- paste(base_path, "cdur_plot_hsapiens_vc.png", sep="")
ggsave(out_file, plot = cdur_vc_plot, dpi = 600, width = 65, height = 60, units = "mm")
out_file <- paste(base_path, "cdur_plot_hsapines_vc_stemloop.png", sep="")
ggsave(out_file, plot = cdur_stemloop_plot, dpi = 600, width = 65, height = 60, units = "mm")



