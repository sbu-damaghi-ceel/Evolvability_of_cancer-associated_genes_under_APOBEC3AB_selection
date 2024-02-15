#install.packages('colormap')
library(ggplot2)
library(colormap)
library(magick)
library(ggraph)
library(igraph)
library(reshape2)
library(plyr)
library(EvoFreq)
library(colorspace)
library(gridExtra)
library(extrafont)

#### all files params ###
ending_timepoint <- 600
nsim <- 100

exp_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Spacial_model/A3_mutation_rate_1e-6/"

exp_names <- c("Single_A3_Off", "Single_A3_On")

times <- c()
genotype_heterogeneity <- c()
fitness_heterogeneity <- c()
groupname <- c()
genotype_nums <- c()

for (exp_name in exp_names){
  setwd(paste(exp_path, exp_name, "/", sep = ""))

  info_file <- paste(exp_path, exp_name, "/", "_infos.csv", sep = "")
  infos_df <- read.csv(info_file, check.names = F, header = T)
  
  basal_name <- names(infos_df)[1]
  basename <- strsplit(basal_name, split="seed")[[1]][1]
  
  all_seeds <- unlist(infos_df[,1])
  
  for (i in 1:nsim) {
    seed <- all_seeds[i]
    
    everything_file <- paste(exp_path, exp_name, "/", basename, "seed", as.character(seed), "sim0.csv", sep = "")
    everything_df <- read.csv(everything_file, check.names = F, header = F)
    
    time <- unlist(everything_df[1,])
    genotypes <- unlist(everything_df[31,])
    fitnesses <- unlist(everything_df[30,])
    nums <- unlist(everything_df[32,])
    exps <- rep(exp_name, length(time))
    
    times <- c(times, time)
    genotype_heterogeneity <- c(genotype_heterogeneity, genotypes)
    fitness_heterogeneity <- c(fitness_heterogeneity, fitnesses)
    groupname <- c(groupname, exps)
    genotype_nums <- c(genotype_nums, nums)
    
  }
  

}

everything_plot_df <- data.frame(
 time = times,
 genotype_heterogeneity = genotype_heterogeneity,
 number_of_genotypes = genotype_nums,
 fitness_heterogeneity = fitness_heterogeneity,
 experiment = groupname
 
)

# pvals heterogeneity
pvals <- c()
for (i in seq(from = 50, to = 600, by=50)) {
  print(i)
  tt1 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Single_A3_Off")
  tt2 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Single_A3_On")
  t_result <- t.test(tt1$genotype_heterogeneity, tt2$genotype_heterogeneity)
  
  pvals <- c(pvals, -log10(t_result$p.value))
  
}

pval_df <- data.frame(
  time = seq(from = 50, to = 600, by=50),
  p = pvals
)

pval_plot <- ggplot(pval_df, aes(x=time, y=p)) + 
  geom_point() +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "red") +
  scale_y_continuous(breaks = c(2,5,10,15,20), labels =c(2,5,10,15,20))  +
  ylab("-log10(P)") + xlab("Time") +
  theme_light()+
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial")
  )
pval_plot
ggsave(paste(exp_path, "Figures/pval_onoff_heterogeneity.png", sep = ""), pval_plot, dpi = 1200, width = 85, height = 60, units = "mm")

# pvals num genotypes
pvals <- c()
for (i in seq(from = 50, to = 600, by=50)) {
  print(i)
  tt1 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Single_A3_Off")
  tt2 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Single_A3_On")
  t_result <- t.test(tt1$number_of_genotypes , tt2$number_of_genotypes)
  
  pvals <- c(pvals, -log10(t_result$p.value))
  
}

pval_df <- data.frame(
  time = seq(from = 50, to = 600, by=50),
  p = pvals
)

pval_plot <- ggplot(pval_df, aes(x=time, y=p)) + 
  geom_point() +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "red") +
  scale_y_continuous(breaks = c(2,5,10,15,20), labels =c(2,5,10,15,20))  +
  ylab("-log10(P)") + xlab("Time") +
  theme_light()+
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial")
  )
pval_plot
ggsave(paste(exp_path, "Figures/pval_onoff_genoyptes.png", sep = ""), pval_plot, dpi = 1200, width = 85, height = 60, units = "mm")

exp_names <- c("Uniform_A3_On", "Skewed_A3_On")

# pvals heterogeneity
pvals <- c()
for (i in seq(from = 50, to = 600, by=50)) {
  print(i)
  tt1 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Uniform_A3_On")
  tt2 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Skewed_A3_On")
  t_result <- t.test(tt1$genotype_heterogeneity, tt2$genotype_heterogeneity)
  
  pvals <- c(pvals, -log10(t_result$p.value))
  
}

pval_df <- data.frame(
  time = seq(from = 50, to = 600, by=50),
  p = pvals
)

pval_plot <- ggplot(pval_df, aes(x=time, y=p)) + 
  geom_point() +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "red") +
  scale_y_continuous(breaks = c(2,5,10,15,20), labels =c(2,5,10,15,20))  +
  ylab("-log10(P)") + xlab("Time") +
  theme_light()+
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial")
  )
pval_plot
ggsave(paste(exp_path, "Figures/pval_uniskew_heterogeneity.png", sep = ""), pval_plot, dpi = 1200, width = 85, height = 60, units = "mm")

# pvals num genotypes
pvals <- c()
for (i in seq(from = 50, to = 600, by=50)) {
  print(i)
  tt1 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Uniform_A3_On")
  tt2 <- everything_plot_df %>% arrange(time) %>% group_by(time, experiment) %>% filter(time == i) %>% filter(experiment == "Skewed_A3_On")
  t_result <- t.test(tt1$number_of_genotypes , tt2$number_of_genotypes)
  
  pvals <- c(pvals, -log10(t_result$p.value))
  
}

pval_df <- data.frame(
  time = seq(from = 50, to = 600, by=50),
  p = pvals
)

pval_plot <- ggplot(pval_df, aes(x=time, y=p)) + 
  geom_point() +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "red") +
  scale_y_continuous(breaks = c(2,5,10,15,20), labels =c(2,5,10,15,20))  +
  ylab("-log10(P)") + xlab("Time") +
  theme_light()+
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial")
  )
pval_plot
ggsave(paste(exp_path, "Figures/pval_uniskew_genoyptes.png", sep = ""), pval_plot, dpi = 1200, width = 85, height = 60, units = "mm")


everything_plot_df <- everything_plot_df %>% 
  arrange(time) %>% group_by(time, experiment) %>% 
  mutate(
    num_avg = mean(number_of_genotypes), 
    num_std = sd(number_of_genotypes),
    gen_avg = mean(genotype_heterogeneity), 
    gen_std = sd(genotype_heterogeneity), 
    fit_avg = mean(fitness_heterogeneity), 
    fit_std = sd(fitness_heterogeneity)) %>% 
  select(c(experiment, time, num_avg, num_std, gen_avg, gen_std, fit_avg, fit_std)) %>%
  distinct(experiment, time, num_avg, num_std, gen_avg, gen_std, fit_avg, fit_std)

num_plot <- ggplot(everything_plot_df, aes(x=time, y=num_avg, ymin=num_avg-num_std, ymax=num_avg+num_std )) + 
  geom_line(aes(color = experiment)) +
  geom_ribbon(aes(x=time, y=num_avg, ymin=num_avg-num_std, ymax=num_avg+num_std, fill=experiment), alpha=0.25, linewidth = 0) + 
  ylab("N. of genotypes") + xlab("Time") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    legend.title =  element_blank(),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom",
    legend.key.size = unit(2,"mm"),
    legend.box.margin=margin(-10,0,0,0),
    panel.grid.minor = element_blank(),
  )+ scale_x_continuous(limits=c(0,ending_timepoint), expand=c(0,0))+ scale_y_continuous(limits=c(0,900), expand=c(0,0)) +
  scale_fill_discrete(labels=c("Without APOBEC3A/B", "With APOBEC3A/B")) +
  scale_color_discrete(labels=c("Without APOBEC3A/B", "With APOBEC3A/B"))
num_plot

ggsave(paste(exp_path, "Figures/number_of_genoyptes_onoff.png", sep = ""), num_plot, dpi = 1200, width = 85, height = 60, units = "mm")


div_plot <- ggplot(everything_plot_df, aes(x=time, y=gen_avg, ymin=gen_avg-gen_std, ymax=gen_avg+gen_std )) + 
  geom_line(aes(color = experiment)) +
  geom_ribbon(aes(x=time, y=gen_avg, ymin=gen_avg-gen_std, ymax=gen_avg+gen_std, fill=experiment), alpha=0.25, linewidth = 0) + 
  ylab("Heterogeneity") + xlab("Time") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    legend.title =  element_blank(),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom",
    legend.key.size = unit(2,"mm"),
    legend.box.margin=margin(-10,0,0,0),
    panel.grid.minor = element_blank(),
  )+ scale_x_continuous(limits=c(0,ending_timepoint), expand=c(0,0))+ scale_y_continuous(limits=c(0,15), expand=c(0,0)) +
  scale_fill_discrete(labels=c("Without APOBEC3A/B", "With APOBEC3A/B")) +
  scale_color_discrete(labels=c("Without APOBEC3A/B", "With APOBEC3A/B"))
div_plot

ggsave(paste(exp_path, "Figures/genetic_heterogeneity_onoff.png", sep = ""), div_plot, dpi = 1200, width = 85, height = 60, units = "mm")


#####################################################
exp_names <- c("Uniform_A3_On", "Skewed_A3_On")

times <- c()
genotype_heterogeneity <- c()
fitness_heterogeneity <- c()
groupname <- c()
genotype_nums <- c()

for (exp_name in exp_names){
  setwd(paste(exp_path, exp_name, "/", sep = ""))
  
  info_file <- paste(exp_path, exp_name, "/", "_infos.csv", sep = "")
  infos_df <- read.csv(info_file, check.names = F, header = T)
  
  basal_name <- names(infos_df)[1]
  basename <- strsplit(basal_name, split="seed")[[1]][1]
  
  all_seeds <- unlist(infos_df[,1])
  
  for (i in 1:nsim) {
    seed <- all_seeds[i]
    
    everything_file <- paste(exp_path, exp_name, "/", basename, "seed", as.character(seed), "sim0.csv", sep = "")
    everything_df <- read.csv(everything_file, check.names = F, header = F)
    
    time <- unlist(everything_df[1,])
    genotypes <- unlist(everything_df[31,])
    fitnesses <- unlist(everything_df[30,])
    nums <- unlist(everything_df[32,])
    exps <- rep(exp_name, length(time))
    
    times <- c(times, time)
    genotype_heterogeneity <- c(genotype_heterogeneity, genotypes)
    fitness_heterogeneity <- c(fitness_heterogeneity, fitnesses)
    groupname <- c(groupname, exps)
    genotype_nums <- c(genotype_nums, nums)
    
  }
  
  
}

everything_plot_df <- data.frame(
  time = times,
  genotype_heterogeneity = genotype_heterogeneity,
  number_of_genotypes = genotype_nums,
  fitness_heterogeneity = fitness_heterogeneity,
  experiment = groupname
  
)

everything_plot_df <- everything_plot_df %>% 
  arrange(time) %>% group_by(time, experiment) %>% 
  mutate(
    num_avg = mean(number_of_genotypes), 
    num_std = sd(number_of_genotypes),
    gen_avg = mean(genotype_heterogeneity), 
    gen_std = sd(genotype_heterogeneity), 
    fit_avg = mean(fitness_heterogeneity), 
    fit_std = sd(fitness_heterogeneity)) %>% 
  select(c(experiment, time, num_avg, num_std, gen_avg, gen_std, fit_avg, fit_std)) %>%
  distinct(experiment, time, num_avg, num_std, gen_avg, gen_std, fit_avg, fit_std)

num_plot <- ggplot(everything_plot_df, aes(x=time, y=num_avg, ymin=num_avg-num_std, ymax=num_avg+num_std )) + 
  geom_line(aes(color = experiment)) +
  geom_ribbon(aes(x=time, y=num_avg, ymin=num_avg-num_std, ymax=num_avg+num_std, fill=experiment), alpha=0.25, linewidth = 0) + 
  ylab("N. of genotypes") + xlab("Time") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    legend.title =  element_blank(),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom",
    legend.key.size = unit(2,"mm"),
    legend.box.margin=margin(-10,0,0,0),
    panel.grid.minor = element_blank(),
  )+ scale_x_continuous(limits=c(0,ending_timepoint), expand=c(0,0))+ scale_y_continuous(limits=c(0,1100), expand=c(0,0)) +
  scale_fill_discrete(labels=c("Skewed distribution", "Uniform distribution")) +
  scale_color_discrete(labels=c("Skewed distribution", "Uniform distribution"))
num_plot

ggsave(paste(exp_path, "Figures/number_of_genoyptes_uniskew.png", sep = ""), num_plot, dpi = 1200, width = 85, height = 60, units = "mm")


div_plot <- ggplot(everything_plot_df, aes(x=time, y=gen_avg, ymin=gen_avg-gen_std, ymax=gen_avg+gen_std )) + 
  geom_line(aes(color = experiment)) +
  geom_ribbon(aes(x=time, y=gen_avg, ymin=gen_avg-gen_std, ymax=gen_avg+gen_std, fill=experiment), alpha=0.25, linewidth = 0) + 
  ylab("Heterogeneity") + xlab("Time") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 8, family = "Arial"), #24
    axis.title.y = element_text(size = 8, family = "Arial"), #24
    axis.text.x = element_text(size = 6, family = "Arial"), # 18
    axis.text.y = element_text(size = 6, family = "Arial"), # 18
    legend.title =  element_blank(),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom",
    legend.key.size = unit(2,"mm"),
    legend.box.margin=margin(-10,0,0,0),
    panel.grid.minor = element_blank(),s
  )+ scale_x_continuous(limits=c(0,ending_timepoint), expand=c(0,0))+ scale_y_continuous(limits=c(0,17), expand=c(0,0)) +
  scale_fill_discrete(labels=c("Skewed distribution", "Uniform distribution")) +
  scale_color_discrete(labels=c("Skewed distribution", "Uniform distribution"))
div_plot

ggsave(paste(exp_path, "Figures/genetic_heterogeneity_uniskew.png", sep = ""), div_plot, dpi = 1200, width = 85, height = 60, units = "mm")

##############################################################################
ending_timepoint <- 600

exp_path <- "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Spacial_model/A3_mutation_rate_1e-6/"
exp_names <- c("Single_A3_Off", "Single_A3_On")
exp_names <- c("Uniform_A3_On", "Skewed_A3_On")

for (exp_name in exp_names){
  setwd(paste(exp_path, exp_name, "/", sep = ""))
  
  info_file <- "./_infos.csv"
  infos_df <- read.csv(info_file, check.names = F, header = T)
  
  basal_name <- names(infos_df)[1]
  basename <- strsplit(basal_name, split="seed")[[1]][1]
  
  all_seeds <- unlist(infos_df[,1])

  for (i in c(77)){
    seed <- all_seeds[i]
    
    clone_history_file <- paste("./", basename, "seed", as.character(seed), "sim0_clones.csv", sep = "")
    clone_parents_file <- paste("./", basename, "seed", as.character(seed), "sim0_parents.csv", sep = "")
    
    clone_df <- read.csv(clone_history_file, row.names = 1, check.names = F, header = T)
    parent_df <- read.csv(clone_parents_file, check.names = F, header = F)
    parent_list <- as.numeric(parent_df$V1)
    clone_list <- as.numeric(row.names(clone_df))
    time_pts <- colnames(clone_df)
    
    pos_df <- get_evofreq(
      size_df = clone_df, 
      clones = clone_list, 
      parents = parent_list, 
      threshold=0.04, 
      scale_by_sizes_at_time = FALSE, 
      clone_cmap = "blackbody", 
      time_pts = as.numeric(time_pts),
      interpolation_steps = 0
    )
    
    if (exp_name == "Uniform_A3_On") { 
      titlename = "Uniform distribution"
    } else {
      titlename = "Skewed distribution"
      }
    
    fp <- plot_evofreq(
      freq_frame = pos_df, 
      bw=0,    # width of lines surrounding polygons 
      bc="black", # color of lines surrounding polygons
      end_time = ending_timepoint
    )
    fp <- fp + scale_x_continuous(limits=c(0,ending_timepoint), expand=c(0,0))+scale_y_continuous(expand=c(0,0)) # remove labels
    fp <- fp + guides(fill="none", color="none", alpha="none") + xlab("Time") + ylab("Population Size") + ggtitle(titlename)+
      theme(
        plot.title = element_text(hjust = 0.5,size = 8, family = "Arial"),
        axis.title.x = element_text(size = 8, family = "Arial"), #24
        axis.title.y = element_text(size = 8, family = "Arial"), #24
        axis.text.x = element_text(size = 6, family = "Arial"), # 18
        axis.text.y = element_text(size = 6, family = "Arial"), # 18
      )
    
    ggsave(paste("./", basename, "seed", as.character(seed), "sim0_fish_plot_Freq_Th004b.png", sep = ""), fp, dpi = 1200, width = 85, height = 40, units = "mm")
    
  }
  
}



