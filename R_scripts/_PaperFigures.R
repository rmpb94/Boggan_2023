##  Figures ####

##  Libraries ####
library(ggplot2)
library(tidyverse)
library(patchwork)

##  Data  ####
col_names = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
enceph = read.table("Output/MerlinRegress/encephalopathy.txt", stringsAsFactors = F, h = F)
diab = read.table("Output/MerlinRegress/diabetes.txt", stringsAsFactors = F, h = F)
hear = read.table("Output/MerlinRegress/hearing.txt", stringsAsFactors = F, h = F)
stroke = read.table("Output/MerlinRegress/sle.txt", stringsAsFactors = F, h = F)

##  Add column names
names(enceph) = col_names
names(diab) = col_names
names(hear) = col_names
names(stroke) = col_names

##  Significance thresholds from simulations  ####
e_sig = 3.67
e_sug = 2.23
d_sig = 2.58
d_sug = 2.06
h_sig = 3.88
h_sug = 2.32
s_sig = 2.03
s_sug = 1.94

##  Formatting for graphing ####
##  Encephalopathy
data_enceph = enceph %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)
##  Add new BP position to new dataframe
linkage_enceph = enceph %>%
  inner_join(data_enceph, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)
enceph_axis_set = linkage_enceph %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))
data_diab = diab %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)
##  Add new BP position to new dataframe
linkage_diab = diab %>%
  inner_join(data_diab, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)
diab_axis_set = linkage_diab %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))
data_hear = hear %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)
##  Add new BP position to new dataframe
linkage_hear = hear %>%
  inner_join(data_hear, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)
hear_axis_set = linkage_hear %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))
data_stroke = stroke %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)
##  Add new BP position to new dataframe
linkage_stroke = stroke %>%
  inner_join(data_stroke, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)
stroke_axis_set = linkage_stroke %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

##  Plots ####
enceph_plot = 
  ggplot(linkage_enceph, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = e_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = e_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = enceph_axis_set$CHR, breaks = enceph_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  scale_colour_manual(values = rep(c("#233d4d", "#fe7f2d"), unique(length(enceph_axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "LOD Score", title = "Encephalopathy") +
  theme( 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, colour = "black"))

enceph7 = subset(linkage_enceph, linkage_enceph$CHR == 7)
enceph_chr_7 = 
  ggplot(enceph7, aes(x = BP/1000000, y = LOD)) + theme_bw() +
  geom_hline(yintercept = e_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = e_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(colour = "#233d4d", linewidth = 1, alpha = 0.6) +
  geom_point(colour = "#233d4d", alpha = 0.9, size = 2, shape = 16) +
    scale_x_continuous(breaks=seq(min(enceph7$BP/1000000), max(enceph7$BP/1000000), length.out = 10), labels = scales::number_format(accuracy = 10)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
    scale_colour_manual(values = rep(c("#233d4d", "#fe7f2d"), unique(length(enceph_axis_set$CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "LOD Score", title = "Chromosome 7", y = "Mbp") +
    theme( 
      legend.position = "none",
      axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 12, colour = "black"),
      axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size=12, colour = "black"))


enceph = enceph_plot + enceph_chr_7 + plot_layout(nrow = 2)

tiff("Output/Figures/Encephalopathy.tiff", units = "in", width = 10, height = 6, res = 600)
enceph
dev.off()

##  SLE ####
sle_plot =  
  ggplot(linkage_stroke, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = s_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = s_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = stroke_axis_set$CHR, breaks = stroke_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.1)) +
  scale_colour_manual(values = rep(c("#233d4d", "#fe7f2d"), unique(length(stroke_axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "LOD Score", title = "Stroke-like episodes") +
  theme( 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, colour = "black"))
 
##  SLE Encephalopathy chr5
chr5 = rbind(subset(linkage_enceph, linkage_enceph$CHR == 5), subset(linkage_stroke, linkage_stroke$CHR == 5))
head(chr5) 

chr5_plot = 
  ggplot(chr5, aes(x = BP/1000000, y = LOD, colour = TRAIT)) + theme_bw() +
  geom_point(alpha = 0.9, size = 2, shape = 16) +
  geom_line(linewidth = 1, alpha = 0.6) +
  geom_hline(yintercept = e_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = e_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  scale_x_continuous(breaks=seq(min(chr5$BP/1000000), max(chr5$BP/1000000), length.out = 10), labels = scales::number_format(accuracy = 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  scale_colour_manual(values = rep(c("#233d4d", "#fe7f2d"), unique(length(enceph_axis_set$CHR))), name = "", labels = c("Encephalopathy", "Stroke-like \nepisodes")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "LOD Score", title = "Chromosome 5", y = "Mbp") +
  theme( 
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, colour = "black"))

sle_figure = sle_plot + chr5_plot + plot_layout(ncol = 1)

tiff("Output/Figures/Chr5.tiff", units = "in", width = 10, height = 6, res = 600)
sle_figure
dev.off()

##  Dibetes and Hearing ####
diabetes_plot =  
  ggplot(linkage_diab, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = d_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = d_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = diab_axis_set$CHR, breaks = diab_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  scale_colour_manual(values = rep(c("#233d4d", "#fe7f2d"), unique(length(diab_axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "LOD Score", title = "Diabetes") +
  theme( 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, colour = "black"))

hearing_plot =  
  ggplot(linkage_hear, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = d_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = d_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = hear_axis_set$CHR, breaks = hear_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  scale_colour_manual(values = rep(c("#233d4d", "#fe7f2d"), unique(length(hear_axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "LOD Score", title = "Hearing impairment") +
  theme( 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, colour = "black"))

midd = diabetes_plot + hearing_plot + plot_layout(ncol = 1)

tiff("Output/Figures/MIDD.tiff", units = "in", width = 10, height = 6, res = 600)
midd
dev.off()
