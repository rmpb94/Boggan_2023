##  Figures ####
##  Figures for the paper

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
  ggplot(linkage_enceph, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_linedraw() +
  geom_hline(yintercept = e_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = e_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = enceph_axis_set$CHR, breaks = enceph_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  scale_colour_manual(values = rep(c("#009988", "#ee7733"), unique(length(enceph_axis_set$CHR)))) +
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

#enceph_chr_7 = 
  ggplot(subset(linkage_enceph, linkage_enceph$CHR == 7), aes(x = bp_cum, y = LOD)) + theme_linedraw() +
  geom_hline(yintercept = e_sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = e_sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(colour = "#009988", linewidth = 1, alpha = 0.6) +
  geom_point(colour = "#009988", alpha = 0.9, size = 2, shape = 16) +
    scale_x_continuous(breaks=seq(min(enceph7$BP/1000000), max(enceph7$BP/1000000), length.out = 10), labels = scales::number_format(accuracy = 10))
    scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
    scale_colour_manual(values = rep(c("#009988", "#ee7733"), unique(length(enceph_axis_set$CHR)))) +
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
  
