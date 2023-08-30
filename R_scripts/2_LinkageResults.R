##  Plot linkage results  ####

##  Libraries ####
library(ggplot2)
library(tidyverse)

##  Encephalopathy  ####
enceph = read.table("Output/MerlinRegress/encephalopathy.txt", stringsAsFactors = F, h = F)
names(enceph) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
head(enceph)
enceph_dat = enceph
sig = 3.67
sug = 2.23
data_enceph = enceph %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)

##  Add new BP position to new dataframe
linkage_enceph = enceph %>%
  inner_join(data_enceph, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

head(linkage_enceph)

enceph_axis_set = linkage_enceph %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

#linkage_plot_enceph = 
  ggplot(linkage_enceph, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
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
  
##  Encephalopathy  ####
enceph = read.table("Output/MerlinRegress/encephalopathy.txt", stringsAsFactors = F, h = F)
names(enceph) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
head(enceph)
enceph_dat = enceph
sig = 3.67
sug = 2.23
data_enceph = enceph %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)

##  Add new BP position to new dataframe
linkage_enceph = enceph %>%
  inner_join(data_enceph, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

head(linkage_enceph)

enceph_axis_set = linkage_enceph %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

linkage_plot_enceph = 
  ggplot(linkage_enceph, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
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
 

##  SLE  ####
sle = read.table("Output/MerlinRegress/sle.txt", stringsAsFactors = F, h = F)
names(sle) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
head(sle)
sle_dat = sle
sig = 2.03
sug = 1.94
data_sle = sle %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)

##  Add new BP position to new dataframe
linkage_sle = sle %>%
  inner_join(data_sle, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

head(linkage_sle)

sle_axis_set = linkage_sle %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

linkage_plot_sle = 
ggplot(linkage_sle, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = sle_axis_set$CHR, breaks = sle_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.3)) +
  scale_colour_manual(values = rep(c("#009988", "#ee7733"), unique(length(sle_axis_set$CHR)))) +
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

##  Diabetes  ####
diabetes = read.table("Output/MerlinRegress/diabetes.txt", stringsAsFactors = F, h = F)
names(diabetes) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
head(diabetes)
diabetes_dat = diabetes
sug = 2.06
sig = 2.58
data_diabetes = diabetes %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)

##  Add new BP position to new dataframe
linkage_diabetes = diabetes %>%
  inner_join(data_diabetes, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

head(linkage_diabetes)

diabetes_axis_set = linkage_diabetes %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

linkage_plot_diabetes = 
  ggplot(linkage_diabetes, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = diabetes_axis_set$CHR, breaks = diabetes_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  scale_colour_manual(values = rep(c("#009988", "#ee7733"), unique(length(diabetes_axis_set$CHR)))) +
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


##  Hearing  ####
hearing = read.table("Output/MerlinRegress/hearing.txt", stringsAsFactors = F, h = F)
names(hearing) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
head(hearing)
hearing_dat = hearing
sig = 3.88
sug = 2.32
data_hearing = hearing %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric((max_bp))), default = 0)) %>%
  select(CHR, bp_add)

##  Add new BP position to new dataframe
linkage_hearing = hearing %>%
  inner_join(data_hearing, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

head(linkage_hearing)

hearing_axis_set = linkage_hearing %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

linkage_plot_hearing = 
  ggplot(linkage_hearing, aes(x = bp_cum, y = LOD, colour = as_factor(CHR))) + theme_bw() +
  geom_hline(yintercept = sig, colour = "#cc3311", linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = sug, colour = "grey48", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_line(aes(colour = as.factor(CHR)), linewidth = 1, alpha = 0.6) +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_x_continuous(label = hearing_axis_set$CHR, breaks = hearing_axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  scale_colour_manual(values = rep(c("#009988", "#ee7733"), unique(length(hearing_axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "LOD Score", title = "Hearing loss") +
  theme( 
    legend.position = "none",
    axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, colour = "black"))
