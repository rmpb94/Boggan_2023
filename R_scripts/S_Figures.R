##  Odds ratios and pseudo-R2 #####

##  Library ####
library(ggplot2)
library(patchwork)

##  Data  ####
or = read.table("Input/OddsRatiosHetAge_261021.txt", stringsAsFactors = F, h = T)
head(or)
mcfr2 = read.table("Input/McFaddenR2_261021.txt", stringsAsFactors = F, h = T)

positions = c("encephalopathy", "sle", "diabetes", "hearing")

or_plot =
  ggplot(or) + theme_light() +
  geom_point(aes(x = phenotype, y = OR_het, colour = "#fe7f2d"), size= 3, position = position_nudge(x = -0.2)) +
  geom_errorbar(aes(x = phenotype, ymax=CI_Lower_het,ymin=CI_Upper_het), position = position_nudge(x = -0.2), colour = "#fe7f2d", width = 0.2, size = 0.8)+
  geom_point(aes(x = phenotype, y = OR_age, colour = "#619b8a"), size = 3, shape = 17) +
  geom_errorbar(aes(x = phenotype, ymax=CI_Lower_age,ymin=CI_Upper_age), colour = "#619b8a", width = 0.2, size = 0.8) +
  geom_point(aes(x = phenotype, y = OR_sex, colour = "#fcca46"), size = 3, shape = 15, position = position_nudge(x = 0.2)) +
  geom_errorbar(aes(x = phenotype, ymax=CI_Lower_sex,ymin=CI_Upper_sex), colour = "#fcca46", width = 0.2, size = 0.8, position = position_nudge(x = 0.2)) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=12)) +
  scale_x_discrete(limits = positions,labels=c("encephalopathy" = "Encephalopathy",
                                                "sle" =  "Stroke-like \nepisodes",
                                                "diabetes" =  "Diabetes",
                                                "hearing" = "Hearing \nimpairment")) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4.5)) +
  scale_colour_identity(guide = "legend", breaks = c( "#fe7f2d", "#619b8a", "#fcca46"), labels = c("m.3243A>G \nlevel", "Age", "Sex")) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, 17, 15)))) +
  labs(colour = "", y = "")

r2_plot = 
  ggplot(mcfr2) + theme_light() +
  geom_point(aes(x = phenotype, y = pseudo_r2), size = 2) +
  scale_x_discrete(limits = positions,labels=c("encephalopathy" = "",
                                               "sle" = "",
                                               "hearing" = "",
                                               "diabetes" = "")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=10)) + 
  labs(x = "", y = expression(paste("Pseudo ", R^2)))

  
or2_plot = r2_plot / or_plot + plot_layout(heights = c(3,7))

tiff("Output/Figures/Covariates.tiff", units = "in", res = 600, height = 5, width = 6)
or2_plot
dev.off()
