##  Maximum LOD analysis  ####

##  Libraries ####
library(ggplot2)
library(reshape2)

##  Set variables ####
columns = paste0("LOD_", 1:1000)
phenotype = c("diabetes", "encephalopathy", "hearing", "sle")

##  Simulated Data  ####
indir = "Input/Simulated/"

##  Results dataframe ####
results = as.data.frame(matrix(ncol = length(phenotype), nrow = 1000))
names(results) = phenotype

number = as.data.frame(matrix(ncol = 2, nrow = length(phenotype)))
names(number) = c("phenotype", "n")
number$phenotype = phenotype  

for(p in phenotype){
  dat = read.table(paste0(indir, p, ".tbl"), stringsAsFactors = F, h = T)
  max_lods = c()
  for(c in columns){
    lods = dat[[c]]
    ml = max(lods)
    max_lods = c(max_lods, ml)
    
  }
  results[[p]] = max_lods
  no_na = max_lods[!is.na(max_lods)]
  number$n = ifelse(number$phenotype == p, length(no_na), number$n)
}

head(results)

##  Boxplot ####
results$sim = columns
results_long = melt(results, id.vars = c("sim"))
head(results_long)
names(results_long) = c("sim", "phenotype", "LOD")
head(results_long)
results_long$phenotype = factor(results_long$phenotype, levels = c("encephalopathy", "sle", "diabetes", "hearing"))
lod_distribution=
  ggplot(subset(results_long, !is.na(LOD)), aes(x=phenotype, y=LOD)) + theme_light() +
    geom_boxplot(na.rm=T, outlier.shape= NA, alpha = 0.1, fill = "#233D4D") +
    geom_jitter(alpha=0.3, height=0, width=0.3, colour = "#009988", size = 1.5) +
    scale_x_discrete("",labels=c("encephalopathy" = "Encephalopathy",
                                 "sle" = "Stroke-like \nepisodes",
                                 "hearing" = "Hearing \nimpairment",
                                 "diabetes" = "Diabetes")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=11)) +
    geom_hline(yintercept = 3.3, colour = "red") +
    geom_hline(yintercept = 1.86, colour = "grey50", size = 1, alpha = 0.8, linetype = "dashed") +
    ylab("Maximum LOD score")

tiff("Output/Figures/LODDistributions.tiff", units = "in", height = 4, width = 8, res = 600)  
lod_distribution
dev.off()
