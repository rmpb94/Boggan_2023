##  Simulation Analysis- Suggestive peaks ####

##  Functions ####
"%!in%" = Negate("%in%")

##  Libraries ####
library(tidyverse)
library(patchwork)

##  Data  ####
peaks = read.table("Input/PeakBoundariesAdjusted.txt", stringsAsFactors = F, h = T)
dim(peaks) ## 4421
head(peaks)

peaks_enceph = subset(peaks, peaks$pheno == "encephalopathy")
head(peaks_enceph)

real_peaks = read.table("Output/PeaksAdjusted.txt", stringsAsFactors = F, h = T)

##  Peaks per phenotype ####
phenotype = c("diabetes", "encephalopathy", "hearing", "sle")
sims = paste0("LOD_", 1:1000)

##  Boxplot ####
peaks2 = peaks[,-which(names(peaks) %in% c("BP", "LOD"))]
peaks2 = peaks2[!duplicated(peaks2),]

peaks2$pheno = factor(peaks2$pheno, levels = c("diabetes", "encephalopathy", "hearing", "sle"))

##  Table of number of peaks predicted
peak_phenos = c("diabetes", "encephalopathy", "hearing", "sle")
results = as.data.frame(matrix(ncol = 3, nrow = 0))

for(p in peak_phenos){
  
  dat = subset(peaks2, peaks2$pheno == p)
  
  peaks_pred = as.data.frame(table(dat$Sim))
  head(peaks_pred)
  
  names(peaks_pred) = c("sim", "Freq")
  
  peaks_pred$phenotype = p
  
  results = rbind(results, peaks_pred)  
}

head(results)
write.table(results, "Output/PredictedPeaksSimulations.txt", col.names = T, row.names = F, quote = F)

##  Encephalopathy  ####
enceph1 = subset(peaks2, peaks2$pheno == "encephalopathy")
enceph2 = as.data.frame(table(enceph1$Sim))
##  Add number of simulations in which 0 peaks were identified
enceph.nopeaks.length = 1000 - length(enceph2$Var1)
enceph.nopeaks = as.data.frame(matrix(ncol = 2, nrow = enceph.nopeaks.length))
names(enceph.nopeaks) = c("Var1", "Freq")
enceph.nopeaks$Var1 = paste0("LOD", 1000 + 1:enceph.nopeaks.length)
enceph.nopeaks$Freq = 0

enceph3 = rbind(enceph2, enceph.nopeaks)
head(enceph3)
enceph.centile = quantile(enceph3$Freq, probs = c(0.99))

##  Plot barchart
enceph.plot = 
  ggplot(data = enceph3) + theme_light() +
  geom_bar(aes(x = Freq),  fill = "#009988", colour="#009988") +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  labs(x = "Number of peaks predicted", y = "Frequency", title = "Encephalopathy")

##  Diabetes  ####
diabetes1 = subset(peaks2, peaks2$pheno == "diabetes")
diabetes2 = as.data.frame(table(diabetes1$Sim))
##  Add number of simulations in which 0 peaks were identified
diabetes.nopeaks.length = 1000 - length(diabetes2$Var1)
diabetes.nopeaks = as.data.frame(matrix(ncol = 2, nrow = diabetes.nopeaks.length))
names(diabetes.nopeaks) = c("Var1", "Freq")
diabetes.nopeaks$Var1 = paste0("LOD", 1000 + 1:diabetes.nopeaks.length)
diabetes.nopeaks$Freq = 0

diabetes3 = rbind(diabetes2, diabetes.nopeaks)
head(diabetes3)
diabetes.centile = quantile(diabetes3$Freq, probs = c(0.99))
diabetes.actual = 1
##  Plot barchart
diabetes.plot = 
  ggplot(data = diabetes3) + theme_light() +
  geom_bar(aes(x = Freq),  fill = "#009988", colour="#009988") +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  labs(x = "Number of peaks predicted", y = "Frequency", title = "Diabetes")

##  Hearing ####
hearing1 = subset(peaks2, peaks2$pheno == "hearing")
hearing2 = as.data.frame(table(hearing1$Sim))
##  Add number of simulations in which 0 peaks were identified
hearing.nopeaks.length = 1000 - length(hearing2$Var1)
hearing.nopeaks = as.data.frame(matrix(ncol = 2, nrow = hearing.nopeaks.length))
names(hearing.nopeaks) = c("Var1", "Freq")
hearing.nopeaks$Var1 = paste0("LOD", 1000 + 1:hearing.nopeaks.length)
hearing.nopeaks$Freq = 0

hearing3 = rbind(hearing2, hearing.nopeaks)
head(hearing3)
hearing.centile = quantile(hearing3$Freq, probs = c(0.99))

##  Plot barchart
hearing.plot = 
  ggplot(data = hearing3) + theme_light() +
  geom_bar(aes(x = Freq),  fill = "#009988", colour="#009988") +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  labs(x = "Number of peaks predicted", y = "Frequency", title = "Hearing impariment")

##  SLE ####
sle1 = subset(peaks2, peaks2$pheno == "sle")
sle2 = as.data.frame(table(sle1$Sim))
##  Add number of simulations in which 0 peaks were identified
sle.nopeaks.length = 1000 - length(sle2$Var1)
sle.nopeaks = as.data.frame(matrix(ncol = 2, nrow = sle.nopeaks.length))
names(sle.nopeaks) = c("Var1", "Freq")
sle.nopeaks$Var1 = paste0("LOD", 1000 + 1:sle.nopeaks.length)
sle.nopeaks$Freq = 0

sle3 = rbind(sle2, sle.nopeaks)
head(sle3)
sle.centile = quantile(sle3$Freq, probs = c(0.99))

##  Plot barchart
sle.plot = 
  ggplot(data = sle3) + theme_light() +
  geom_bar(aes(x = Freq),  fill = "#009988", colour="#009988") +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  labs(x = "Number of peaks predicted", y = "Frequency", title = "Stroke-like episodes")

pl = enceph.plot + sle.plot + diabetes.plot + hearing.plot + plot_layout(ncol = 2)

tiff("Output/Figures/SimulationPeaks.tiff", units = 'in', width = 6, height = 4, res = 600)
pl
dev.off()
