##  Simulation Analysis ####

##  Functions ####
"%!in%" = Negate("%in%")

##  Libraries ####
library(tidyverse)
library(patchwork)

##  Data  ####
peaks = read.table("Input/PeaksAdjusted.txt", stringsAsFactors = F, h = T)
dim(peaks) ## 6584
head(peaks)

peaks_enceph = subset(peaks, peaks$pheno == "encephalopathy")
peaks_enceph = subset(peaks_enceph, peaks_enceph$LOD >= 2.23)
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

enceph.centile = quantile(enceph3$Freq, probs = c(0.95))

##  Plot barchart
enceph.plot = 
  ggplot(data = enceph3) + theme_light() +
  geom_bar(aes(x = Freq),  fill = "#009988", colour="#009988") +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  geom_point(aes(x = 6, y = 50), fill = "#233D4D", shape = 24, size = 3) +
  geom_point(aes(x = enceph.centile, y = 100), colour = "#ee7733", shape = 8, size = 3) +
  labs(x = "Number of peaks predicted", y = "Frequency")

