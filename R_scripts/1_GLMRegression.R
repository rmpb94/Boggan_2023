##  GLM Phenotype data  ####

##  This uses the binary phenotype data for the full cohort
##  Libraries ####
library(arm)
library(plyr)
library(ggplot2)

##  Fucntions ####
"%notin%" = Negate("%in%")

##  Data  ####
ped = read.table("Input/FullCohortBinaryPedigree.txt", stringsAsFactors = F, h = T)
dim(ped)##  796 46

##  Formatting  ####
##  Fix missing sex value 
ped[is.na(ped$sex),]
ped$sex = ifelse(ped$pid == "EX1607228", "M", ped$sex)

##  Missing pheno info, add in
ped[ped$pid == "EXLD24562-01",]
ped$diabetes = ifelse(ped$pid == "EXLD24562-01", 1, ped$diabetes)
ped$diabetes_age = ifelse(ped$pid == "EXLD24562-01", 48, ped$diabetes_age)
ped$hearing = ifelse(ped$pid == "EXLD24562-01", 1, ped$hearing)
ped$hearing_age = ifelse(ped$pid == "EXLD24562-01", 48, ped$hearing_age)
ped$het = ifelse(ped$pid == "EXLD24562-01", 0.807897329429264, ped$het)
ped[ped$pid == "EXLD24562-01",]

## Sex formatting
ped$sex = ifelse(ped$sex == "M", 1, 2)
table(ped$sex)

##  Heteroplasmy info ####
het_data = ped[,which(names(ped) %in% c("famid", "pid", "het"))]
head(het_data)

het_data_miss = subset(het_data, is.na(het_data$het))
het_data_miss = subset(het_data_miss, het_data_miss$pid %notin% c(1:20))
dim(het_data_miss) ## 71 individuals without heteroplasmy values
write.table(het_data_miss, "Output/MissingHeteroplasmy.txt", col.names = T, row.names = F, quote = F)

##  Showing the distribution of heteroplasmy values in the data
ggplot(het_data) + theme_light() +
  geom_histogram(aes(x = het), bins = 10, fill = "#027020") +
  labs(x = "m.3243A>G heteroplasmy level", y = "Count")

##  Create a base pedigree from this initial data file  ####
head(ped)
ped = ped[,-which(names(ped) %in% c("mztwin", "alt_id"))]
pedigree = ped[,1:5]
head(pedigree)

##  Make data GLM-ready ####
##  Remove ANYONE with no heteroplasmy information
ped2 = subset(ped, !is.na(ped$het))
dim(ped2)## 482 44

head(ped2)
ped2 = ped2[,which(names(ped2) %in% c("pid", "famid", "fid", "mid", "sex", "het", "diabetes", "diabetes_age", "hearing", "hearing_age", "encephalopathy", "encephalopathy_age", "sle", "sle_age"))]
head(ped2)

##  Define variables
phenotype = c("diabetes", "hearing", "sle", "encephalopathy")

##  Model 1- Age and Heteroplasmy ####
for(i in phenotype){
  age = paste0(i, "_age")
  data = ped2[,c("pid", "het", age, i, "sex")]
  print(dim(data))
  
  ##  Remove inds with no pheno info
  data[[i]] = as.numeric(data[[i]])
  data = subset(data, !is.na(data[[i]]))
  print(dim(data))
  
  data[[age]] = as.numeric(data[[age]])
  print(dim(data))
  
  ##  Take out anyone with no age of diagnosis
  data = subset(data, !is.na(data[[age]]))
  print(dim(data))
  
  ##  Response variable needs to be coded as a factor
  data[[i]] = as.factor(data[[i]])
  
  ##  Create model
  if(i %notin% c("encephalopathy")){
    data_model = glm(data[[i]] ~ het + data[[age]], family = binomial(link = 'logit'), data = data)
    summary(data_model)
  } else{
    data_model = glm(data[[i]] ~ het + data[[age]] + sex, family = binomial(link = 'logit'), data = data)
    summary(data_model)
  }
  ##  Extract the residuals from the model
  data$dev_res = resid(data_model)
  
  ##  Plot residuals
  tiff(paste0("Output/GLMResidualsGraphs/", i, ".tiff")) 
  binnedplot(fitted(data_model), 
             residuals(data_model, type = "response"), 
             nclass = NULL, 
             xlab = "Expected Values", 
             ylab = "Average residual", 
             main = "Binned residual plot", 
             cex.pts = 0.8, 
             col.pts = 1, 
             col.int = "gray")
  dev.off()
  
  ##  Merge residuals information back into a pedigree file
  tempped = ped2[,c("famid", "pid", "fid", "mid", "sex", "het")]
  dim(tempped)
  data_2 = merge(tempped, data, by = c("pid", "het", "sex"), all.x = T)
  dim(data_2)
  data_ped = data_2[,-which(names(data_2) %in% i)]
  data_ped = data_ped[,c("famid", "pid", "fid", "mid", "sex", "dev_res", age, "het")]
  names(data_ped) = c("famid", "pid", "fid", "mid", "sex", i, age, "het")
  dim(data_ped)
  
  ##  Save the ped file
  write.table(data_ped, paste0("Output/GLMOutput/", i, "_glm_ped.txt"), col.names = T, row.names = F, quote = F)
}

##  Combine into one file
##  Define path for glm files
path = "Output/GLMOutput/"
dat_2 = pedigree[,1:5]

for(i in phenotype){
  dat = read.table(paste0(path, i, "_glm_ped.txt"), stringsAsFactors = F, h = T)
  print(dim(dat))
  dat = dat[,-which(names(dat) %in% c(paste0(i,"_age"), "het"))]
  print(dim(dat))
  print(head(dat))
  
  dat_2 = merge(dat_2, dat, by = c("pid", "famid", "fid", "mid", "sex"), all.x = T)
  print(dim(dat_2))
  print(head(dat_2))
}
head(dat_2)

col_order = c("famid", "pid", "fid", "mid", "sex", "diabetes", "encephalopathy", "hearing", "sle")

dat_2 = dat_2[,col_order]
head(dat_2)
tail(dat_2)

##  Need to make missing values X
dat_2[is.na(dat_2)] = "X"
head(dat_2)

##  Back into pedigree file ####
dat_3 = merge(pedigree, dat_2, by = names(pedigree), all.x = T)
head(dat_3)
dat_3[1:10,1:7]

##  Make missing fid and mid values Xs
dat_3$fid = ifelse(dat_3$fid == 0, "X", dat_3$fid)
dat_3$mid = ifelse(dat_3$mid == 0, "X", dat_3$mid)
dat_3[1:10,1:7]

dat_3[is.na(dat_3)] = "X"
dat_3 = dat_3[,col_order]
head(dat_3)

write.table(dat_3, "Output/FullCohortPhenotypeResiduals.txt", col.names = T, row.names = F, quote = F)

##  Remove un-genotyped individuals ####
##  Also need to remove the phenotypic info of anyone who failed QC

##  Data  ####
ped = read.table("Output/FullCohortPhenotypeResiduals.txt", stringsAsFactors = F, h = T)
dim(ped)##  796 9

keep = read.table("Input/408_inds.txt", stringsAsFactors = F, h = F)
head(keep)
dim(keep)## 408, 5

##  Formatting  ####
tail(keep)
names(keep) = c("famid", "pid", "fid", "mid", "sex")
head(keep)

##  Make a list of the individuals who need their phenotype info kept
inds = keep$pid

##  Remove the phenotype information of anyone not in inds list
phenotypes = c("diabetes", "encephalopathy", "hearing", "sle")

dat = ped

for(p in phenotypes){
  dat[[p]] = ifelse(dat$pid %in% inds, dat[[p]], "X")
}

head(dat)
for(c in phenotypes){
  dat2 = dat[,which(names(dat) %in% c("pid", c))]
  dat2[[c]] = ifelse(dat2[[c]] == "X", 0, 1)
  print(sum(dat2[[c]]))
}

dat = dat[,c("famid", "pid", "fid", "mid", "sex", "diabetes", "encephalopathy", "hearing", "sle")]
dat[1:5,]

##  Save the file ####
write.table(dat, "Output/FullCohortResidualsForAnalysis.ped", col.names = T, row.names = F, quote = F)

##  Dat file  ####
datfile = as.data.frame(matrix(ncol = 2, nrow = 4))
datfile$V1 = "T"
datfile$V2 = c("diabetes", "encephalopathy", "hearing", "sle")
datfile
write.table(datfile, "Output/FullCohortResidualsDatFile.dat", col.names = F, row.names = F, quote = F)

##  Trait file  ####
pheno = read.table("Output/FullCohortPhenotypeResiduals.txt", stringsAsFactors = F, h = T)
names(pheno)
pheno = pheno[,-which(names(pheno) %in% c("famid", "fid", "mid", "sex"))]
phenotypes = names(pheno[,2:5])

##  Make an empty dataframe for the trait file  ####
trait = as.data.frame(matrix(ncol = 5, nrow = 4))
names(trait) = c("trait", "mean", "variance", "heritability", "label")
trait$trait = phenotypes
trait$label = phenotypes

##  Insert binary heritability values from SP 2018  ####
##  If no binary heritability, use NMDAS heritability
h2 = c("0.37", "0.64", "0.52", "0.51")
trait$heritability = h2

##  Loop to calculate the mean and variance ####
##  First part is to extact the phenotype column
##  Then make it numeric, remove all the NA values
##  Calculate the mean and variance
##  Put this into the dataframe
##  There will be warnings, as NAs introcuced by the as.numeric() 
for(p in phenotypes){
  t = as.numeric(pheno[[p]])
  t1 = t[!is.na(t)]
  
  m = mean(t1)
  v = var(t1)
  
  trait$mean = ifelse(trait$trait == p, m, trait$mean)  
  trait$variance = ifelse(trait$trait == p, v, trait$variance)
  
}

head(trait)

##  Save the traitfile  ####
write.table(trait, "Output/TraitFileFullCohort.txt", col.names = F, row.names = F, quote = F)
