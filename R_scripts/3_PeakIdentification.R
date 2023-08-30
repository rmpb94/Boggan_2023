##  Peak identification and adjustment  ####
##  Define peak boundaries for actual data
##  I want to have BP positions for the peaks I've identified in the linkage results
##  And the width of the peak
##  +/- 15Mbp is the window that was decided on

##  Functions ####
'%notin%' = Negate('%in%')

##  Data  ####
indir = c("Input/Results/")
phenotypes = c("diabetes", "encephalopathy", "hearing", "sle")

##  Read in all of the results files
file_list = list.files("Input/Results/", pattern = ".txt")
data_files = lapply(paste0("Input/Results/", file_list), read.table, header=F)

##  Define variables used for storage in loop
chr_list = c()
pheno_list = c()
bp_list = c()
up_list = c()
low_list = c()
lod_list = c()

for(file in 1:length(data_files)){
  pheno_val = phenotypes[[file]]
  
  names(data_files[[file]]) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
  
  data = data_files[[file]]
  data_sub = subset(data, data$LOD >= 1.86)
  
  lods = c()
  pos = c()
  window = 15000000
  
  ##  Do the next step per chromosome
  for(chr in 1:22){
    data_chr = subset(data_sub, data_sub$CHR == chr)
    
    if(length(data_chr$LOD) > 0){
      repeat{
        max_lod = max(data_chr$LOD)
        lod_pos = subset(data_chr, data_chr$LOD == max_lod)
        chromosome = chr
        lod_pos = lod_pos$BP[1]
        
        upper_lim = lod_pos + window
        lower_lim = lod_pos - window
        bp_remove = c(lower_lim:upper_lim)
        
        bp_list = c(bp_list, lod_pos)
        up_list = c(up_list, upper_lim)
        low_list = c(low_list, lower_lim)
        
        lods = c(lods, max_lod)
        pos = c(pos, lod_pos)
        
        chr_list = c(chr_list, chromosome)
        pheno_list = c(pheno_list, pheno_val)
        lod_list = c(lod_list, max_lod)
        
        data_chr = subset(data_chr, data_chr$BP %notin% bp_remove)
        
        if(length(data_chr$LOD) == 0) break
      }
    }
  }
}

length(bp_list)
length(chr_list)

results = as.data.frame(matrix(ncol = 6, nrow = 15))
names(results) = c("pheno", "chr", "LOD", "BP", "BP_Lower", "BP_Upper")
results$pheno = pheno_list
results$chr = chr_list
results$LOD = lod_list
results$BP = bp_list
results$BP_Lower = low_list
results$BP_Upper = up_list
results
##  Replace any negative values of lower limits with 0s
##  This is because I can't have a negative value of BPs
results$BP_Lower = ifelse(results$BP_Lower < 0, 0, results$BP_Lower)
results

##  Save this table
write.table(results, "Output/PeakBoundries.txt", col.names = T, row.names = F, quote = F)

##  Peak Adjustment ####
peaks = read.table("Output/PeakBoundries.txt", stringsAsFactors = F, h = T)
head(peaks)

##  Replace negative BP values with 0s
peaks$BP_Lower = ifelse(peaks$BP_Lower < 0, 0, peaks$BP_Lower)
peaks$Peak_LOD = peaks$LOD
new_data = as.data.frame(matrix(ncol = 7, nrow = 0))
names(new_data) = c("pheno", "chr", "LOD", "BP", "BP_Lower", "BP_Upper", "Peak_LOD")

phenotypes = c("diabetes", "encephalopathy", "hearing", "sle")

for(p in phenotypes){
  dat1 = subset(peaks, peaks$pheno == p)
  for(c in 1:22){
    dat2 = dat1
    dat3 = subset(dat2, dat2$chr == c)
    ##  Order the dataframe smallest to largest
    dat3 = dat3[order(dat3$BP_Lower),]
    
    ##  If there is 1 row
    if(length(dat3$chr) == 1){
      new_data = merge(new_data, dat3, by = names(new_data), all = T)
    }
    
    ##  If there is more than 1 row
    else if(length(dat3$chr) > 1){
      ##  Get a list of the number of rows
      rows = c(1:length(dat3$chr))
      
      ##  Compare the rows
      for(r in rows){
        rowA = r
        rowB = r + 1
        
        ##  If the value of rowB is less than or equal to the length of the df
        ##  I can can compare both rows
        if(rowB <= length(dat3$chr)){
          ##  In this case, there are 2 outcomes
          ##  1. Peaks don't overlap
          ##  2. Peaks overlap
          
          ##  Peaks don't overlap
          if(dat3$BP_Upper[rowA] < dat3$BP_Lower[rowB]){
            ##  Take row A, and put this into new_data
            #row_data = dat3[rowA,]
            #new_data = rbind(new_data, row_data)
          }
          
          ##  Peaks overlap
          else if(dat3$BP_Upper[rowA] > dat3$BP_Lower[rowB]){
            
            ##  The upper boundary of row A is greater than the lower boundary of row B
            ##  This means that the Upper of row A becomes the upper of row B
            ##  And the lower of row B becomes the lower of row A
            ##  ANYWHERE IN DAT3
            dat3$BP_Upper = ifelse(dat3$BP_Upper == dat3$BP_Upper[rowA], dat3$BP_Upper[rowB], dat3$BP_Upper)
            dat3$BP_Lower = ifelse(dat3$BP_Lower == dat3$BP_Lower[rowB], dat3$BP_Lower[rowA], dat3$BP_Lower)
            
            ##  Also want to show the max LOD score across these windows
            lods = subset(dat3, dat3$BP_Lower == dat3$BP_Lower[rowA])
            lods = lods$LOD
            newlod = max(lods)
            dat3$Peak_LOD = ifelse(dat3$BP_Lower == dat3$BP_Lower[rowA], newlod, dat3$LOD)
          }
        }
        if(rowB == length(dat3$chr)){
          ##  If rowB indicates we are at the final row, merge the data with the new_data df
          new_data = merge(new_data, dat3, by = names(dat3), all = T, sort = F)
        }
        ##  If the value of rowB is greater than the length of the df
        else if(rowB > length(dat3$chr)){
          break
        }
      }
    }
  }
}
new_data
new_data2 = new_data[,-which(names(new_data) %in% c("LOD", "BP"))]
new_data2 = new_data2[!duplicated(new_data2),]

##  This gives the adjusted peaks for the actual data
write.table(new_data2, "Output/PeaksAdjusted.txt", col.names = T, row.names = F, quote = F)
