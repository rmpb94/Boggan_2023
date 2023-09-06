##############################################################################################
##                                            Simulations                                   ##
##############################################################################################
##  Simulations to calculate the empirical significance of the simulation results

mkdir /nobackup/proj/spnmmd/Roisin/3243/Simulation
cd /nobackup/proj/spnmmd/Roisin/3243/Simulation/

mkdir Input
mkdir Output
mkdir /nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedData
mkdir /nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis

##  Into input file need to put the full-cohort ped, dat, and map files
cp /nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/full-data* /nobackup/proj/spnmmd/Roisin/3243/Simulation/Input/.

##  Put phenotype files and trait files into Input directory from R project
cp /nobackup/proj/spnmmd/Roisin/3243/Linkage/Input/FullCohortResidualsForAnalysis.ped /nobackup/proj/spnmmd/Roisin/3243/Simulation/Input/.
cp /nobackup/proj/spnmmd/Roisin/3243/Linkage/Input/FullCohortResidualsDatFile.dat /nobackup/proj/spnmmd/Roisin/3243/Simulation/Input/.
cp /nobackup/proj/spnmmd/Roisin/3243/Linkage/Input/TraitFileFullCohort.txt /nobackup/proj/spnmmd/Roisin/3243/Simulation/Input/.

mkdir Jobs
cd Jobs

##  Single job  ####
SimulationGeneration.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load merlin
out_dir="/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedData/"
indir="/nobackup/proj/spnmmd/Roisin/3243/Simulation/Input/"

##  Random seed generation
SEED=$(shuf -i 100000-999999 -n 1)

##  Merlin generate simulated files
merlin -p ${indir}full-data.ped -d ${indir}full-data.dat -m ${indir}full-data.map -x X --simulate -r ${SEED} --save --prefix ${out_dir}${SEED}

##  Analysis 
raw_analysis_dir="/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Raw/"
analysis_dir="/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Processed/"
##  Analysis files
pedin="/nobackup/proj/spnmmd/Roisin/3243/Simulation/Input/FullCohortResidualsForAnalysis.ped"
datin="/nobackup/proj/spnmmd/Roisin/3243/Simulation/Input/FullCohortResidualsDatFile.dat"
traitin="/nobackup/proj/spnmmd/Roisin/3243/Linkage/Input/TraitFileFullCohort.txt"
simdir="/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedData/"

merlin-regress -p ${pedin},${simdir}${SEED}-replicate.ped -d ${datin},${simdir}${SEED}-replicate.dat -m ${simdir}${SEED}-replicate.map -x X -t ${traitin} --markerNames --tabulate --prefix ${raw_analysis_dir}${SEED}

##  Merge all analysis files into one file
cat /nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Raw/${SEED}-regress-chr*.tbl >> ${analysis_dir}simulation_${SEED}.tbl

##  Launch job  ####
##  Need to remember to change the permissions to make this executable
chmod +rwx SimulationGeneration.sh
##  Run the job in batches of 100
SimulationLaunch.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
srun /nobackup/proj/spnmmd/Roisin/3243/Simulation/Jobs/SimulationGeneration.sh

##  Test
sbatch --ntasks=2 /nobackup/proj/spnmmd/Roisin/3243/Simulation/Jobs/SimulationLaunch.sh

##  This is how to launch the script on ROCKET 100 times
sbatch --ntasks=100 /nobackup/proj/spnmmd/Roisin/3243/Simulation/Jobs/SimulationLaunch.sh

##############################################################################################
##                                    Simulation Processing                                 ##
##############################################################################################
PhenotypeCombine.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH -p bigmem
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load R
Rscript PhenotypeCombine.R

library(plyr)
file_list = list.files("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Processed/", pattern = ".tbl")
sim_files = lapply(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Processed/", file_list), read.table, header=T, sep = "\t")

##  Map/BP data
##  File needs to have BP poisitions in 
map = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/bpfull-data.map", stringsAsFactors = F, h = F)
names(map) = c("CHR", "LABEL", "POS", "BP")
map = map[,-which(names(map) %in% c("POS"))]

for(i in 1:length(sim_files)){
    ##  Specify which dataframe to use
    data = sim_files[[i]]
    data$PHENOTYPE = sub("Trait: ", "", data$PHENOTYPE)
    ##  Get simulation number for the end filename
    sim_number = file_list[[i]]
    sim_number = substr(sim_number, 1, 16)

    ##  Add the BP column
    ##  Rename the POS col in the data
    data = plyr::rename(data, c("POS" = "LABEL"))
    data = merge(data, map, by = c("CHR", "LABEL"), all.x = T, sort = F)
    ##  Make sure LOD and BP cols are numeric
    data$LOD = as.numeric(data$LOD)
    data$BP = as.numeric(data$BP)

    ##  Split files by phenotype
    ##  Save resulting dataframe
    pheno_list = c("diabetes", "encephalopathy", "hearing", "sle")
    for(p in 1:length(pheno_list)){
        ph = pheno_list[[p]]

        pheno_dat = subset(data, data$PHENOTYPE == ph)
        write.table(pheno_dat, paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/ByPhenotype/", ph, "/", ph, "-", sim_number, ".tbl"), col.names = F, row.names = F, quote = F)
    }
}

##############################################################################################
##                                      Simulation Analysis                                 ##
##############################################################################################
##  Per phenotype, I need to look at the results per 1000 simulations

sbatch AllLODs.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH -p bigmem
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load R
Rscript AllLODs.R

##  This makes one file per phenotype that stores all the LOD scores in one master file
##  I can now use these files for analysis, and not have to worry about individual phenotype files
phenotypes = c("diabetes", "encephalopathy", "hearing", "sle")

for(p in phenotypes){

    indir = paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/ByPhenotype/", p, "/")

    file_list = list.files(indir, pattern = ".tbl")
    data_files = lapply(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/ByPhenotype/", p, "/", file_list), read.table, header=F)

    head(data_files[[1]])
    ##  Rename the LOD columns to make them unique to each simulation 
    ##  This is so they won't overwrite each other in a merge
    for(i in 1:length(data_files)){
        names(data_files[[i]]) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE", "BP")
    }
    ##  Combine all simulated LOD output columns together
    ##  Define a results dataframe
    results = as.data.frame(matrix(ncol = 3, nrow = 8214))
    results = plyr::rename(results, c("V1" = "CHR", "V2" = "LABEL", "V3" = "BP"))
    results$CHR = data_files[[1]]$CHR
    results$LABEL = data_files[[1]]$LABEL
    results$BP = data_files[[1]]$BP

    ##  Results needs a BP_ADJ column which gives the adjusted positions for genome wide positions, not just chromosomal
    results$BP_ADJ = results$BP + 150000000

    ##  All LOD results into one dataframe
    for(i in 1:length(data_files)){
        dat = data_files[[i]]
        print(file_list[i])
        lod_col = dat[,which(names(dat) %in% c("LOD"))]

        results$LOD = lod_col
        results = plyr::rename(results, c("LOD" = paste0("LOD_", i)))
    }
    dim(results)
    results[1:10,1:15]

    ##  Save the results dataframe output
    outdir = paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results/", p, ".tbl")
    write.table(results, outdir, col.names = T, row.names = F, quote = F)
}

##############################################################################################
##                                   Simulation Peak Boundaries                             ##
##############################################################################################
##  Identify linakge peaks in simulations, for each phenotype use empirical suggestive threshold
SimulatedPeakBoundariesDiabetes.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH -p bigmem
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load R
Rscript SimulatedPeakBoundariesDiabetes.R

##  Functions ####
'%notin%' = Negate('%in%')
phenotypes = c("encephalopathy")
columns = paste0("LOD_", 1:1000)
window = 15000000

file_list = list.files("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results", pattern = ".tbl")
data_files = lapply(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results/", file_list), read.table, header=T)

## Diabetes
chr_list = c()
pheno_list = c()
bp_list = c()
up_list = c()
low_list = c()
lod_list = c()
pos = c()
sim_list = c()

dat1 = data_files[[1]]
p = "diabetes"
lod_suggestive = 2.06
  
for(c in columns){
    ##  Look at one LOD column at a time
    dat2 = dat1[,which(names(dat1) %in% c("CHR", "LABEL", "BP", c))]
    ##  Remove any LOD scores lower than suggestive threshold
    dat3 = subset(dat2, dat2[[c]] >= lod_suggestive)
    for(chr in 1:22){
        ##  Subset data to look per chromosome
        dat4 = subset(dat3, dat3$CHR == chr)
        if(length(dat4[[c]]) > 0){
            repeat{
              max_lod = max(dat4[[c]])
              lod_pos = subset(dat4, dat4[[c]] == max_lod)
              chromosome = chr
              lod_pos = lod_pos$BP[1]
              sim = c
              upper_lim = lod_pos + window
              lower_lim = lod_pos - window
              bp_remove = c(lower_lim:upper_lim)
              bp_list = c(bp_list, lod_pos)
              up_list = c(up_list, upper_lim)
              low_list = c(low_list, lower_lim)
              lod_list = c(lod_list, max_lod)
              pos = c(pos, lod_pos)
              sim_list = c(sim_list, c)
              chr_list = c(chr_list, chromosome)
              pheno_list = c(pheno_list, p)
              dat4 = subset(dat4, dat4$BP %notin% bp_remove)
            if(length(dat4$LOD) == 0) break
            }
        }
    }
}

length(bp_list)
length(chr_list)

results = as.data.frame(matrix(ncol = 7, nrow = length(bp_list)))
names(results) = c("pheno", "chr", "LOD", "BP", "BP_Lower", "BP_Upper", "Sim")
results$pheno = pheno_list
results$chr = chr_list
results$LOD = lod_list
results$BP = bp_list
results$BP_Lower = low_list
results$BP_Upper = up_list
results$Sim = sim_list
results
write.table(results, "/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/PeakBoundariesDiabetes.txt", col.names = T, row.names = F, quote = F)

##  Encephalopathy  ####
SimulatedPeakBoundariesEncephaloapthy.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH -p bigmem
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load R
Rscript SimulatedPeakBoundariesEncephalopathy.R

##  Functions ####
'%notin%' = Negate('%in%')
phenotypes = c("encephalopathy")
columns = paste0("LOD_", 1:1000)
window = 15000000

file_list = list.files("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results", pattern = ".tbl")
data_files = lapply(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results/", file_list), read.table, header=T)

## Diabetes
chr_list = c()
pheno_list = c()
bp_list = c()
up_list = c()
low_list = c()
lod_list = c()
pos = c()
sim_list = c()

dat1 = data_files[[2]]
p = "encephalopathy"
lod_suggestive = 2.23
  
for(c in columns){
    ##  Look at one LOD column at a time
    dat2 = dat1[,which(names(dat1) %in% c("CHR", "LABEL", "BP", c))]
    ##  Remove any LOD scores lower than suggestive threshold
    dat3 = subset(dat2, dat2[[c]] >= lod_suggestive)
    for(chr in 1:22){
        ##  Subset data to look per chromosome
        dat4 = subset(dat3, dat3$CHR == chr)
        if(length(dat4[[c]]) > 0){
            repeat{
              max_lod = max(dat4[[c]])
              lod_pos = subset(dat4, dat4[[c]] == max_lod)
              chromosome = chr
              lod_pos = lod_pos$BP[1]
              sim = c
              upper_lim = lod_pos + window
              lower_lim = lod_pos - window
              bp_remove = c(lower_lim:upper_lim)
              bp_list = c(bp_list, lod_pos)
              up_list = c(up_list, upper_lim)
              low_list = c(low_list, lower_lim)
              lod_list = c(lod_list, max_lod)
              pos = c(pos, lod_pos)
              sim_list = c(sim_list, c)
              chr_list = c(chr_list, chromosome)
              pheno_list = c(pheno_list, p)
              dat4 = subset(dat4, dat4$BP %notin% bp_remove)
            if(length(dat4$LOD) == 0) break
            }
        }
    }
}

length(bp_list)
length(chr_list)

results = as.data.frame(matrix(ncol = 7, nrow = length(bp_list)))
names(results) = c("pheno", "chr", "LOD", "BP", "BP_Lower", "BP_Upper", "Sim")
results$pheno = pheno_list
results$chr = chr_list
results$LOD = lod_list
results$BP = bp_list
results$BP_Lower = low_list
results$BP_Upper = up_list
results$Sim = sim_list
results
write.table(results, "/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/PeakBoundariesEncephalopathy.txt", col.names = T, row.names = F, quote = F)

##  Hearing Impairment  ####
SimulatedPeakBoundariesHearing.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH -p bigmem
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load R
Rscript SimulatedPeakBoundariesHearing.R

##  Functions ####
'%notin%' = Negate('%in%')
phenotypes = c("encephalopathy")
columns = paste0("LOD_", 1:1000)
window = 15000000

file_list = list.files("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results", pattern = ".tbl")
data_files = lapply(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results/", file_list), read.table, header=T)

chr_list = c()
pheno_list = c()
bp_list = c()
up_list = c()
low_list = c()
lod_list = c()
pos = c()
sim_list = c()

dat1 = data_files[[3]]
p = "hearing"
lod_suggestive = 2.32
  
for(c in columns){
    ##  Look at one LOD column at a time
    dat2 = dat1[,which(names(dat1) %in% c("CHR", "LABEL", "BP", c))]
    ##  Remove any LOD scores lower than suggestive threshold
    dat3 = subset(dat2, dat2[[c]] >= lod_suggestive)
    for(chr in 1:22){
        ##  Subset data to look per chromosome
        dat4 = subset(dat3, dat3$CHR == chr)
        if(length(dat4[[c]]) > 0){
            repeat{
              max_lod = max(dat4[[c]])
              lod_pos = subset(dat4, dat4[[c]] == max_lod)
              chromosome = chr
              lod_pos = lod_pos$BP[1]
              sim = c
              upper_lim = lod_pos + window
              lower_lim = lod_pos - window
              bp_remove = c(lower_lim:upper_lim)
              bp_list = c(bp_list, lod_pos)
              up_list = c(up_list, upper_lim)
              low_list = c(low_list, lower_lim)
              lod_list = c(lod_list, max_lod)
              pos = c(pos, lod_pos)
              sim_list = c(sim_list, c)
              chr_list = c(chr_list, chromosome)
              pheno_list = c(pheno_list, p)
              dat4 = subset(dat4, dat4$BP %notin% bp_remove)
            if(length(dat4$LOD) == 0) break
            }
        }
    }
}

length(bp_list)
length(chr_list)

results = as.data.frame(matrix(ncol = 7, nrow = length(bp_list)))
names(results) = c("pheno", "chr", "LOD", "BP", "BP_Lower", "BP_Upper", "Sim")
results$pheno = pheno_list
results$chr = chr_list
results$LOD = lod_list
results$BP = bp_list
results$BP_Lower = low_list
results$BP_Upper = up_list
results$Sim = sim_list
results
write.table(results, "/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/PeakBoundariesHearing.txt", col.names = T, row.names = F, quote = F)

##  Stoke-like episodes  ####
SimulatedPeakBoundariesSLE.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH -p bigmem
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load R
Rscript SimulatedPeakBoundariesSLE.R

##  Functions ####
'%notin%' = Negate('%in%')
phenotypes = c("encephalopathy")
columns = paste0("LOD_", 1:1000)
window = 15000000

file_list = list.files("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results", pattern = ".tbl")
data_files = lapply(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/Results/", file_list), read.table, header=T)

chr_list = c()
pheno_list = c()
bp_list = c()
up_list = c()
low_list = c()
lod_list = c()
pos = c()
sim_list = c()

dat1 = data_files[[4]]
p = "sle"
lod_suggestive = 1.94
  
for(c in columns){
    ##  Look at one LOD column at a time
    dat2 = dat1[,which(names(dat1) %in% c("CHR", "LABEL", "BP", c))]
    ##  Remove any LOD scores lower than suggestive threshold
    dat3 = subset(dat2, dat2[[c]] >= lod_suggestive)
    for(chr in 1:22){
        ##  Subset data to look per chromosome
        dat4 = subset(dat3, dat3$CHR == chr)
        if(length(dat4[[c]]) > 0){
            repeat{
              max_lod = max(dat4[[c]])
              lod_pos = subset(dat4, dat4[[c]] == max_lod)
              chromosome = chr
              lod_pos = lod_pos$BP[1]
              sim = c
              upper_lim = lod_pos + window
              lower_lim = lod_pos - window
              bp_remove = c(lower_lim:upper_lim)
              bp_list = c(bp_list, lod_pos)
              up_list = c(up_list, upper_lim)
              low_list = c(low_list, lower_lim)
              lod_list = c(lod_list, max_lod)
              pos = c(pos, lod_pos)
              sim_list = c(sim_list, c)
              chr_list = c(chr_list, chromosome)
              pheno_list = c(pheno_list, p)
              dat4 = subset(dat4, dat4$BP %notin% bp_remove)
            if(length(dat4$LOD) == 0) break
            }
        }
    }
}

length(bp_list)
length(chr_list)

results = as.data.frame(matrix(ncol = 7, nrow = length(bp_list)))
names(results) = c("pheno", "chr", "LOD", "BP", "BP_Lower", "BP_Upper", "Sim")
results$pheno = pheno_list
results$chr = chr_list
results$LOD = lod_list
results$BP = bp_list
results$BP_Lower = low_list
results$BP_Upper = up_list
results$Sim = sim_list
results
write.table(results, "/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/PeakBoundariesSLE.txt", col.names = T, row.names = F, quote = F)


##############################################################################################
##                                   Simulation Peak Adjustment                             ##
##############################################################################################
##  Combine all simulated files together
R
file_list = list.files("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/", pattern = ".txt")
data_files = lapply(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/", file_list), read.table, header=T)
all_data = rbind(data_files[[1]], data_files[[2]], data_files[[3]], data_files[[4]])
write.table(all_data, "/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/PeakBoundariesAll.txt", col.names = T, row.names = F, quote = F)

PeakAdjustment.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load R
Rscript PeakAdjustment.R

##  Data  ####
peaks = read.table("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/PeakBoundariesAll.txt", stringsAsFactors = F, h = T)
head(peaks)
##  Create a Peak_LOD column in current data or I cannot use it later on
peaks$Peak_LOD = peaks$LOD
head(peaks)

##  Replace negative BP values with 0s
peaks$BP_Lower = ifelse(peaks$BP_Lower < 0, 0, peaks$BP_Lower)

simulations = paste0("LOD_", 1:1000)

new_data = as.data.frame(matrix(ncol = 8, nrow = 0))
names(new_data) = c("pheno", "chr", "LOD", "BP", "BP_Lower", "BP_Upper", "Sim", "Peak_LOD")
phenotypes = unique(peaks$pheno)

for(p in phenotypes){
  dat1 = subset(peaks, peaks$pheno == p)
  for(sim in simulations){

      dat2 = subset(dat1, dat1$Sim == sim)
    for(c in 1:22){

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
}

write.table(new_data, "/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/SimulatedAnalysis/PeakBoundariesAdjusted.txt", col.names = T, row.names = F, quote = F)

##############################################################################################
##                                Empirical p-value calculations                            ##
##############################################################################################
sbatch EmpiricalP.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
#SBATCH --mem=10G
#SBATCH -p bigmem
module load R
Rscript EmpP.R

library(reshape2)
phenotype =  c("diabetes", "encephalopathy", "hearing", "sle")
for(p in phenotype){
  dat = read.table(paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/AllLODs/", p, ".tbl"), stringsAsFactors = F, h = T)
  dat = melt(dat, id.vars = c("CHR", "LABEL", "BP", "BP_ADJ"))
  names(dat) = c("CHR", "LABEL", "BP", "BP_ADJ", "Sim", "LOD")

  list1 = unique(dat$LOD)

  results = data.frame()
  print(p)
  for(i in 1:length(list1)){

    print(i)

    dat2 = subset(dat, LOD == list1[i])

    dat2$above = sum(dat$LOD >= list1[i])

    results = rbind(results, dat2)
  }

  write.table(results, paste0("/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/", p, "_empirical_p.txt"), col.names = T, row.names = F, quote = F)
}