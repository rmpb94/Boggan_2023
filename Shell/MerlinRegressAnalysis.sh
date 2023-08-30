###############################################################################################################################
##                                               LINKAGE ANALYSIS                                                            ##
###############################################################################################################################
module load merlin
module load R

cd /nobackup/proj/spnmmd/Roisin/3243/Linkage
##  Copy QC files over to the new directory
cd /nobackup/proj/spnmmd/Roisin/3243/Linkage/Input
cp /nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/full-data* .

cd /nobackup/proj/spnmmd/Roisin/3243/Linkage/
##  Copy over the trait file from the R project output files
out_dir="/nobackup/proj/spnmmd/Roisin/3243/Linkage/Output/"
in_dir="/nobackup/proj/spnmmd/Roisin/3243/Linkage/Input"

##  Run pedstats
pedstats -p ${in_dir}FullCohortResidualsForAnalysis.ped,${in_dir}full-data.ped -d ${in_dir}FullCohortResidualsDatFile.dat,${in_dir}full-data.dat -x X --pdf > ${out_dir}full-cohort-pedstats.txt

##  Copy over the phenotype files from the LOCAL OneDrive file into this new directory
cd Jobs
##  Run Merlin
sbatch MerlinAnalysis.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
module load merlin
in_dir="/nobackup/proj/spnmmd/Roisin/3243/Linkage/Input/"
out_dir="/nobackup/proj/spnmmd/Roisin/3243/Linkage/Output/"

merlin-regress -p ${in_dir}FullCohortResidualsForAnalysis.ped,${in_dir}full-data.ped -d ${in_dir}FullCohortResidualsDatFile.dat,${in_dir}full-data.dat -m ${in_dir}full-data.map -x X -t ${in_dir}TraitFileFullCohort.txt --markerNames --tabulate --pdf --prefix ${out_dir}merlin


##  Combine the results files together
R
##  List all files in directory
chrno = c("01", "02", "03", "04", "05", "06", "07", "08", "09", 10:22)
file_list = paste0("/nobackup/proj/spnmmd/Roisin/3243/Linkage/Output/merlin-regress-chr", chrno , ".tbl")
##  Combine the analysis files per chromosome

##  CHR POS LABEL TRAIT H2 LOD PVALUE
##  Make a base dataframe that can be used to merge these single chromosome files together
data = as.data.frame(matrix(ncol = 8, nrow = 0))
names(data) = c("CHR", "LABEL", "PHENOTYPE", "H2", "SD", "INFO", "LOD", "PVALUE")

for(file in file_list){

    temp_dat = read.table(file, stringsAsFactors = F, h = T, row.names = NULL)
    data = rbind(data, temp_dat)
    print(dim(data))
}
names(data) = c("CHR", "LABEL", "TRAIT", "H2", "SD", "INFO", "LOD", "PVALUE")

##  Need to add base pair positions into this file
map = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/bpfull-data.map", stringsAsFactors = F, h = F)
names(map) = c("CHR", "LABEL", "POS", "BP")
map = map[,-which(names(map) %in% c("POS"))]

data = merge(data, map, by = c("CHR", "LABEL"), all.T = T, sort = F)

data$LOD = as.numeric(data$LOD)
data$BP = as.numeric(data$BP)

diabetes = subset(data, data$TRAIT == "diabetes")
encephalopathy = subset(data, data$TRAIT == "encephalopathy")
hearing = subset(data, data$TRAIT == "hearing")
sle = subset(data, data$TRAIT == "sle")

##  Save the individual files
write.table(diabetes, "/nobackup/proj/spnmmd/Roisin/3243/Linkage/Output/FullPhenotypeFiles/diabetes.txt", col.names = F, row.names = F, quote = F)
write.table(encephalopathy, "/nobackup/proj/spnmmd/Roisin/3243/Linkage/Output/FullPhenotypeFiles/encephalopathy.txt", col.names = F, row.names = F, quote = F)
write.table(hearing, "/nobackup/proj/spnmmd/Roisin/3243/Linkage/Output/FullPhenotypeFiles/hearing.txt", col.names = F, row.names = F, quote = F)
write.table(sle, "/nobackup/proj/spnmmd/Roisin/3243/Linkage/Output/FullPhenotypeFiles/sle.txt", col.names = F, row.names = F, quote = F)