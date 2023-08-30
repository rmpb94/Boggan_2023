##  Full cohort QC

##  The following pipeline is QC for data from the full cohort
##  Initial QC steps have come from SP_QC_PIPELINE
##  SP has already merged the two sets of data into one cohort

##  Modules
module load PLINK
module load R
module load merlin

##  Directory
cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/Input

##  Copy over files I know I'll need from previous QC work
cp /nobackup/proj/spnmmd/Axiom_Aug19/merged/inds_to_remove.txt .
cp /nobackup/proj/spnmmd/Roisin/QC_DATA/inds-to-remove-post-qc.txt .
cp /nobackup/proj/spnmmd/Roisin/QC_DATA/range0_23-26.txt .

cp /nobackup/proj/spnmmd/Roisin/OCT20/SNP/INPUT/CENTIMORGAN_MAP_FILE.txt .
##  Edit range of range file
##  Remove 23
##  Becomes range0_24-26.txt

##  SNP IDs code and files
cp /nobackup/proj/spnmmd/Roisin/OCT19/FILES_FOR_SARAH/SNP_IDS.txt .
cp /nobackup/proj/spnmmd/Roisin/OCT19/FILES_FOR_SARAH/AFFY_IDS.txt .
cp /nobackup/proj/spnmmd/Axiom_Aug19/merged2/batch1_Affy_IDs.txt .
cp /nobackup/proj/spnmmd/Axiom_Aug19/merged2/batch2_Affy_IDs.txt .


##  Data

##  BATCH 1 DATA    
##  Get data that hasn't been through marker QC and rename
cp /nobackup/proj/spnmmd/Roisin/QC_DATA/clean-data-export.bed ./batch1.bed
cp /nobackup/proj/spnmmd/Roisin/QC_DATA/clean-data-export.bim ./batch1.bim
cp /nobackup/proj/spnmmd/Roisin/QC_DATA/clean-data-export.fam ./batch1.fam

##  Recode data to get map file
cd ../Output
plink --bfile ../Input/batch1 --recode --out batch1_recoded

##  BATCH 2 DATA
cd ../Input
##  Following data have individuals trimmed out (those that failed sample QC) but have not undergone marker QC
cp /nobackup/proj/spnmmd/Axiom_Aug19/export2/trimmed_inds.bed ./batch2.bed
cp /nobackup/proj/spnmmd/Axiom_Aug19/export2/trimmed_inds.bim ./batch2.bim
cp /nobackup/proj/spnmmd/Axiom_Aug19/export2/trimmed_inds.fam ./batch2.fam

##  Recode data to get map file
cd ../Output
plink --bfile ../Input/batch2 --recode --out batch2_recoded

##  So that both BATCH 1 and BATCH 2 datasets can be merged
##  Rename SNPs from probe IDs to SNP IDs
cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/RScripts/
Rscript: RENAME_SNPS_B1.R
    ##  Data
    path1 = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Input/")
    path2 = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/")

    affy_ids2 = read.table(paste0(path1, "batch1_Affy_IDs.txt"), stringsAsFactors = F, header = T)
    map = read.table(paste0(path2, "batch1_recoded.map"), stringsAsFactors = F, h = F)
    head(map)
    dim(map)    ##  677351 4 

    map$order = c(1:dim(map)[1])

    ##  Check that all map probe ids are in master list
    sum(map$V2 %in% affy_ids2$probeset_id)
    ##  Yes - 677351
    dim(subset(map, !(V2 %in% affy_ids2$probeset_id))) 
    names(map) = c("chr", "probeset_id", "cM", "bp", "order")
    map2 = merge(map, affy_ids2, all.x=T)
    dim(map)    ##  677351  
    dim(map2)## 677351  

    ##  Check to see if any are missing
    subset(map2, is.na(affy_snp_id))    

    ##  None are missing, so can re-order and export
    map3 = map2[,c("chr", "affy_snp_id", "cM", "bp", "order")]
    map3 = map3[order(map3$order),]
    ##  Remove the order column
    map3 = map3[,c("chr", "affy_snp_id", "cM", "bp")]
    # write out map  - overwrite imported file as I can always generate it again from bam files
    write.table(map3, file=paste0(path2, "batch1_recoded.map"), row.names=F, col.names=F, quote=F)

##  Now do the same for BATCH 2
Rscript: RENAME_SNPS_B2.R
    path1 = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Input/")
    path2 = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/")
    
    affy_ids2 = read.table(paste0(path1, "batch2_Affy_IDs.txt"), stringsAsFactors = F, header = T)
    map = read.table(paste0(path2, "batch2_recoded.map"), stringsAsFactors = F, h = F)
    head(map)
    dim(map)    ##  690777 4 

    map$order = c(1:dim(map)[1])

    ##  Check that all map probe ids are in master list
    sum(map$V2 %in% affy_ids2$probeset_id)
    ##  Yes - 690777
    dim(subset(map, !(V2 %in% affy_ids2$probeset_id))) 
    names(map) = c("chr", "probeset_id", "cM", "bp", "order")
    map2 = merge(map, affy_ids2, all.x=T)
    dim(map)
    dim(map2)

    ##  Check to see if any are missing
    subset(map2, is.na(affy_snp_id))

    ##  None are missing, so can re-order and export
    map3 = map2[,c("chr", "affy_snp_id", "cM", "bp", "order")]
    map3 = map3[order(map3$order),]
    ##  Remove the order column
    map3 = map3[,c("chr", "affy_snp_id", "cM", "bp")]
    # write out map  - overwrite imported file as I can always generate it again from bam files
    write.table(map3, file=paste0(path2,"batch2_recoded.map"), row.names=F, col.names=F, quote=F)

##  Now update both sets of PLINK files with new map files, remake the bed files
cd ../Output
plink --file batch1_recoded --make-bed --out batch1_affy
plink --file batch2_recoded --make-bed --out batch2_affy

plink --bfile batch2_affy --bmerge batch1_affy.bed batch1_affy.bim batch1_affy.fam --make-bed --out merged
##  Next step is to merge files and remove any individuals who are duplicated

    ##  181 people loaded from batch2_affy.fam.
    ##  236 people to be merged from batch1_affy.fam.
    ##  Of these, 236 are new, while 0 are present in the base dataset.
    ##  690777 markers loaded from batch2_affy.bim.
    ##  677351 markers to be merged from batch1_affy.bim.
    ##  Of these, 30487 are new, while 646864 are present in the base dataset.
    ##  Warning: Variants 'Affx-6291744' and 'Affx-52322710' have the same position.
    ##  Warning: Variants 'Affx-8217984' and 'Affx-79397138' have the same position.
    ##  Warning: Variants 'Affx-9279652' and 'Affx-89013337' have the same position.
    ##  215 more same-position warnings: see log file.
    ##  Performing single-pass merge (417 people, 721264 variants).
    ##  Merged fileset written to merged-merge.bed + merged-merge.bim +
    ##  merged-merge.fam .
    ##  721264 variants loaded from .bim file.
    ##  417 people (157 males, 260 females) loaded from .fam.
    ##  417 phenotype values loaded from .fam.
    ##  Using 1 thread (no multithreaded calculations invoked).
    ##  Before main variant filters, 320 founders and 97 nonfounders present.
    ##  Calculating allele frequencies... done.
    ##  Total genotyping rate is 0.943866.
    ##  721264 variants and 417 people pass filters and QC.
    ##  Phenotype data is quantitative.
    ##  --make-bed to merged.bed + merged.bim + merged.fam ... done.

##  Check same-position warnings
more merged.log
wc -l merged.log
##  These have the same position in Axiom, but are often different alleles
##  Only one tends to be flagged as 'best and recommended'
##  There are only 252 warnings, so probably best to just strip these out?
##  Try to list these duplicates. They are likely to be removed when looking at missingness anyway
plink -bfile merged --list-duplicate-vars suppress-first 
more plink.dupvar

##  This only lists those with same alleles
##  Only 2 duplicates
    CHR     POS     ALLELES IDS
    14      23898488        0,G     Affx-10272592 Affx-92042208
    17      41243835        0,G     Affx-13890964 Affx-79376531
# These look like multiallelic SNPs

##  Strip out individuals who have been genotyped twice
##  Previous work identified the sample that has the lowest genotyping rate - remove this one 
##  Quick comparison of different files from SP inds to remove and RB inds remove post qc
##  SP file
more ../Input/inds_to_remove.txt
    ##  FID IID
    ##  803127518F 803127518F
    ##  A3243G0019 e8b2
    ##  A3243G0035 381c
    ##  A3243G0055 f271
    ##  A3243G0062 b7a7
    ##  A3243G0109 c1c0
    ##  A3243G0109 fbe2
    ##  A3243G0119 c63d
    ##  E12 EX1211638

##  RB file
more ../Output/inds-to-remove-post-qc.txt
    ##  A3243G0082      feab
    ##  sing_UCL085     UCL085
    ##  UCL_PED_016     UCL034
    ##  UCL_PED_006     UCL024

##  Just work with SPs file for now
##  These are inds SP has already identified
plink --bfile merged --remove ../Input/inds_to_remove.txt --make-bed --out merged_dups_removed
# --remove: 408 people remaining

##  Marker QC   ####
plink --bfile merged_dups_removed --hardy --out merged_dups_removed
plink --bfile merged_dups_removed --freq --out merged_dups_removed
plink --bfile merged_dups_removed --missing --out merged_dups_removed

Rscript: MISS_MAF_HWE.R
    ##  Explore the SNP QC parameters
    ##  Missingness
    ##  MAF
    ##  HWE
    ##  Libraries
    library(Cairo)
    ##  Data
    path = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/")
    hwe = read.table(file=paste0(path,"merged_dups_removed.hwe"), header=T)
    freq = read.table(file=paste0(path, "merged_dups_removed.frq"), header=T)
    miss = read.table(file=paste0(path,"merged_dups_removed.lmiss"), header=T)
    ##  Merge into one dataset
    dat = merge(hwe, freq)
    dat = merge(dat, miss)
    dim(dat) #   721264   markers
    dat$HWElogP = log(dat$P, base=10)
    dat1 = dat[,c("MAF", "F_MISS", "HWElogP")]
    CairoJPEG(file=paste0(path, "Graphs/HWE_miss_freq.jpeg"))
    plot(dat1)
    dev.off()

    ##  How many missing?
    dim(subset(dat, F_MISS > 0.03)) # 76149
    ##  How many with low freq genos?
    dim(subset(dat, MAF < 0.01)) # 68697    15

##  Remove SNPs with F_MISS>0.03
##  Leave SNPs below HWE and MAF thresholds in for now 
##  Remove once additional family members have been merged in

plink --bfile merged_dups_removed --make-bed --exclude ../Input/range0_23-26.txt --geno 0.03 --out merged_filtered
    ##  Total genotyping rate is 0.944085.
    ##  76149 variants removed due to missing genotype data (--geno).
    ##  645115 variants and 408 people pass filters and QC.

##  Check again for duplicates
plink -bfile merged_filtered --list-duplicate-vars
##  None

##  Now check HWE, MAF, and MISSINGNESS again
plink --bfile merged_filtered --hardy --out merged_filtered
plink --bfile merged_filtered --freq --out merged_filtered
plink --bfile merged_filtered --missing --out merged_filtered

Rscript: MISS_MAF_HWE_CHECK.R
    ##  Libraries
    library(Cairo)
    ##  Data
    path = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/")
    hwe = read.table(file=paste0(path,"merged_filtered.hwe"), header=T)
    freq = read.table(file=paste0(path,"merged_filtered.frq"), header=T)
    miss = read.table(file=paste0(path,"merged_filtered.lmiss"), header=T)

    ##  Merge into one dataset
    dat = merge(hwe, freq)
    dat = merge(dat, miss)
    dim(dat) # 645115 markers

    dat$HWElogP = log(dat$P, base=10)

    dat1 = dat[,c("MAF", "F_MISS", "HWElogP")]

    CairoJPEG(file=paste0(path, "Graphs/HWE_miss_freq_filtered.jpeg"))
    plot(dat1)
    dev.off()

    ##  How many missing?
    dim(subset(dat, F_MISS > 0.03)) # None- filtering has worked!

    dim(subset(dat, P<0.00000001)) # 115 markers
    dim(subset(dat, MAF<0.01))  #   41529

##  Individual QC has already been performed but check missingness and heterozygosity
plink --bfile merged_filtered --exclude ../Input/range0_23-26.txt --range --het --out ind_het

Rscript: HET_MISS.R
    ##  Checking heterozygosity and missingness in full cohort
    ##  Individual QC has already been performed, so this should just confirm this
    ##  Libraries
    library("ggplot2")
    library("ggrepel")

    ##  Data
    path = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/")
    miss = read.table(file=paste0(path,"merged_filtered.imiss"), header=T)
    het = read.table(file=paste0(path,"ind_het.het"), header=T)
    dat = merge(miss, het)

    ##  Caclualte het
    dat$het = (dat$N.NM. - dat$O.HOM.)/dat$N.NM.

    ##  Calculate mean and sd of heterozygosity (het)
    m = mean(dat$het)   ##  0.1949068
    sd = sd(dat$het)    ##  0.002630719

    dat$label = ifelse(dat$het >= (m+(3*sd)) | dat$het <= (m-(3*sd)) | dat$F_MISS > 0.03, as.character(dat$IID), "")
    dat$label2 = ifelse(dat$FID=="COMB001", as.character(dat$IID), "")

    path2 = c("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/Graphs/")

    pdf(paste0(path2, "miss_vs_het.pdf"))
    ggplot(dat, aes(x=F_MISS, y=het)) + theme_bw() +
     geom_point()+
     xlab("Frequency of missing SNP calls") +
     ylab("Mean heterozygosity") +
     geom_vline(xintercept=0.03, linetype="dashed", color="red") +
     geom_hline(yintercept=m, color="green") +
     geom_hline(yintercept=(m + (3*sd)), linetype="dashed", color="red") +
     geom_hline(yintercept=(m - (3*sd)), linetype="dashed", color="red") +
     geom_text_repel(label=dat$label, nudge_y=0.001, nudge_x=0.001, size=2) 
    dev.off()

    pdf(paste0(path2, "miss_vs_het_COMB001.pdf"))
    ggplot(dat, aes(x=F_MISS, y=het)) + theme_bw() +
     geom_point()+
     xlab("Frequency of missing SNP calls") +
     ylab("Mean heterozygosity") +
     geom_vline(xintercept=0.03, linetype="dashed", color="red") +
     geom_hline(yintercept=m, color="green") +
     geom_hline(yintercept=(m + (3*sd)), linetype="dashed", color="red") +
     geom_hline(yintercept=(m - (3*sd)), linetype="dashed", color="red") +
     geom_text_repel(label=dat$label2, nudge_y=0.001, nudge_x=0.001, size=2) 
    dev.off()

    ##  Now also label inds who were problematic in original cohort QC
    ##  Has their position changed?

    dat$label2 = ifelse(dat$het >= (m+(3*sd)) | dat$het <= (m-(3*sd)) | dat$F_MISS > 0.03 |
        as.character(dat$IID)=="UCL024" | as.character(dat$IID)=="UCL034" | as.character(dat$IID)=="UCL022" | as.character(dat$IID)=="UCL058"
        | as.character(dat$IID)=="UCL053" | as.character(dat$IID)=="UCL085", as.character(dat$IID), "")

    pdf(paste0(path2, "miss_vs_het.pdf"))
    ggplot(dat, aes(x=F_MISS, y=het)) + theme_bw() +
     geom_point()+
     xlab("Frequency of missing SNP calls") +
     ylab("Mean heterozygosity") +
     geom_vline(xintercept=0.03, linetype="dashed", color="red") +
     geom_hline(yintercept=m, color="green") +
     geom_hline(yintercept=(m + (3*sd)), linetype="dashed", color="red") +
     geom_hline(yintercept=(m - (3*sd)), linetype="dashed", color="red") +
     geom_text_repel(label=dat$label2, nudge_y=0.001, nudge_x=0.001, size=2) 
    dev.off()

    ##  Leave all inds in for now
    ##  UCL053, 58, 59 - all in same pedigree
    ##  912173544H and EX1204175 - same pedigree and also African according to PCA
    ##  EX1704871 - also African, although a singleton

plink --bfile merged_filtered --recode --out merged_recoded
##  Just want a list of individuals
    R
    genos = fread("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/merged_recoded.ped", stringsAsFactors = F, header = F)
    genos = genos[,1:5]
    write.table(genos,"/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/408_inds.txt", col.names = F, row.names = F, quote = F)


Rscript: FULL_COHORT_PED_FILES.R

##  Need to make a pedigree file for the whole cohort
##  SP has already done this in an R project on the shared drive
##  Just want to run this again to check
##  Transfer files from R project 
##  This is the location on the shared drive
    #!/bin/bash
    #SBATCH -A spnmmd
    #SBATCH --mem=20G # set to about 20G if not on big memory node
    #$ -oe
    #SBATCH --workdir=/nobackup/proj/spnmmd/Roisin/3243/Simulation/Output/
    module load R
    Rscript full_cohort_ped_files.R

    ##  Libraries
    library(ggplot2)
    library(plyr)
    library(data.table)

    ##  Function
    '%notin%' = Negate('%in%')

    ##  Data
    ped = read.table(file="/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Input/FullCohortPedigree_260821.txt", stringsAsFactors = F, header = T)
    head(ped)
    ped$sex = ifelse(ped$sex == "M", 1, 2)
    dim(ped)##  806 7
    ped = ped[,1:6]

    ##  Check the UCL_PED_007 is as expected
    ped[ped$famid == "UCL_PED_007",]
    ##  No, this hasn't been updated in this file, so I need to do that
    ped$famid = ifelse(ped$pid == "UCL018", "sing_UCL018", ped$famid)
    ped$mid = ifelse(ped$pid == "UCL018", 0, ped$mid)
    ped$fid = ifelse(ped$pid == "UCL018", 0, ped$fid)
    ped$fid = ifelse(ped$pid == "UCL030", 0, ped$fid)
    ped$mid = ifelse(ped$pid == "UCL030", 0, ped$mid)
    dim(ped)## 806 6
    ped[ped$famid == "UCL_PED_007",]
    ped = ped[c(-116,-187),]
    dim(ped)##  804 6

    genos = fread("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/merged_recoded.ped", stringsAsFactors = F, header = F)
    genos[400:408, 1:10]
    dim(genos) ##  408 1290236
    ##  Genos is just genetic info for individuals, not whole pedigrees
    ##  Seems to be sequential numbers in affstat for all of the second batch of genotypes
    table(genos$V6)
    ##  Overwrite these to "1"s so that they represent peoples samples that have genotype info available 
    ##  Not sure this is true, but can do the next step anyway
    genos$V6 = 1

    ped_geno = genos[,c(2,5,6)]
    head(ped_geno)

    names(ped_geno) = c("pid", "sex","DNA_avail")
    head(ped_geno)
    dim(ped_geno)## 408 3
    ##  ped_geno just shows the DNA availability for individuals, so anyone who was in genos has a DNA_avail of 1

    names(ped_geno)
    names(ped)
    names(ped) = c("pid", "famid", "fid", "mid", "sex", "aff")
    ##  Aff for all individuals is 1, its the 'phenotype' dummy column
    ##  Merge updated ped with ped_geno and check
    merged_peds = merge(ped_geno, ped, by = "pid", all.y = T)
    head(merged_peds)

    merged_peds$sex_final = ifelse(!is.na(merged_peds$sex.x), merged_peds$sex.x, merged_peds$sex.y)
    merged_peds = merged_peds[,c(-2, -7)]

    merged_peds = plyr::rename(merged_peds, c('sex_final'='sex'))

    merged_peds$DNA_avail = ifelse(is.na(merged_peds$DNA_avail), 0, 1)  
    sum(merged_peds$DNA_avail)
    dim(ped)    ##  804
    dim(merged_peds) ## 804    ##  Both ped and merged peds are the same length

    ##  Put columns and rows into correct order
    merged_peds = merged_peds[order(merged_peds$famid),c("famid", "pid", "fid", "mid", "sex", "aff", "DNA_avail")]

    ##  Write PLINK files for extra individuals not in genos dataset (ie not in merged PLINK files)
    ##  This needs to be done in PLINK as dataset too big otherwise! 
    ##  First isolate those inds who have no genotype information
    extra_inds = subset(merged_peds, DNA_avail==0)
    dim(extra_inds) ##  395 7
    ##  Extract just ped info
    extra_inds = extra_inds[,1:6]
    extra_inds$aff

    ##  Create empty ped and map files for the missing inds
    ##  What should map file look like?
    map = read.table(file="merged_recoded.map")
    ##  Build map and data files with just the first marker
    map_empty = map[1,]
    ped_empty = extra_inds
    head(ped_empty)
    names(ped_empty) = c("V1", "V2", "V3", "V4", "V5", "V6")
    ped_empty$V7 = 0
    ped_empty$V8 = 0

    ##  Write out these new ped and map files and then try to merge them in PLINK!
    write.table(ped_empty, file="extra_inds.ped", row.names=F, col.names=F, quote=F)
    write.table(map_empty, file="extra_inds.map", row.names=F, col.names=F, quote=F)

    ##  Update genos file from PLINK so that this has the correct pedigree infomation in it
    genos[1:10, 1:10]
    ##  Merge in updated pedigree info
    ped_genos_updated = subset(merged_peds, DNA_avail==1)

    ped_genos_updated = ped_genos_updated[,c("famid", "pid", "fid", "mid", "sex", "DNA_avail")]
    ped_genos_updated = plyr::rename(ped_genos_updated, c("pid" = "V2"))
    ##  Exclude columns not wanted from genos
    genos = genos[,c("V1", "V3", "V4", "V5", "V6"):=NULL]

    genos_updated = merge(ped_genos_updated, genos, by = "V2", all.x = T)
    genos_updated[, 1:16]

    ##  Check UCL_PED_007
    genos_updated[genos_updated$famid == "UCL_PED_007", 1:10]
    dim(genos_updated)## 409 1290236
    ##  6650 is in there twice
    genos_updated[genos_updated$V2 == "6650",1:10]
    ##  Remove one of these records
    ind = genos_updated[genos_updated$V2 == "6650", ]
    ind = ind[1,]

    genos_updated = subset(genos_updated, genos_updated$V2 != "6650")
    dim(genos_updated)  ##  4-7

    ## Add 1 6650 record back in
    genos_updated = rbind(genos_updated, ind)

    ##  Change order of columns
    setcolorder(genos_updated, c(2, 1, 3:ncol(genos_updated)))
    ##  Now write this file out
    fwrite(genos_updated, file="genos_updated.ped", row.names=F, col.names=F, quote=F, sep=" ")
    write.table(map, file="genos_updated.map", row.names=F, col.names=F, quote=F)   ## 408
    
##  Something different with founders post this R script
##  Extra individual files are extra_inds.ped and extra_inds.map
##  SP has written these with only one marker present and genotypes set to 0
##  SP has also updated the merged_recoded files in R - called genos_updated
plink --file genos_updated --merge extra_inds.ped extra_inds.map --make-bed --output-missing-phenotype 0 --out merged_extra_inds
    ##  408 people loaded from merged_extra_inds-temporary.fam
    ##  395 people to be merged from extra_inds.ped
    ##  Of these, 395 are new, while 0 are present in the base dataset
    ##  645115 variants loaded from .bim file.
    ##  803 people (390 males, 413 females) loaded from .fam.
    ##  803 phenotype values loaded from .fam.
    ##  Using 1 thread (no multithreaded calculations invoked).
    ##  Before main variant filters, 500 founders and 303 nonfounders present.
    ##  Calculating allele frequencies... done.
    ##  Total genotyping rate is 0.506578.
    ##  645115 variants and 803 people pass filters and QC.
    ##  Among remaining phenotypes, 0 are cases and 803 are controls.

plink --bfile merged_extra_inds --make-bed --hwe 0.00000001 --exclude ../Input/range0_24-26.txt --output-missing-phenotype 0 --out merged_hwe

##  Filter on HWE and then check HWE and allele freq
##  --hwe: 50 variants removed due to Hardy-Weinberg exact test.
##  645065 variants and 805 people pass filters and QC.

plink --bfile merged_hwe --hardy --out merged_hwe
plink --bfile merged_hwe --freq --out merged_hwe
plink --bfile merged_hwe --missing --out merged_hwe

Rscript: MISS_MAF_HWE_CHECK2.R
    ##  Libraries
    library("Cairo")

    ##  Data
    hwe = read.table(file="merged_hwe.hwe", header=T)
    hwe = subset(hwe, TEST=="ALL")
    freq = read.table(file="merged_hwe.frq", header=T)
    miss = read.table(file="merged_hwe.lmiss", header=T)

    dat = merge(hwe, freq)
    dat = merge(dat, miss)
    dim(dat)    ##  645115 14

    dat$HWElogP = log(dat$P, base=10)

    dat1 = dat[,c("MAF", "F_MISS", "HWElogP")]

    CairoJPEG(file="/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/Graphs/HWE_miss_freq_filtered_extra_inds.jpeg")
    plot(dat1)
    dev.off()

    ##  How many are now below thresholds previously used?
    ##  --maf 0.01 --hwe 0.00000001
    dim(subset(dat, MAF<0.01)) ##   48725
    dim(subset(dat, P<0.00000001))  ##  0

##  Check data in PEDSTATS
plink --bfile merged_hwe --recode --output-missing-phenotype 0 --out merged_hwe_recoded

##  Make a data (.dat) file for pedstats/merlin and a map file that is in cM, not Morgans
plink --file merged_hwe_recoded --write-snplist --out snps

Rscript: MAP_FILE.R
    ##  Testing the map file and centimorgan distances
    axiom = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Input/CENTIMORGAN_MAP_FILE.txt", stringsAsFactors = F, h = T, sep = "\t")
    dim(axiom)##    836727  5
    head(axiom)

    ##  Does every snp in the file have a cM position?
    temp1 = subset(axiom, is.na(axiom$cm))
    dim(temp1)##    1169    5
    ##  There are 1169 markers without cm positions

    datfile = read.table("snps.snplist")
    dim(datfile)##  645065  1
    ##  This is a list of markers in the data
    ##  Add marker column 
    datfile[,2] = "M"
    affstat = matrix(0, nrow= 1, ncol = 2)
    affstat = as.data.frame(affstat)
    affstat$V2 = "A"
    affstat$V1 = "DUMMY"
    datfile = rbind(affstat, datfile)
    datfile = datfile[,c(2,1)]
    write.table(datfile, "full-data.dat", col.names = F, row.names = F, quote = F)

    ##  There are less SNPs in the dat file, than in the axiom file
    ##  Subset the axiom file to just include those in the dat file
    datsnps = datfile$V1
    datsnps = datsnps[-1]
    length(datsnps)##   645064
    axiom2 = subset(axiom, axiom$affy_snp_id %in% datsnps)
    dim(axiom2)##   645887
    length(datsnps) - length(axiom2$affy_snp_id)
    ##  There are more markers in the axiom data than in the datsnps
    ##  This might be due to multiple records of the same marker
    temp2 = as.data.frame(table(axiom2$affy_snp_id))
    temp2 = subset(temp2, temp2$Freq > 1)## 808 2
    ##  This temp2 shows that there are 808 markers that have multiple records in the axiom data
    ##  Check the positional info, are they the same?
    dupsnps = temp2$Var1
    temp3 = subset(axiom, axiom$affy_snp_id %in% dupsnps)
    temp3 = temp3[order(temp3$affy_snp_id),]
    ##  Yes, positional information is the same, but they are duplicated as the probeset_id has 2 different ids
    ##  Which set of IDs have been used in each dataset?
    ##  This will affect their ability to be placed

    map = read.table("merged_hwe_recoded.map")
    dim(map)##  645064
    ##  Which IDs are used in the map file?
    ##  This uses the Affx-# ID

    names(map) = c("chr", "snp", "cm", "bp")
    map = map[,-3]
    head(axiom)
    names(axiom) = c("label", "affy_id", "chr", "bp", "cm")

    dim(axiom)  ##  836727  5
    dim(map)    ##  645065  3

    ##  Merge by bp
    dat1 = merge(map, axiom, by = c("chr", "bp"), all.x = T)
    dim(dat1)## 646366  6

    dat1 =  dat1[,which(names(dat1) %in% c("chr", "snp", "cm", "bp"))]
    ##  Need to keep the Affx-# IDs as the labels for the markers

    ##  Who is missing cm distances?
    temp4 = subset(dat1, is.na(dat1$cm))
    dim(temp4)##   13271   4
    head(temp4)
    table(temp4$chr)
    ##  The missing cM values are all from chromosomes 23 and 25
    ##  For now, use the old adjustment of *100 to get a rough cM measurement
    map = read.table("merged_hwe_recoded.map")
    names(map) = c("chr", "snp", "cm", "bp")
    map2 = subset(map, map$snp %in% temp4$snp)
    head(map2)
    map2$cm = map2$cm * 100
    map2 = map2[,which(names(map2) %in% c("snp", "cm"))]
    names(map2) = c("snp", "cm2")

    head(dat1)
    dim(dat1)## 646366  4

    ##  Add new cM distances to the data
    dat2 = dat1

    dat2 = merge(dat1, map2, by = c("snp"), all.x = T)
    dim(dat2)

    ##  Create final cm column
    dat2$finalcm = ifelse(!is.na(dat2$cm), dat2$cm, dat2$cm2)
    temp = subset(dat2, is.na(dat2$finalcm))
    ##  Now there are no markers missing cm distances
    dat2 = dat2[,which(names(dat2) %in% c("snp", "chr", "bp", "finalcm"))]

    names(dat2) = c("snp", "chr", "bp", "cm")

    dat1 = dat1[,c("chr", "snp", "cm", "bp")]
    dat1 = dat1[,-4]
    write.table(dat1, "full-data.map", col.names = F, row.names = F, quote = F)
    # R script to produce file for leaving only 10 SNPs in ped file
    datfile = datfile[1:11,]
    write.table(datfile, "remove_all_SNPs_but_10.txt", col.names = F, row.names = F, quote = F)

##  Run pedstats to check parental information is correct
##  First create ped file with no marker data
plink --file merged_hwe_recoded --extract remove_all_SNPs_but_10.txt --output-missing-phenotype 0 --recode --out merged_empty

##  Run this through pedstats
pedstats -p merged_empty.ped -d remove_all_SNPs_but_10.txt -x0 --verbose -- byFamily --familyPDF 

##  All individuals are connected

##  Need to identify mendelian inconsistancies
##  Do this with a script on Rocket
sbatch MERLIN_ERROR.sh
    #!/bin/bash
    #SBATCH -A spnmmd
    #SBATCH -p bigmem # may not be necessary
    #SBATCH --mem=100G # set to about 20G if not on big memory node - job actually used about 25G
    #$ -oe
    cd ../Output
    module load merlin
    merlin -p merged_hwe_recoded.ped -m full-data.map -d full-data.dat -x0 --error 
    echo Finishing job

##  Find out how much memory was actually used:
sacct -o reqmem,maxrss,averss,elapsed -j 
##  Ran for ~ 3.15 hrs

##  Now use R to look at the error file that this rocket script outputs
Rscript: MENDELIAN_INCONSISTANCY.R

    dat = read.table(file="merlin.err", header=T)
    table = table(dat$MARKER)

    pdf(file="/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/Graphs/hist.pdf")
    hist(table)
    dev.off()

    ##  Most markers have <6 errors, look at top end
    table = data.frame(table)
    pdf(file="/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/Graphs/hist2.pdf")
    hist(subset(table, Freq>=5)$Freq)
    dev.off()

    ##  Looks like Freq>6 might be a good cutoff 
    dim(subset(table, Freq>6)) ##   189 markers in total
    ##  Write out a list of these and strip them out of the genotype files
    snps_to_remove = subset(table, Freq > 6)
    snp_list = snps_to_remove[,1]

    write.table(snp_list, "snps_to_remove.txt", col.names = F, row.names = F, quote = F)

##  Now run through pedstats for mendelian errors
sbatch PEDSTATS.sh
    #!/bin/bash
    #SBATCH -A spnmmd
    #SBATCH --mem=20G # set to about 20G if not on big memory node
    #$ -oe
    #SBATCH --mail-user=b8035576@ncl.ac.uk
    module load merlin
    cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/
    pedstats -p ./merged_hwe_recoded.ped -d ./full-data.dat -x0

##  Mendelian errors written to slurm-[JOBID].out
##  13378177 
##  Transfer to desktop to open - so that can be reformatted (in excel!)
##  On OneDrive m3243/LinkageAnalysis/
##  Open and reformat
##  Saved as mendelian_errors.txt in same dir as above
##  Look at in R on rocket
##  Are there any markers (or pedigrees / individuals) that appear more often?
Rscript: MERLIN_ERROR_V2.R
    mend = read.table("mendelian_errors_030921.txt", header = F, stringsAsFactors = F, sep = "\t")
    dim(mend)  ##   7108    3
    names(mend) = c("Marker", "famid", "pid")

    ##  Is there any bias?
    markers = data.frame(table(mend$Marker))    ##  5511 2
    peds = data.frame(table(mend$famid))    ##  53 2

    ##  Graph
    pdf("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/Graphs/hist_markers_peds.pdf")
    hist(markers$Freq, breaks = 25)
    hist(peds$Freq, breaks = 25)
    dev.off()

    ##  Look at ind peds
    subset(peds, Freq>250)  ##  COMB001
    subset(mend, famid == "COMB001")

    table(data.frame(table(subset(mend, famid == "COMB001")$Marker))$Freq)
    ##    1   2   3
    ##  262  52  20

    ##  Which markers have 3 errors? Are these the same  as the ones identified above

    subset(peds, Freq > 150)
    ##               Var1 Freq
    ##  9   A3243G0016  249
    ##  12  A3243G0025  184
    ##  13  A3243G0030  232
    ##  14  A3243G0035  171
    ##  19  A3243G0053  200
    ##  22  A3243G0091  224
    ##  23  A3243G0097  210
    ##  24  A3243G0103  249
    ##  25  A3243G0105  155
    ##  27  A3243G0109  163
    ##  29  A3243G0119  193
    ##  30  A3243G0126  247
    ##  31  A3243G0127  157
    ##  32     COMB001  426
    ##  36         E08  223
    ##  44         G06  235
    ##  45         G07  159
    ##  47 UCL_PED_001  188
    ##  48 UCL_PED_003  161
    ##  51 UCL_PED_011  151
    ##  
    ##  SP had a quick look at these peds
    ##  It suggested that peds with 3 sibs and a mother are more likely to have more errors
    ##  Less informative families have fewer errors
    ##  But, no pedigree seems to have an excessive number
    subset(data.frame(table(subset(mend, famid == "COMB001")$pid)), Freq>0)
    #     Var1 Freq
    #1   cebd   95
    #2   df4b  102
    #3 UCL044  229

##########

##  Now need IBD analysis

##  First prune SNPs for LD
plink  --file merged_hwe_recoded --exclude ../Input/range0_24-26.txt --range --indep 50 5 2 --out prunedsnplist
##  342797 of 643846 variants removed.

##  Run extract command to create bed files with only pruned SNPs present
plink -file merged_hwe_recoded --extract prunedsnplist.prune.in --make-bed --out prunedsnps

plink --bfile  prunedsnps --Z-genome --out relationships

Rscript: IBD.R
    ##  Libraries
    library(ggplot2)
    library(ggrepel)
    library(Cairo)

    ##  Data
    ibd = read.table(gzfile("relationships.genome.gz"), header=T)
    dim(ibd)
    # 322003 pairs of people

    ibd$label = ifelse(ibd$FID1 == "COMB001" & ibd$FID2 =="COMB001" & ibd$IID1!="1"& ibd$IID2!="1",
    paste(as.character(ibd$IID1),"-",as.character(ibd$IID2), sep=""), "")

    CairoTIFF("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/Graphs/ibd.tiff")
    ggplot(ibd, aes(x=Z0, y=Z1)) + theme_bw() +
    geom_point(aes(colour = RT)) +
    scale_colour_manual(values = c("FS" = "#619B8A", "HS" = "#FCCA46","PO" = "#A1C181", "OT" = "#FE7F2D", "UN" = "#233D4D"), labels = c("FS" = "Full siblings", "HS" = "Half siblings", "PO" = "Parent-offspring", "OT" = "Avuncular", "UN" = "Unrelated")) +
    labs(fill = "Relationship")
    dev.off()

##  These look fine on IBD graph. What about missingness and heterozygosity?
##  Added this in above and they look totally fine

##  Error removal
##  Combine info about mendelian and unlikely errors to produce a list of likely problematic markers
Rscript: ERROR_REMOVAL.R
    ##  Data
    mend = read.table(file="mendelian_errors_030921.txt", header=F, sep = "\t")
    unlikely = read.table(file="merlin.err", header=T)

    names(mend) = c("MARKER", "FAMILY", "PERSON")

    ##  Combine lists of markers and get a frequency distribution of these
    markers = c(as.character(mend$MARKER), as.character(unlikely$MARKER)) 
    markers = data.frame(table(markers))
    head(markers)
    table(markers$Freq)

    ##  Most markers have very few errors - maybe strip out those with 6+ errors?
    ##  Look at cluster plots first
    subset(markers, Freq>5)
    ##    1     2     3     4     5     6     7     8     9    10    11    12    13 
    ##   7275 19496  8350  1531   844   301   103    55    18    22    10    10    12 
    ##     14    15    16    17    18    19    20    21    22    23    26    27 
    ##      7     8     2     6     5     5     2     3     1     1     1     1 

    ##  SP visually inspected the cluster plots for these and they look fine - keep them in
    length(unique(mend$MARKER)) # 5511
    dat = data.frame(table(mend$MARKER))

    table(dat$Freq)
    ##     1    2    3    4    5    6    7    8    9   10   11 
    ##  4461  755  171   61   35   12    6    5    1    2    2 

    ##  Look at markers with >=5
    dim(subset(dat, Freq>=5)) # 63 in total
    ##  Write a list of these out and exclude them

    ##  Do similar with unlikely genos
    length(unique(unlikely$MARKER)) # 33664 markers in total
    dat_u = data.frame(table(unlikely$MARKER))
    table(dat_u$Freq)

    ##  How much overlap? How about looking at only those with >=5
    markers_m = subset(dat, Freq>=5)$Var1 #  63
    markers_u = subset(dat_u, Freq>=5)$Var1 # 1178
    sum(markers_m %in% markers_u) # 16
    markers_m[markers_m %in% markers_u]
    ##   [1] Affx-10497960 Affx-12087970 Affx-12302706 Affx-12554348 Affx-13367540
    ##   [6] Affx-16067476 Affx-16162568 Affx-20538928 Affx-21555970 Affx-24860594
    ##  [11] Affx-25961937 Affx-2736686  Affx-29475144 Affx-30498389 Affx-5356591 
    ##  [16] Affx-5440769 
    ##  Probably best to exlude all of these - some cluster polts do look a bit strange
    ##  Write out list of markers to exclude - all those with either >5 mendelian or >5 unlikely errors
    markers_all = unique(c(as.character(markers_m), as.character(markers_u)))# 346 markers
    write.table(as.character(markers_all), file="exclude_mendelian_errors.txt", row.names=F, col.names=F, quote=F)

##  Run PLINK to remove SNPs
plink --bfile merged_hwe --exclude exclude_mendelian_errors.txt --make-bed --out merged_cleaner
##  643848 variants and 805 people pass filters and QC

plink --bfile merged_cleaner --hardy --out merged_cleaner
plink --bfile merged_cleaner --freq --out merged_cleaner
plink --bfile merged_cleaner --missing --out merged_cleaner

##  Re-run missing, HWE and freq
Rscript: MISS_MAF_HWE_CHECK3.R
    ##  Libraries
    library("Cairo")
    
    ##  Data
    hwe = read.table(file="merged_cleaner.hwe", header=T)
    hwe = subset(hwe, TEST=="ALL")
    freq = read.table(file="merged_cleaner.frq", header=T)
    miss = read.table(file="merged_cleaner.lmiss", header=T)
    
    dat = merge(hwe, freq)
    dat = merge(dat, miss)
    dim(dat)    ##  644728
    
    dat$HWElogP = log(dat$P, base=10)
    
    dat1 = dat[,c("MAF", "F_MISS", "HWElogP")]
    
    CairoJPEG(file="/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/Graphs/HWE_miss_freq_cleaner.jpeg")
    plot(dat1)
    dev.off()

    ##  How many are now below thresholds previously used?
    ##  --maf 0.01 --hwe 0.00000001
    dim(subset(dat, MAF<0.01)) ##   45776
    dim(subset(dat, P<0.00000001))

##  Linkage Analysis Data Prep  ####
##  Need ped and map files
plink --bfile merged_cleaner --recode --out merged_cleaner

Rscript: UPDATE_CM.R
    ##  Libraries   ####
    ##  Data    ####
    new_map = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Input/CENTIMORGAN_MAP_FILE.txt", stringsAsFactors = F, h = T, sep = "\t")
    map = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/merged_cleaner.map", stringsAsFactors = F, h = F)

    ##  Formatting  ####
    ##  File has been reformatted in excel
    names(new_map) = c("label", "snp", "chr", "bp", "cm")
    dim(new_map)## 836727 5     

    ##  Names for map
    names(map) = c("chr", "snp_old", "cm_old", "bp")

    ##  Merge together  ####
    dat1 = merge(map, new_map, by = c("chr", "bp"))

    ##  Need Affy-# formats ####
    dat1 = dat1[,which(names(dat1) %in% c("snp", "cm"))]

    ##  Check for duplicates
    length(unique(dat1$snp))
    ##  There are 631916 variants, but only 630725 are unique
    ##  Which ones are these?
    temp = data.frame(table(dat1$snp))
    temp2 = subset(temp, temp$Freq > 1)
    dim(temp2)  ##  1031
    table(temp2$Freq)
    #  2   3   4   6   9
    #990  10  11  17   3

    ##  Make a list of these
    dups = temp2$Var1
    ##  Extract from dat1
    dat2 = subset(dat1, dat1$snp %in% dups)
    ##  It seems like the dups all have the same cM positions
    ##  So, just need to extract one copy of everything
    dat3 = unique(dat1)
    dim(dat3)   ##  630725

    write.table(dat3, "NEW_MAP_030921.txt", col.names = F, row.names = F, quote = F)


plink --file merged_cleaner --update-cm NEW_MAP_030921.txt 2 1 --recode --out cm_updated

##  Now prepare data for Linkage Analysis   ####

##  MAF filtering
plink --file cm_updated --maf 0.35 --recode --out maf_filtered
##  563049 variants removed due to minor allele threshold(s)
##  80799 variants and 803 people pass filters and QC.

##  LD Pruning
plink  --file maf_filtered --exclude ../Input/range0_24-26.txt --range --indep 50 5 2 --out prunedsnplist
##  Pruning complete.  38975 of 80545 variants removed.
##  Extract the pruned SNPs
plink -file maf_filtered --extract prunedsnplist.prune.in --make-bed --out prunedsnps

plink --bfile prunedsnps --recode --out prunedsnps
##  41570 variants and 803 people pass filters and QC

##  MapThin
/nobackup/proj/spnmmd/mapthin/mapthin-v1.11-linux-x86_64/mapthin -t 2.4 prunedsnps.map mapthinmap.map
    ##  Parameters:
    ##  Input file: prunedsnps.map
    ##  Output file: mapthinmap.map
    ##  SNPs per cM: 2.4
    ##  
    ##  Statistics: 
    ##  Total number of SNPs in original file: 41570
    ##  Number of SNPs in thinned file: 8214 (19.7594%)
    ##  
    ##  Mean genetic distance between SNPs: 0.428768 cM
    ##  St. dev. of genetic distance between SNPs: 0.192281 cM
    ##  Range of genetic distances between SNPs: (0.0028, 4.1123)
 
##  8214 SNPs in new map file
##  Now need to create the SNP file with the thinned SNPs in
##  Also set mendelian errors to missing
plink --file prunedsnps --extract mapthinmap.map --set-me-missing --make-bed --out merged_thinned
## -set-me-missing: 167 errors addressed.
##  Get ped and map files
plink --bfile merged_thinned --recode --out merged_thinned

mkdir /nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/
cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/
cp /nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/merged_thinned* .

cp merged_thinned.map bpfull-data.map
cp merged_thinned.ped full-data.ped


##  Need to make a full-data dat file for the markers
R
dat1 = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/full-data.map", stringsAsFactors = F, h = F)
dim(dat1)## 8214    4
dat2 = dat1[,1:2]
dat2$V1 = "M"
write.table(dat2, "/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/full-data.dat", col.names = F, row.names = F, quote = F)

##  Remove affection status column from the pedigree file
R
ped1 = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/full-data.ped", stringsAsFactors = F, h = F)
ped1[1:10,1:8]
ped2 = ped1[,-6]
ped2[1:10,1:8]
##  In pedigree also change 0s for missing to "X"a
ped2[ped2==0] = "X"
write.table(ped2, "/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/full-data.ped", col.names = F, row.names = F, quote = F)

##  Map file remove BP column
R
map = read.table("/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/bpfull-data.map", stringsAsFactors = F, h = F)
map2 = map[,-4]
write.table(map2, "/nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/full-data.map", col.names = F, row.names = F, quote = F)

##  PCA analysis    ####
cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/PCA/Input
module load BCFtools
module load PLINK

sbatch bcf.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
#SBATCH --mem=10G
#SBATCH -p bigmem
cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/PCA/Input
module load BCFtools
for chr in {1..22}; do
    bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
      ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
      bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both \
          > ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;

    bcftools index ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;
done

sbatch plinkformat.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
#SBATCH --mem=10G
#SBATCH -p bigmem
cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/PCA/Input
module load PLINK
for chr in {1..22}; do
    plink --noweb \
      --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
done

sbatch prune.sh
#!/bin/bash
#SBATCH -A spnmmd
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nrb177@ncl.ac.uk
#SBATCH --mem=10G
#SBATCH -p bigmem
cd /nobackup/proj/spnmmd/Roisin/3243/SNPQC/PCA/Input
mkdir Pruned

module load PLINK
for chr in {1..22}; do
    plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
      --maf 0.10 --indep 50 5 1.5 \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;

    plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
      --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.prune.in \
      --make-bed \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
done

##  Combine into one file
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;
sed -i 's/.bim//g' ForMerge.list ;
plink --merge-list ForMerge.list --out Merge ;

##  Find common variants between 1000G and my data
cp /nobackup/proj/spnmmd/Roisin/3243/SNPQC/FinalData/merged_thinned* .
cp /nobackup/proj/spnmmd/Roisin/3243/SNPQC/Output/prunedsnps* .

##  Convert 1000G files to map and ped
plink --bfile Merge --recode --out Merge

##  Convert my data to map and ped
plink --bfile prunedsnps --recode --out mydata

##  Extract SNP list from actual data
plink --file mydata --maf 0.05 --write-snplist --out mydata
##  Extract SNP list from 1000Gs
plink --bfile Merge --maf 0.05 --write-snplist --out 1000G

##  In R make file to update snps
##  Transfer to rocket

##  Now update snps in my data files
plink --file mydata --update-name updatesnpspca.txt 2 1 --make-bed --out mydataupdated

##  Extract matching SNPs from 1000G data
plink --bfile mydataupdated --maf 0.05 --write-snplist --out snpsforpca
plink --bfile Merge --extract snpsforpca.snplist --make-bed --out pcafiles
plink --bfile pcafiles --maf 0.05 --write-snplist --out newsnplist
plink --bfile mydataupdated --extract newsnplist.snplist --make-bed --out mydatapca

##  Combine data
plink --bfile pcafiles --bmerge mydatapca.bed mydatapca.bim mydatapca.fam --make-bed --out combined1000g

##  Then run PCA on combined
plink --bfile combined1000g --pca --mind

##  Transfer to R