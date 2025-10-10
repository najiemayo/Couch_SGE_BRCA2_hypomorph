#!/usr/local/biotools/r/R-4.2.2/bin/Rscript

## where the r file reside
rpath <- "/rpath/rpgm/pipeline/"
### check here before rerun
## results path for Xexon
res_path <- "/rpath/results/12_2024/pipeline/"
## Gene
Gene <- "BRCA2"
## EventType file chunling provided (some of the variants were annotated differently from the results directly off bioinformatics pipeline)
ETfile <- NULL ##"
## exon with different replicates and pam site information
datalist <- "rpath/exonDataList_121924.txt"
## target file so to control the region
targetfile <- "rpath/target_ind_BRCA2_hypo.csv"
## HDR files from lab
HDRfile <- NULL
## variants to exclude from model building
vexclude <- NULL
## clinvar file for model building
clinvardir <- NULL
clinvarfile <- NULL
## dtdir where the xexon run to fetch individual exon results and individual report run output
dtdir <- "rpath/results/12_2024/pipeline/Ind_full/"
## intermediate path for xexon runs
outdir <- paste0(res_path , "/xexon_inputfiles/")
## fn shortened file name for xexon results
fns <- data.frame(exons = c("15,16,17,18C,18N,19,20,21,22,23,24,25C,25N", "15_drug,16_drug,17_drug,18C_drug,18N_drug,19_drug,20_drug,21_drug,22_drug,23_drug,24_drug,25C_drug,25N_drug"),
                  fn = c("nodrug", "drug"))
vset <- 'allv'
## ind exons run sets
# from run set config.

## cross exons run sets
Score_methods <- c("B") ##c("A", "B", "C", "D") 
Lab_curations <- c("No")
MMs <- c("Normal") 
Outlierremoves <- c("No") 
exfs <- "No" ## c("No","Yes") 
