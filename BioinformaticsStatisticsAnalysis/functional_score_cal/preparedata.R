## read in each file and split to two files. 
##
## if a file has no lib, copy D5 data to be lib 


library(dplyr)
library(tidyr)
library(readr)
library(readxl)

rdir <- 'rpath/pipeline/'
dtdir <- 'rpath/raw_data/'
resdir <- '/rpath/fakedata/'

filenames <- list.files(dtdir) 

fakedata17 <- function(filenamei, dtdir, resdir){
  ## read in data
  test1 <-read_delim(paste0(dtdir, filenamei), 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
  
  print(unique(test1$sample_id))
  snames <- unique(test1$sample_id)
  
  D5name <- snames[grep("_D5$", snames)]
  D8name <- snames[grep("_ola", tolower(snames))]
  
  Libname <- snames[grep("pUCSSMmu_lib", snames)]
  

  test1.drug <- test1 %>% filter(sample_id %in% c(Libname, D5name, D8name[1]))
  test1.drug$sample_id[which(test1.drug$sample_id == D8name[1])] <- "D14"
  
  D8.nodrug.name <- setdiff(snames, c(Libname, D5name, D8name))
  test1.nodrug <- test1 %>% filter(sample_id %in% c(Libname, D5name, D8.nodrug.name))
  test1.nodrug$sample_id[which(test1.nodrug$sample_id == D8.nodrug.name)] <- "D14"
  
  
  write_delim(test1.drug, gsub("tsv", "drug.tsv", paste0(resdir, filenamei)), delim = "\t")
  write_delim(test1.nodrug, gsub("tsv", "nosdrug.tsv", paste0(resdir, filenamei)), delim = "\t")
}

for(fi in filenames[grep("E17", filenames)]){
  fakedata17(fi, dtdir, resdir)
}

## for files except E17.
fakedata <- function(filenamei, dtdir, resdir){
  ## read in data
  test1 <-read_delim(paste0(dtdir, filenamei), 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
  
  print(unique(test1$sample_id))
  snames <- unique(test1$sample_id)
  
  D5name <- snames[grep("_D5$", snames)]
  D8name <- snames[grep("_ola", tolower(snames))]
  
  if(length(grep("pUCSSMmu_lib", snames)) > 0){
    Libname <- snames[grep("pUCSSMmu_lib", snames)]
    test1.drug <- test1 %>% filter(sample_id %in% c(Libname, D5name, D8name))
    test1.drug$sample_id[which(test1.drug$sample_id == D8name)] <- "D14"
    D8.nodrug.name <- setdiff(snames, c(Libname, D5name, D8name))
    test1.nodrug <- test1 %>% filter(sample_id %in% c(Libname, D5name, D8.nodrug.name))
    test1.nodrug$sample_id[which(test1.nodrug$sample_id == D8.nodrug.name)] <- "D14"
  } else{ ## add lib 
    Libname <- "pUCSSMmu_lib"
    temp <- test1 %>% filter(sample_id == D5name)
    temp$EventCount <- temp$EventCount + sample(c(1, 2), size = nrow(temp), replace = TRUE)
    temp$sample_id <- Libname
    test1 <- bind_rows(test1, temp)
    test1.drug <- test1 %>% filter(sample_id %in% c(Libname, D5name, D8name))
    test1.drug$sample_id[which(test1.drug$sample_id == D8name)] <- "D14"
    D8.nodrug.name <- setdiff(snames, c(D5name, D8name))
    test1.nodrug <- test1 %>% filter(sample_id %in% c(Libname, D5name, D8.nodrug.name))
    test1.nodrug$sample_id[which(test1.nodrug$sample_id == D8.nodrug.name)] <- "D14"
  }
  write_delim(test1.drug, gsub("tsv", "drug.tsv", paste0(resdir, filenamei)), delim = "\t")
  write_delim(test1.nodrug, gsub("tsv", "nosdrug.tsv", paste0(resdir, filenamei)), delim = "\t")
}

for(fi in filenames[-grep("E17", filenames)]){
  fakedata(fi, dtdir, resdir)
}
