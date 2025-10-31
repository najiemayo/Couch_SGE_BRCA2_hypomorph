library(dplyr)
library(tidyr)
library(data.table)

loadOneExperimentSNV = function(inputFile, SNVOnly = T) {
  # load one experiment and keep only SNV, add SNV ID, and keep counts 
  # for D5 >= 10 and lib >= 10, and add sampleID2, librarySize and frequency
  
  dataRaw = as.data.frame(fread(inputFile))
  print(dim(dataRaw))
  
  if (SNVOnly) {
    # keep the SNV
    index = (dataRaw$REF %in% c("A", "C", "T", "G")) & 
      (dataRaw$ALT %in% c("A", "C", "T", "G")) 
    data = dataRaw[index, , drop = F]
  } else {
    data = dataRaw
  }
  SNVID = paste0(data$X.CHROM, "-", data$POS, "-", data$REF, "-", data$ALT)
  data$SNVID = SNVID
  
  # keep SNVs at least 10 in lib and at least 10 in D5
  sampleID2 = data$sample_id
  sampleID2 = sub(".*_", "", sampleID2)
  sampleID2[sampleID2 == "drug"] = "0nM"
  data$sampleID2 = sampleID2
  
  SNVs = intersect(data$SNVID[data$sampleID2 == "D5" & data$EventCount >= 10],
                   data$SNVID[data$sampleID2 == "lib" & data$EventCount >= 10])
  dataFiltered = data[data$SNVID %in% SNVs, , drop = F]
  
  # calculate the library size for each sample
  librarySize = matrix(NA, dim(dataFiltered)[1], 1)
  frequency = matrix(NA, dim(dataFiltered)[1], 1)
  uniqueSampleID = unique(data$sample_id)
  for (i in 1:length(uniqueSampleID)) {
    index = (dataFiltered$sample_id == uniqueSampleID[i])
    size = sum(dataFiltered$EventCount[index])
    librarySize[index] = size
    frequency[index] = dataFiltered$EventCount[index]/size
  }
  dataFiltered$librarySize = c(librarySize)
  dataFiltered$frequency = c(frequency)
  
  return(dataFiltered)
}

loadOneExperimentSNVMinilib = function(inputFile, SNVOnly = F) {
  # load one experiment and keep only SNV, add SNV ID, and keep counts 
  # for D5 >= 10 and lib >= 10, and add sampleID2, librarySize and frequency
  
  dataRaw = as.data.frame(fread(inputFile))
  print(dim(dataRaw))
  
  if (SNVOnly) {
    # keep the SNV
    index = (dataRaw$REF %in% c("A", "C", "T", "G")) & 
      (dataRaw$ALT %in% c("A", "C", "T", "G")) 
    data = dataRaw[index, , drop = F]
  } else {
    data = dataRaw
  }
  SNVID = paste0(data$X.CHROM, "-", data$POS, "-", data$REF, "-", data$ALT)
  data$SNVID = SNVID
  
  # keep SNVs at least 10 in lib and at least 10 in D5
  # After processing, we will have lib D5 D8 Olaparib_D8
  # for E17NEW_hypomini_BRCA2.KH584_R2.tsv with _a and _b with drug
  # choose Olaparib_D8_a
  data = data[!grepl("Olaparib_D8_b", data$sample_id), , drop = F]
  
  sampleID2 = data$sample_id
  print(unique(sampleID2))
  sampleID2 = sub(".*_D5", "D5", sampleID2)
  sampleID2[grepl("Olaparib|ola", sampleID2) & 
              !grepl("Olaparib_D8_b", sampleID2)] = "Olaparib_D8"
  sampleID2[grepl("D8", sampleID2) & !grepl("Olaparib|ola", sampleID2)] = 
      "D8"
  sampleID2[grepl("lib", sampleID2)] = "lib"
  print(unique(sampleID2))
  
  data$sampleID2 = sampleID2
  
  if ("lib" %in% data$sampleID2) {
    SNVs = intersect(data$SNVID[data$sampleID2 == "D5" & data$EventCount >= 10],
                   data$SNVID[data$sampleID2 == "lib" & data$EventCount >= 10])
  } else { # no lib
    SNVs = data$SNVID[data$sampleID2 == "D5" & data$EventCount >= 10]
  }
  dataFiltered = data[data$SNVID %in% SNVs, , drop = F]
  
  # calculate the library size for each sample
  librarySize = matrix(NA, dim(dataFiltered)[1], 1)
  frequency = matrix(NA, dim(dataFiltered)[1], 1)
  uniqueSampleID = unique(data$sample_id)
  for (i in 1:length(uniqueSampleID)) {
    index = (dataFiltered$sample_id == uniqueSampleID[i])
    size = sum(dataFiltered$EventCount[index])
    librarySize[index] = size
    frequency[index] = dataFiltered$EventCount[index]/size
  }
  dataFiltered$librarySize = c(librarySize)
  dataFiltered$frequency = c(frequency)
  
  return(dataFiltered)
}

plotData = function(data, outputFile) {
  # plot data
  
  # convert data to a wide format
  dataToWide = data[, c("X.CHROM", "POS", "REF", "ALT", "SNVID", 
                        "EventType", "sampleID2", "frequency")]
  dataWide = dataToWide %>% pivot_wider(names_from = sampleID2, 
                                        values_from = frequency)
  
  # make sure there is no obvious position effect by plotting
  # plot the D10 vs D5
  columnDrug = c("0nM", "1nM", "5nM", "25nM", "100nM")
  eventColor = character(dim(dataWide)[1])
  eventColor[] = "gray"
  eventColor[dataToWide$EventType == "Missense"] = "blue"
  eventColor[dataToWide$EventType == "Synonymous"] = "green"
  eventColor[dataToWide$EventType == "StopGain"] = "red"
  
  pdf(outputFile, 
      height = 3*length(columnDrug), width = 10)
  par(mfrow = c(length(columnDrug) + 1, 1))
  # plot D5 vs lib
  D5lib = log2(as.matrix(dataWide$D5) / dataWide$lib)
  plot(dataWide$POS, D5lib, col = eventColor, 
       xlab = "position", ylab = "D5 vs lib", main = "no drug")
  # plot D10 vs D5
  for (i in 1:length(columnDrug)) {
    D10D5 = log2(as.matrix(dataWide[, columnDrug[i]]) / dataWide$D5)
    plot(dataWide$POS, D10D5, col = eventColor, 
         xlab = "position", ylab = "D10 vs D5", main = columnDrug[i])
  }
  dev.off()
  
}

loadData = function(inputFileList, outputFolder, 
                    dataName = "MAVE", SNVOnly = T) {
  # load multiple experiments
  
  # input
  # inputFileList: a file with header, 
  #                the first column is the input file, 
  #                the second column is exon ID,
  #                the third column is the replicate ID
  # outputFolder: 
  # dataName: "MAVE" or "MiniLib"
  
  fileInfo = read.table(inputFileList, header = T, as.is = T, sep = "\t")
  dataAll = list()
  for (i in 1:dim(fileInfo)[1]) {
    inputFile = fileInfo[i, 1]
    exonID = fileInfo[i, 2]
    replicateID = fileInfo[i, 3]
    print(exonID)
    if (dataName == "MAVE") {
      data = loadOneExperimentSNV(inputFile, SNVOnly)
    } else if (dataName == "MiniLib") {
      data = loadOneExperimentSNVMinilib(inputFile, SNVOnly)
    }
    
    # plot data
    outputFile = paste0(outputFolder, "/", exonID, "/", exonID, 
                        ".", replicateID, ".pdf")
    if (!file.exists(paste0(outputFolder, "/", exonID))) {
      dir.create(paste0(outputFolder, "/", exonID), recursive = T)
    }
    
    if (dataName == "MAVE") {
      plotData(data, outputFile)
    }
    
    # add the replicate ID to the end of the sample_id
    data$sampleID2 = paste0(data$sampleID2, "_", replicateID)
    dataAll[[i]] = data
  }
    
  # make it a wide format
  spliceColumnMax = max(grep("SpliceAI", colnames(data)))
  columnsToUse = c(colnames(data)[1:spliceColumnMax], "SNVID", "sampleID2")
  
  exons = fileInfo[, 2]
  uniqueExons = unique(exons)
  dataWide = list()
  for (j in 1:length(uniqueExons)) {
    exonID = uniqueExons[j]
    data = NULL
    index = which(exons == exonID)
    for (i in 1:length(index)) {
      data = rbind(data, dataAll[[ index[i] ]])
    }
    
    # convert data to a wide format
    dataToWide = data[, columnsToUse]
    wide = dataToWide %>% pivot_wider(names_from = sampleID2, 
                                          values_from = EventCount)
    dataWide[[exonID]] = wide
  }
  
  
  return(dataWide)
}

compareResult = function() {
  dataFolder = "/research/bsi/projects/breast/s108235.tripneg_Couch/projects/Wenan/crispr/data/Drug/analysis/E17/"
  defaultFile = paste0(dataFolder, "/E17_DESeq2_all_results.default.csv")
  synonymousFile = paste0(dataFolder, "/E17_DESeq2_all_results.Synonymous.csv")
  
  result1 = read.csv(defaultFile)
  result2 = read.csv(synonymousFile)
  
  # make sure SNVID are identical
  identical(result1$SNVID, result2$SNVID)
  
  pdf(paste0(dataFolder, "/pvalue_comparison.pdf"), onefile = T)
  par(mfrow = c(2, 2))
  
  metrics = c("log2FoldChange", "pvalue", "padj")
  for (j in 1:length(metrics)) {
    columns = paste0(metrics[j], c(".treat_D5_vs_X0nM", ".treat_X5nM_vs_X0nM", 
                ".treat_X25nM_vs_X0nM", ".treat_X100nM_vs_X0nM"))
    for (i in 1:length(columns)) {
      if (j == 1) {
        x = result2[, columns[i]]
        y = result1[, columns[i]]
      } else {
        x = -log10(result2[, columns[i]])
        y = -log10(result1[, columns[i]])
      }
      plot(x, y, 
           xlab = "synonymous", ylab = "default", main = columns[i])
      abline(a = 0, b =1, col = "red")
    }
  }
  
  dev.off()
  
}

main = function() {
  # identify the outliers from the drug experiment
  
  args = commandArgs(T)
  
  inputFileList = args[1]
  outputFolder = args[2]
  dataName = args[3]
  SNVOnly = as.logical(args[4])
  
  # # for MAVE
  # refType = "Synonymous"
  # # dataName = "MAVE"
  
  # for MiniLib
  refType = "default"
  # refType = "Synonymous"
  SNVOnly = F
  dataName = "MiniLib"
  inputFileList="/research/bsi/projects/breast/s108235.tripneg_Couch/projects/Wenan/crispr/data/Hypomorph_mini_drug/analysis/input.txt"
  outputFolder="/research/bsi/projects/breast/s108235.tripneg_Couch/projects/Wenan/crispr/data/Hypomorph_mini_drug/analysis/"
  
  if (!file.exists(outputFolder)) {
    dir.create(outputFolder, recursive = T)
  }
  
  # browser()
  data = loadData(inputFileList, outputFolder, dataName, SNVOnly)
  
  exonIDAll = names(data)
  for (j in 1:length(exonIDAll)) {
    exonID = exonIDAll[j]
    print(exonID)
    
    if (!file.exists(paste0(outputFolder, "/", exonID))) {
      dir.create(paste0(outputFolder, "/", exonID), recursive = T)
    }
    
    # browser()
    # run DE analysis
    dataDE = data[[exonID]]
    indexStart = grep("SNVID", colnames(dataDE))
    libIndex = grep("lib_", colnames(dataDE))
    if (refType != "default") {
      indexRefType = which(dataDE$EventType == refType)
    } else {
      indexRefType = NULL
    }
    index = setdiff(indexStart:dim(dataDE)[2], libIndex)
    
    dataToUse = data.frame(dataDE[, index])

    # remove rows with NA
    countData = dataToUse[, -1]
    toRemove = which(rowSums(is.na(countData)) > 0)
    if (length(toRemove) > 0) {
      dataDE = dataDE[-toRemove, ]
      dataToUse = dataToUse[-toRemove, ]
    }
    
    # prepare data
    countData = dataToUse[, -1]
    rownames(countData) = dataToUse[, 1]
    
    countColumns = colnames(countData)
    replicateID = factor(sub(".*_", "", countColumns))
    treat = sub("_.*", "", countColumns)
    
    if (dataName == "MAVE") {
      treat = factor(treat, levels = c("X0nM", "D5", "X1nM", 
                                     "X5nM", "X25nM", "X100nM"))
    } else if (dataName == "MiniLib") {
      treat = factor(treat, levels = c("D8", "D5", "Olaparib"))
    }
    
    colData = data.frame(treat, replicateID)
    
    # use synonymous variants as the library size
    librarySize = NULL
    if (!is.null(indexRefType)) {
      librarySize = colSums(countData[indexRefType, ])
    } 
    
    library(DESeq2)
    sizeFactor = NULL
    if (!is.null(librarySize)) {
      if (is.matrix(librarySize)) {
        sizeFactor = librarySize / exp(rowMeans(log(librarySize)))
      } else {
        sizeFactor = librarySize / exp(mean(log(librarySize)))
        # expand to a matrix
        sizeFactor = matrix(c(sizeFactor), byrow = T, 
                            nrow = dim(countData)[1], ncol = dim(countData)[2])
      }
    }
    
    dds = DESeqDataSetFromMatrix(countData = countData,
                                 colData = colData,
                                 design = ~ treat + replicateID)
    if (!is.null(sizeFactor)) {
      normalizationFactors(dds) = sizeFactor
    }
    dds = DESeq(dds)
    
    # output all comparisons
    contrasts = grep("treat_", resultsNames(dds), value = T)
    resultAll = NULL
    for (i in 1:length(contrasts)) {
      resultData = results(dds, name=contrasts[i])
      colnames(resultData) = paste0(colnames(resultData), ".", contrasts[i])
      if (i > 1) {
        resultAll = cbind(resultAll, resultData[, -1])
      } else {
        resultAll = resultData
      }
    }
    
    # plot the candidate hypomorphic variants 
    pdf(paste0(outputFolder, "/", exonID, "/", exonID, "_DESeq2_all_results.",  
               refType, "_SNVOnly", SNVOnly, ".log2FC.pdf"))
    
    if (dataName == "MAVE") {
      conditions = c(".treat_X1nM_vs_X0nM", ".treat_X5nM_vs_X0nM", 
                     ".treat_X25nM_vs_X0nM", ".treat_X100nM_vs_X0nM")
      names = c("1nM", "5nM", "25nM", "100nM")
      par(mfrow = c(2, 2))
    } else {
      conditions = c(".treat_Olaparib_vs_D8")
      names = c("Olaparib")
      par(mfrow = c(1, 1))
    }
    
    for (i in 1:length(conditions)) {
      # shape
      indexSig = (resultAll[, paste0("padj", conditions[i])] < 0.05)
      pch = rep(1, dim(resultAll)[1])
      pch[indexSig] = 3
      
      # color
      type = dataDE$EventType
      color = rep("gray", dim(resultAll)[1])
      color[type == "Synonymous"] = "green"
      color[type == "StopGain"] = "red"
      color[type == "Missense"] = "blue"
      
      if (dataName == "MAVE") {
        x = -resultAll[, "log2FoldChange.treat_D5_vs_X0nM"]
        y = resultAll[, paste0("log2FoldChange", conditions[i])]
        plot(x, y, col = color, pch = pch, xlab = "log2FC D10 vs D5 (no drug)", 
             ylab = "log2FC D10 with vs without drug", 
             main = names[i])
      } else if (dataName == "MiniLib") {
        x = -resultAll[, "log2FoldChange.treat_D5_vs_D8"]
        y = resultAll[, paste0("log2FoldChange", conditions[i])]
        plot(x, y, col = color, pch = pch, xlab = "log2FC D8 vs D5 (no drug)", 
             ylab = "log2FC D8 with vs without drug", 
             main = names[i])
      }
    }
    dev.off()
    
    dataOut = cbind(dataDE, resultAll)
    write.csv(as.data.frame(dataOut), 
              file = paste0(outputFolder, "/", exonID, "/", 
                            exonID, "_DESeq2_all_results.",  
                            refType, "_SNVOnly", SNVOnly, ".csv"))
  }
}

main()