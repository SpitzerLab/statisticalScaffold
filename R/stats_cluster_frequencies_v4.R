options(stringsAsFactors = FALSE)

##Functions for cluster frequency analysis
getTotalCellNumbers = function(file, wd) {
  setwd(wd)
  table = read.csv(file, header=TRUE)
  files = list.files(pattern = "*.clustered.txt$")
  indecies = rep.int(0, length(files))
  
  for (sample in 1:length(table[,1])) {
    sampleName = table[sample,1]
    if (!(all(grepl(sampleName, files)==FALSE))) {
      order = which(grepl(sampleName, files)==TRUE)
      indecies[order] = sample
    }
  }
  
  if (!all(indecies != 0)) {
    print("There is a missing index")
  }
  
  totalCellNumbers = table[indecies,2]
  return(totalCellNumbers)
}

getSampleIDs = function(group1, group2, wd) {
  setwd(wd)
  files = list.files(pattern = "*clustered.txt$")
  sampleID = rep.int(0, length(files))
  
  for (variable in 1:length(group1)){
    group = group1[variable]
    indecies = grep(group, files)
    sampleID[indecies] = 1
  }
  
  for (variable in 1:length(group2)){
    group = group2[variable]
    indecies = grep(group, files)
    sampleID[indecies] = 2
  }
  
  if (!all(sampleID!=0)) {
    print("There is an unassigned file")
  }
  
  return(sampleID)
}

getFrequencies = function(wd, totalCellNumbers) {
  setwd(wd)
  clusteredFiles = as.matrix(list.files(pattern = "*.clustered.txt$"))
  
  extractPopSizes = function(f) {
    newFile = read.table(f, header = TRUE, stringsAsFactors = FALSE, sep="\t")
    popSizes = as.matrix(newFile$popsize)
    return(popSizes)
  }
  
  popSizeMatrix = apply(clusteredFiles, 1, extractPopSizes)
  
  if (totalCellNumbers[1] == "NA") {
    freqMatrix = t(t(popSizeMatrix) / colSums(popSizeMatrix))
  } else {
    freqMatrix = t(t(popSizeMatrix) / totalCellNumbers)
  }
  
  rownames(freqMatrix) = 1:length(freqMatrix[,1])
  colnames(freqMatrix) = clusteredFiles
  write.csv(freqMatrix, file = paste(wd, "frequencyMatrix.csv", sep="/"))
  return(freqMatrix)
}


##This function extracts the column names from the first "clustered.txt" file and 
##stores them. This is to avoid changes caused during read-write to punctuation.
extractColnames = function(wd) {
  return(read.table(list.files(pattern = "*clustered.txt$")[1], nrows = 1, sep = "\t"))
}


##This funciton reads in all files in the working directory, extracts the frequency for
##each cluster and stores this in a matrix. Then this matrix is passed into SAM along 
##with the sampleID codes to build the model of different features by sample type.
run_SAM_analysis = function(freqMatrix, sampleID, nperms) {
  if (length(unique(sampleID)) == 2) {
    family = "Two class unpaired"
  } else {
    family = "Multiclass"
  }
  
  fullFreqMatrix = as.data.frame(freqMatrix)
  
  row_without_values = apply(fullFreqMatrix, 1, function(row) all(row==0 || row==1))
  row_without_values_indecies = which(row_without_values == TRUE)
  
  if (length(row_without_values_indecies) > 0) {
    freqMatrix = fullFreqMatrix[-row_without_values_indecies,]
  } else {
    freqMatrix = fullFreqMatrix
  }
  
  freqMatrix = as.matrix(freqMatrix) ## Fix for new version of R
  samResults = SAM(x=freqMatrix,y=sampleID,resp.type=family,
                   genenames=rownames(freqMatrix), geneid=rownames(freqMatrix), nperms=nperms)
  return(samResults)
}

##This function manually calculates the fold chnage between groups for plotting
##regardless of significance cutoff             FIND ME
calculate_FoldChange = function(freqMatrix, group1, group2) {
  result_FC <- numeric()
  
  group1_mean <- rowMeans(freqMatrix[,which(grepl(group1, colnames(freqMatrix)))])
  group2_mean <- rowMeans(freqMatrix[,which(grepl(group2, colnames(freqMatrix)))])
  # 
  # for(a in 1:length(group1_mean)){
  #   if(group1_mean[a]|group2_mean[a] == 0) {group1_mean[a] = group1_mean[a] + 1/200; group2_mean[a] = group2_mean[a] + 1/200}
  # }
  
  result_FC <- group2_mean/group1_mean
  result_FC[is.na(result_FC)] <- 1; result_FC[is.infinite(result_FC)] <- 1
  
  result_FC <- round(log2(result_FC), digits = 3)

  return(result_FC)
}


##This funciton checks each cluster to see if it is in the list of 
##higher or lower populations and saves the transformed q-values as a matrix
getSignifMatrix = function(model, num_clusters, qValue_cutoff, include_foldC) {
  getClusterSignif = function(cluster) {
    if (is.element(cluster, pop_higher$`Gene Name`)) {
      entry = which(pop_higher$`Gene Name` == cluster)
      if (as.numeric(as.character(pop_higher$`q-value(%)`[entry])) <= qValue_cutoff) {
        return((100 - as.numeric(as.character(pop_higher$`q-value(%)`[entry]))) / 100)
      } else {
        return(0.5)
      }
    } else if (is.element(cluster, pop_lower$`Gene Name`)) {
      entry = which(pop_lower$`Gene Name` == cluster)
      if (as.numeric(as.character(pop_lower$`q-value(%)`[entry])) <= qValue_cutoff) {
        return((as.numeric(as.character(pop_lower$`q-value(%)`[entry]))) / 100)
      } else {
        return(0.5)
      }
    } else {
      return(0.5)
    }
  }
  
  getClusterFoldChange = function(cluster) {
    if (is.element(cluster, pop_higher$`Gene Name`)) {
      entry = which(pop_higher$`Gene Name` == cluster)
      if (as.numeric(as.character(pop_higher$`q-value(%)`[entry])) <= qValue_cutoff) {
        return(as.numeric(as.character(pop_higher$`Fold Change`[entry]))) 
      } else {
        return(1)
      }
    } else if (is.element(cluster, pop_lower$`Gene Name`)) {
      entry = which(pop_lower$`Gene Name` == cluster)
      if (as.numeric(as.character(pop_lower$`q-value(%)`[entry])) <= qValue_cutoff) {
        return(as.numeric(as.character(pop_lower$`Fold Change`[entry])))
      } else {
        return(1)
      }
    } else {
      return(1)
    }
  }
  
  ##If there is only one significant cluster, SAM will return a vector. 
  ##If there are multiple, SAM will return a matrix.
  if (is.null(model$siggenes.table$genes.up)) {
    pop_higher = c()
  } else if (is.vector(model$siggenes.table$genes.up)) {
    pop_higher = as.data.frame(t(model$siggenes.table$genes.up), stringsAsFactors = FALSE)
  } else if (is.matrix(model$siggenes.table$genes.up)) {
    pop_higher = as.data.frame(model$siggenes.table$genes.up, stringsAsFactors = FALSE)
  }
  
  if (is.null(model$siggenes.table$genes.lo)) {
    pop_lower = c()
  } else if (is.vector(model$siggenes.table$genes.lo)) {
    pop_lower = as.data.frame(t(model$siggenes.table$genes.lo), stringsAsFactors = FALSE)
  } else if (is.matrix(model$siggenes.table$genes.lo)) {
    pop_lower = as.data.frame(model$siggenes.table$genes.lo, stringsAsFactors = FALSE)
  }
  
  
  cluster_signif = mat.or.vec(num_clusters, 3)
  colnames(cluster_signif) = c("cellType", "Significance", "SAMFoldChange")
  cluster_signif[,1] = 1:num_clusters
  
  clusters = as.matrix(1:num_clusters)
  cluster_signif[,2] = apply(clusters, 1, getClusterSignif)
  if(include_foldC) {
    cluster_signif[,3] = apply(clusters, 1, getClusterFoldChange)
    if(any(is.infinite(cluster_signif[,3]))) {cluster_signif[which(is.infinite(cluster_signif[,3])),3] = 1}
  } else {cluster_signif = cluster_signif[,-3]}
    
  return(cluster_signif)
}

##This function appends the q-value matrix to each "clustered.txt" file in the wd.
append_freq_signif = function(wd, colNames, signif_matrix, include_foldC = FALSE, set_max_Val = 2) {
  setwd(wd)
  files = as.matrix(list.files(pattern = "*clustered.txt$"))
  
  bind_signif_by_file = function(newFile) {
    currFile = read.table(newFile, header = TRUE, sep = "\t")
    colnames(currFile) = as.matrix(colNames)
    
    if(include_foldC) {
      appendedFile = cbind(currFile, signif_matrix[,-1])
      colnames(appendedFile)[grep("Significance", colnames(appendedFile))] = "frequencySignif"
      colnames(appendedFile)[grep("SAMFoldChange", colnames(appendedFile))] = "frequencyFoldChange"
      colnames(appendedFile)[grep("All_Log2_FoldChange", colnames(appendedFile))] = "frequency_ALLFoldChange"
      
      #Log2 Normalize SAM output
      appendedFile$frequencyFoldChange = log2(as.numeric(appendedFile$frequencyFoldChange))
      
      # Normalize Fold Changes
        set_min_Val = -set_max_Val
        FCs <- c("frequencyFoldChange","frequency_ALLFoldChange")
        for (i in 1:2){
          FC_dat <- appendedFile[,which(grepl(FCs[i], colnames(appendedFile)))]
          
          if(any(FC_dat > set_max_Val)) {FC_dat[which(FC_dat > set_max_Val)] = set_max_Val}
          if(any(FC_dat < set_min_Val)) {FC_dat[which(FC_dat < set_min_Val)] = set_min_Val}
          normalized = (FC_dat-set_min_Val)/(set_max_Val-set_min_Val)
        
          if(any(normalized>1)) {normalized[which(normalized>1)] = 0.99999} ##represents greatest increase, 1 used for black landmarks
          if(any(normalized==1)) {normalized[which(normalized==1)] = 0.99999}
        
          appendedFile[,which(grepl(FCs[i], colnames(appendedFile)))] = normalized
        }
    } else{
      appendedFile = cbind(currFile, signif_matrix[,grep("Significance", colnames(signif_matrix))])
      colnames(appendedFile)[ncol(appendedFile)] = "frequencySignif"
    }
    
    write.table(appendedFile, file = newFile, row.names = F, sep = "\t", quote = F)
  }
  
  apply(files, 1, bind_signif_by_file)
  
}

##Run the cluster frequency significance analysis 
analyze_cluster_frequencies = function(wd, group1, group2, qValue_cutoff, nperms, total_cells_in_file, total_cell_numbers_csv, include_foldC,set_max_Val) {
  totalCellNumbers = c()
  
  if (total_cells_in_file == TRUE) {
    totalCellNumbers = "NA"
  } else {
    totalCellNumbers = getTotalCellNumbers(file = paste(wd, total_cell_numbers_csv, sep = "/"), wd = wd)
  }
  
  freqMatrix = getFrequencies(wd = wd, totalCellNumbers = totalCellNumbers)
  
  sampleID = getSampleIDs(group1, group2, wd)
  
  model = run_SAM_analysis(freqMatrix, sampleID = sampleID, nperms = nperms)
  
  getFC <- calculate_FoldChange(freqMatrix, group1, group2)
  
  signif_matrix = getSignifMatrix(model, num_clusters = length(freqMatrix[,1]), qValue_cutoff = qValue_cutoff, include_foldC)
    colNames = extractColnames(wd)
  
  signif_matrix <- as.data.frame(cbind(signif_matrix,as.numeric(getFC)))
    colnames(signif_matrix)[ncol(signif_matrix)] = "All_Log2_FoldChange"
  append_freq_signif(wd, colNames = colNames, signif_matrix = signif_matrix, include_foldC = include_foldC, set_max_Val = set_max_Val)
  
  write.csv(signif_matrix, paste(wd, "freqSignif.csv", sep="/"), row.names = FALSE)
  
  return(list.files(pattern = "*.clustered.txt$"))
}


remove_freqsignif_columns = function(wd) {
  options(stringsAsFactors = F)
  
  setwd(wd)
  
  files = list.files(pattern = "*clustered.txt$")
  for (i in 1:length(files)) {
    origColNames = as.vector(read.table(files[i], nrows = 1, sep="\t"))
    signifColumns = c(grep(pattern="frequency", origColNames))
    
    if (length(signifColumns) > 0) {
      newColNames = origColNames[-signifColumns]
      currFile = read.table(files[i], header=TRUE, sep="\t")
      newFile = currFile[,-signifColumns]
      colnames(newFile) = newColNames
      write.table(newFile, file = files[i], row.names = F, sep = "\t", quote = F)
    }
  } 
  return("Frequency significance columns removed.")
}


print_Feature_Template = function(wd){
  setwd(wd)
  files = list.files(pattern = "*.clustered.txt$")
  temp_mat = data.frame(fileName = files, feature_value = rep(NA, length = length(files)))
  write.csv(temp_mat, paste(wd,"Feature_Template.csv", sep="/"), row.names = FALSE)
}



##Functions for Boolean experession analysis
get_num_clusters = function(wd){
  files = list.files(pattern = "*clustered.txt$")
  testFile = read.table(files[1], header=TRUE, sep="\t")
  return(length(testFile[,1]))
}

get_rdata_col_names <- function(working.directory, f.name)
{
  rdata.file <- readRDS(paste(working.directory, f.name, sep = "/"))
  
  ret <- as.vector(colnames(rdata.file))

  if(any(is.na(ret)))
  {
    w <- is.na(ret)
    ret[w] <- as.vector(w)
  }
  
  return(ret)
}

buildBooleanMatrix = function(wd, feature, booleanThreshold, asinh.cofactor, num_clusters){
  setwd(wd)
  dataFiles = as.matrix(list.files(pattern = "*.RData$"))
  booleanMatrix = mat.or.vec(num_clusters, length(dataFiles))
  
  populateBooleanMatrixColumn = function(dataFile) {
    currFile = readRDS(file = dataFile)
    getBooleanFrequency = function(cluster) {
      currVec = currFile[which(currFile$cellType == cluster), which(colnames(currFile) == feature)]
      booleanFrequency = as.vector(sum(currVec >= asinh(booleanThreshold / asinh.cofactor)) / length(currVec))
      return(booleanFrequency)
    }
    
    clusters = as.matrix(1:num_clusters)
    booleanFreqVec = apply(clusters, 1, getBooleanFrequency)
    return(booleanFreqVec)
  }
  
  fullBooleanMatrix = apply(dataFiles, 1, populateBooleanMatrixColumn)
  rownames(fullBooleanMatrix) = 1:num_clusters
  colnames(fullBooleanMatrix) = dataFiles
  
  #write.csv(fullBooleanMatrix, paste(wd, paste(feature, "BooleanFreq_FULL.csv", sep=""), sep="/"), row.names = TRUE)
  
  row_without_values = apply(fullBooleanMatrix, 1, function(row) all(is.na(row)))
  row_without_values_indecies = which(row_without_values == TRUE)
  
  
  if (length(row_without_values_indecies) > 0) {
    culledBooleanMatrix = fullBooleanMatrix[-row_without_values_indecies,]
  } else {
    culledBooleanMatrix = fullBooleanMatrix
  }
  
  culledBooleanMatrix[is.na(culledBooleanMatrix)] = 0
  
  row_all_zero = apply(culledBooleanMatrix, 1, function(row) all(row==0))
  row_all_zero_indecies = which(row_all_zero == TRUE)
  
  
  if (length(row_all_zero_indecies) > 0) {
    booleanMatrix = culledBooleanMatrix[-row_all_zero_indecies,]
  } else {
    booleanMatrix = culledBooleanMatrix
  }
  
  
  return(booleanMatrix)
}

append_expr_signif = function(wd, colNames, signif_matrix, feature, include_foldC = FALSE, set_max_Val = 2) {
  setwd(wd)
  files = as.matrix(list.files(pattern = "*clustered.txt$"))
  
  bind_signif_by_file = function(newFile) {
    currFile = read.table(newFile, header = TRUE, sep = "\t")
    colnames(currFile) = as.matrix(colNames)
    
    if(include_foldC) {
      appendedFile = cbind(currFile, signif_matrix[,-1])
      colnames(appendedFile)[grep("Significance", colnames(appendedFile))] = paste(feature, "BooleanSignif", sep="")
      colnames(appendedFile)[grep("SAMFoldChange", colnames(appendedFile))] = paste(feature, "BooleanFoldChange", sep="")
      colnames(appendedFile)[grep("All_Log2_FoldChange", colnames(appendedFile))] = paste(feature, "Boolean_ALLFoldChange", sep="")
        
      #Log2 Normalize SAM output
      SAM_Index <- which(grepl(paste(feature, "Boolean_FoldChange", sep=""), colnames(appendedFile)))
        appendedFile[,SAM_Index] = log2(as.numeric(appendedFile[,SAM_Index]))
      
      # Normalize Fold Changes
      set_min_Val = -set_max_Val
      FCs <- c(paste(feature, "Boolean_FoldChange", sep=""),paste(feature, "Boolean_ALLFoldChange", sep=""))
      for (i in 1:2){
        FC_dat <- appendedFile[,which(grepl(FCs[i], colnames(appendedFile)))]
        
        if(any(FC_dat > set_max_Val)) {FC_dat[which(FC_dat > set_max_Val)] = set_max_Val}
        if(any(FC_dat < set_min_Val)) {FC_dat[which(FC_dat < set_min_Val)] = set_min_Val}
        normalized = (FC_dat-set_min_Val)/(set_max_Val-set_min_Val)
        
        if(any(normalized>1)) {normalized[which(normalized>1)] = 0.99999} ##represents greatest increase, 1 used for black landmarks
        if(any(normalized==1)) {normalized[which(normalized==1)] = 0.99999}
        
        appendedFile[,which(grepl(FCs[i], colnames(appendedFile)))] = normalized
      }
    } else{
      appendedFile = cbind(currFile, signif_matrix[,grep("Significance", colnames(signif_matrix))])
      colnames(appendedFile)[ncol(appendedFile)] = paste(feature, "BooleanSignif", sep="")
    }
    
    write.table(appendedFile, file = newFile, row.names = F, sep = "\t", quote = F)
  }
  
  apply(files, 1, bind_signif_by_file)
}


analyze_cluster_expression = function(wd, group1, group2, qValue_cutoff, nperms, feature, booleanThreshold, 
                                      asinh.cofactor, include_foldC, set_max_Val) {
  setwd(wd)
  num_clusters = get_num_clusters(wd)
  booleanMatrix = buildBooleanMatrix(wd, feature = feature, 
                                     booleanThreshold = booleanThreshold, 
                                     asinh.cofactor = asinh.cofactor, 
                                     num_clusters = num_clusters)
  write.csv(booleanMatrix, paste(wd, paste(feature, "BooleanFreq.csv", sep=""), sep="/"), row.names = TRUE)
  sampleID = getSampleIDs(group1, group2, wd)
  model = run_SAM_analysis(booleanMatrix, sampleID = sampleID, nperms = nperms)
  
  getFC_express <- calculate_FoldChange(booleanMatrix, group1, group2)
  
  signif_matrix = getSignifMatrix(model, num_clusters = num_clusters, qValue_cutoff = qValue_cutoff, include_foldC)
    signif_matrix = cbind(signif_matrix, getFC_express)
    colnames(signif_matrix)[ncol(signif_matrix)] = "All_Log2_FoldChange"
  
  write.csv(signif_matrix, paste(wd, paste(feature, "BooleanSignif.csv", sep=""), sep="/"), row.names = FALSE)
  colNames = extractColnames(wd)
  append_expr_signif(wd, colNames = colNames, signif_matrix = signif_matrix, feature = feature, 
                     include_foldC = include_foldC, set_max_Val = set_max_Val)
  return(list.files(pattern = "*.clustered.txt$"))
}


remove_exprsignif_columns = function(wd) {
  options(stringsAsFactors = F)
  
  setwd(wd)
  
  files = list.files(pattern = "*clustered.txt$")
  for (i in 1:length(files)) {
    origColNames = as.vector(read.table(files[i], nrows = 1, sep="\t"))
    signifColumns = c(grep(pattern="Boolean", origColNames))
    
    if (length(signifColumns) > 0) {
      newColNames = origColNames[-signifColumns]
      currFile = read.table(files[i], header=TRUE, sep="\t")
      newFile = currFile[,-signifColumns]
      colnames(newFile) = newColNames
      write.table(newFile, file = files[i], row.names = F, sep = "\t", quote = F)
    }
  }
  return("Expression significance columns removed.")
}


########
## Added 20190.10.10 to plot features correlated with cluster abundances

analyze_cluster_correlation = function(wd, corFeature, corTest, corPlotTypes, qVal_cutoff) {
  setwd(wd)
  
  num_clusters = get_num_clusters(wd)
  corr_matrix = data.frame(cellType = 1:num_clusters, Correlation = rep(NA, length = num_clusters), pVal = rep(NA, length = num_clusters))
  
  freqMatrix = getFrequencies(wd = wd, totalCellNumbers = "NA")
    freqMatrix = freqMatrix[,order(colnames(freqMatrix))]
  
  corr_featureMatrix = read.csv(paste(wd, corFeature,sep="/"), stringsAsFactors = FALSE)
    corr_featureMatrix = corr_featureMatrix[order(corr_featureMatrix$fileName),]
    
  if(all(colnames(freqMatrix) == corr_featureMatrix$fileName)){
    feature_values <- corr_featureMatrix$feature_value
    
    for(i in 1:num_clusters) {
      this_clus = freqMatrix[i,]
      corr_result = cor.test(feature_values, this_clus, method = corTest, exact = FALSE)
      
      corr_matrix$Correlation[i] = corr_result$estimate
      corr_matrix$pVal[i] = corr_result$p.value
      if(is.na(corr_matrix$Correlation[i])) {
        corr_matrix$Correlation[i] = 0
        corr_matrix$pVal[i] = 1
      }
    }
    corr_matrix$pVal <- p.adjust(corr_matrix$pVal, method = "hochberg")
  
    colNames = extractColnames(wd)
    
        corfeatureName <- as.character(corFeature); corfeatureName <- substr(".csv","",corfeatureName)

    append_corr_signif(wd, corfeatureName, colNames = colNames, corr_matrix, corPlotTypes, qVal_cutoff)
  
    write.csv(corr_matrix, paste(wd,"/",corfeatureName,"_correlationSignif.csv", sep=""), row.names = FALSE)

    return(list.files(pattern = "*.clustered.txt$"))
  } else {
    return("Error, please include a feature value for each file in your directory")
  }
}


append_corr_signif = function(wd, corfeatureName, colNames, corr_matrix, corPlotTypes, qVal_cutoff) {
  setwd(wd)
  files = as.matrix(list.files(pattern = "*clustered.txt$"))
  
  bind_signif_by_file = function(newFile) {
    currFile = read.table(newFile, header = TRUE, sep = "\t")
    colnames(currFile) = as.matrix(colNames)
    
    significant_Corrs = corr_matrix$Correlation
      significant_Corrs[which(corr_matrix$pVal > (qVal_cutoff/100))] = 0
    
    appendedFile = cbind(currFile, significant_Corrs)
      colnames(appendedFile)[ncol(appendedFile)] <- paste(corfeatureName,"correlationSIG",sep="_")
    
    if(corPlotTypes){
      appendedFile = cbind(appendedFile, corr_matrix$Correlation)
        colnames(appendedFile)[ncol(appendedFile)] <- paste(corfeatureName,"correlationALL",sep="_")
    }
      
    write.table(appendedFile, file = newFile, row.names = F, sep = "\t", quote = F)
  }
  
  apply(files, 1, bind_signif_by_file)
}


remove_correlation_columns = function(wd) {
  options(stringsAsFactors = F)
  
  setwd(wd)
  
  files = list.files(pattern = "*clustered.txt$")
  for (i in 1:length(files)) {
    origColNames = as.vector(read.table(files[i], nrows = 1, sep="\t"))
    signifColumns = c(grep(pattern="correlation", origColNames))
    
    if (length(signifColumns) > 0) {
      newColNames = origColNames[-signifColumns]
      currFile = read.table(files[i], header=TRUE, sep="\t")
      newFile = currFile[,-signifColumns]
      colnames(newFile) = newColNames
      write.table(newFile, file = files[i], row.names = F, sep = "\t", quote = F)
    }
  }
  return("Expression significance columns removed.")
}

