options(stringsAsFactors = FALSE)

##Functions for cluster frequency analysis
getTotalCellNumbers = function(file, wd) {
    setwd(wd)
    table = read.csv(file, header=TRUE)
    files = list.files(pattern = "*clustered.txt$")
    indecies = rep.int(0, length(files))
    
    for (sample in 1:length(table[,1])) {
        sampleName = table[sample,1]
        sampleName = strsplit(sampleName, ".fcs")
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
        newFile = read.table(f, header = TRUE, stringsAsFactors = FALSE)
        popSizes = as.matrix(newFile$popsize)
        return(popSizes)
    }
    
    popSizeMatrix = apply(clusteredFiles, 1, extractPopSizes)
    
    if (totalCellNumbers[1] == "NA") {
        freqMatrix = t(t(popSizeMatrix) / colSums(popSizeMatrix))
    } else {
        freqMatrix = t(t(popSizeMatrix)/totalCellNumbers)
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
    
    row_without_values = apply(fullFreqMatrix, 1, function(row) all(row==0))
    row_without_values_indecies = which(row_without_values == TRUE)
    
    if (length(row_without_values_indecies) > 0) {
        freqMatrix = fullFreqMatrix[-row_without_values_indecies,]
    } else {
        freqMatrix = fullFreqMatrix
    }
    
    samResults = SAM(x=freqMatrix,y=sampleID,resp.type=family,
                     genenames=rownames(freqMatrix), geneid=rownames(freqMatrix), nperms=nperms)
    return(samResults)
}



##This funciton checks each cluster to see if it is in the list of 
##higher or lower populations and saves the transformed q-values as a matrix
getSignifMatrix = function(model, num_clusters, qValue_cutoff) {
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
    
    ##If there is only one significant cluster, SAM will return a vector. 
    ##If there are multiple, SAM will return a matrix.
    if (is.null(model$siggenes.table$genes.up)) {
        pop_lower = c()
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
    
    
    cluster_signif = mat.or.vec(num_clusters, 2)
    colnames(cluster_signif) = c("cellType", "frequencySignif")
    cluster_signif[,1] = 1:num_clusters
    
    clusters = as.matrix(1:num_clusters)
    cluster_signif[,2] = apply(clusters, 1, getClusterSignif)
    return(cluster_signif)
}

##This function appends the q-value matrix to each "clustered.txt" file in the wd.
append_freq_signif = function(wd, colNames, signif_matrix) {
    setwd(wd)
    files = as.matrix(list.files(pattern = "*clustered.txt$"))
    
    bind_signif_by_file = function(newFile) {
        currFile = read.table(newFile, header = TRUE, sep = "\t")
        colnames(currFile) = as.matrix(colNames)
        appendedFile = cbind(currFile, signif_matrix[,2])
        colnames(appendedFile)[grep("signif_matrix", colnames(appendedFile))] = "frequencySignif"
        write.table(appendedFile, file = newFile, row.names = F, sep = "\t", quote = F)
    }
    
    apply(files, 1, bind_signif_by_file)
    
}

##Run the cluster frequency significance analysis 
analyze_cluster_frequencies = function(wd, group1, group2, qValue_cutoff, nperms, total_cells_in_file, total_cell_numbers_csv) {
    totalCellNumbers = c()
    
    if (total_cells_in_file == TRUE) {
        totalCellNumbers = "NA"
    } else {
        totalCellNumbers = getTotalCellNumbers(file = paste(working.directory, total_cell_numbers_csv, sep = "/"))
    }
    
    freqMatrix = getFrequencies(wd = wd, totalCellNumbers = totalCellNumbers)
    
    sampleID = getSampleIDs(group1, group2, wd)
    
    model = run_SAM_analysis(freqMatrix, sampleID = sampleID, nperms = nperms)
    
    signif_matrix = getSignifMatrix(model, num_clusters = length(freqMatrix[,1]), qValue_cutoff = qValue_cutoff)
    write.csv(signif_matrix, paste(wd, "freqSignif.csv", sep="/"), row.names = FALSE)
    colNames = extractColnames(wd)
    append_freq_signif(wd, colNames = colNames, signif_matrix = signif_matrix)
    return(list.files(pattern = "*.clustered.txt$"))
}

##Remove the frequency significance column
remove_freqsignif_columns = function(wd) {
    options(stringsAsFactors = F)
    
    setwd(wd)
    
    files = list.files(pattern = "*clustered.txt$")
    origColNames = as.vector(read.table(files[1], nrows = 1))
    signifColumns = grep(pattern="frequencySignif", origColNames)
    
    if (length(signifColumns) > 0) {
        newColNames = origColNames[-signifColumns]
        
        for (i in 1:length(files)) {
            currFile = read.table(files[i], header=TRUE)
            newFile = currFile[,-signifColumns]
            colnames(newFile) = newColNames
            write.table(newFile, file = files[i], row.names = F, sep = "\t", quote = F)
        }
    } 
    return("Frequency significance columns removed.")
}


##Functions for Boolean experession analysis
get_num_clusters = function(wd){
    files = list.files(pattern = "*clustered.txt$")
    testFile = read.table(files[1], header=TRUE)
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

append_expr_signif = function(wd, colNames, signif_matrix, feature) {
    setwd(wd)
    files = as.matrix(list.files(pattern = "*clustered.txt$"))
    
    bind_signif_by_file = function(newFile) {
        currFile = read.table(newFile, header = TRUE, sep = "\t")
        colnames(currFile) = as.matrix(colNames)
        appendedFile = cbind(currFile, signif_matrix[,2])
        colnames(appendedFile)[grep("signif_matrix", colnames(appendedFile))] = paste(feature, "BooleanSignif", sep="")
        write.table(appendedFile, file = newFile, row.names = F, sep = "\t", quote = F)
    }
    
    apply(files, 1, bind_signif_by_file)
}


analyze_cluster_expression = function(wd, group1, group2, qValue_cutoff, nperms, feature, booleanThreshold, asinh.cofactor) {
    setwd(wd)
    num_clusters = get_num_clusters(wd)
    booleanMatrix = buildBooleanMatrix(wd, feature = feature, 
                                       booleanThreshold = booleanThreshold, 
                                       asinh.cofactor = asinh.cofactor, 
                                       num_clusters = num_clusters)
    write.csv(booleanMatrix, paste(wd, paste(feature, "BooleanFreq.csv", sep=""), sep="/"), row.names = TRUE)
    sampleID = getSampleIDs(group1, group2, wd)
    model = run_SAM_analysis(booleanMatrix, sampleID = sampleID, nperms = nperms)
    signif_matrix = getSignifMatrix(model, num_clusters = num_clusters, qValue_cutoff = qValue_cutoff)
    write.csv(signif_matrix, paste(wd, paste(feature, "BooleanSignif.csv", sep=""), sep="/"), row.names = FALSE)
    
    colNames = extractColnames(wd)
    append_expr_signif(wd, colNames = colNames, signif_matrix = signif_matrix, feature = feature)
    return(list.files(pattern = "*.clustered.txt$"))
}

##Remove the expression significance columns
remove_exprsignif_columns = function(wd) {
    options(stringsAsFactors = F)
    
    setwd(wd)
    
    files = list.files(pattern = "*clustered.txt$")
    origColNames = as.vector(read.table(files[1], nrows = 1))
    signifColumns = grep(pattern="Signif", origColNames)
    
    if (any(grepl(pattern="frequencySignif", origColNames))) {
        signifColumns = signifColumns[-which(signifColumns == grep(pattern="frequencySignif", origColNames))]
    }
    
    
    if (length(signifColumns) > 0) {
        newColNames = origColNames[-signifColumns]
        
        for (i in 1:length(files)) {
            currFile = read.table(files[i], header=TRUE)
            newFile = currFile[,-signifColumns]
            colnames(newFile) = newColNames
            write.table(newFile, file = files[i], row.names = F, sep = "\t", quote = F)
        }
    } 
    return("Expression significance columns removed.")
}