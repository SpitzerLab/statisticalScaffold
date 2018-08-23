options(stringsAsFactors = F)


get_stats_by_sample <- function(tab) 
{
    tab.medians <- ddply(tab, ~cellType, colwise(median, is.numeric))
    tab.medians.by.sample <- ddply(tab, ~cellType * sample, colwise(median, is.numeric))
    pop.size <- ddply(tab, ~cellType, nrow)
    names(pop.size) <- gsub("V1", "popsize", names(pop.size))
    pop.size.by.sample <- ddply(tab, ~cellType * sample, nrow)
    names(pop.size.by.sample) <- gsub("V1", "popsize", names(pop.size.by.sample))
    tab.medians <- merge(tab.medians, pop.size, by = "cellType")
    tab.medians.by.sample <- merge(tab.medians.by.sample, pop.size.by.sample, by = c("cellType", "sample"), all.x = T)
    
    
    #Rotate the by.sample table
    temp <- melt(tab.medians.by.sample, id = c("cellType", "sample"))       
    temp$variable <- paste(temp$variable, temp$sample, sep = "@")
    temp$sample <- NULL
    temp <- cast(temp, cellType~variable)
    
    
    ret <- merge(tab.medians, temp, by = "cellType", all.x = T)
    
    return(ret)
    
}

process_files_groups <- function(files, wd, col.names, num_clusters, num_samples, asinh.cofactor)
{
    setwd(wd)
    
    cluster_data <- function(tab, col.names, k, algorithm = "", ...)
    {
        m <- as.matrix(tab[, col.names])
        if(algorithm == "clara")
        {
            print("Performing clara clustering")
            groups <- clara(m, k, ...)$clustering
        }
        
        else if(algorithm == "hierarchical")
        {
            print("Performing hierarchical clustering")
            dend <- hclust(dist(m), ...)
            groups <- cutree(dend, k)
        }
        
        print("Clustering done")
        tab <- cbind(tab, groups)
        return(tab)
    }
    tab <- NULL
    orig.data <- NULL
    
    for(f in files)
    {
        print(f)
        fcs.file <- read.FCS(f)
        temp.orig.data <- exprs(fcs.file)
        temp.tab <- convert_fcs(fcs.file, asinh.cofactor)
        colnames(temp.tab) <- pData(parameters(fcs.file))$desc
        
        if(any(is.na(colnames(temp.tab))))
        {
            w <- is.na(colnames(temp.tab))
            colnames(temp.tab)[w] <- pData(parameters(fcs.file))$name[w]
        }
        
        temp.tab <- as.matrix(temp.tab)
        temp.tab[temp.tab < 0] <- 0
        temp.tab <- as.data.frame(temp.tab)
        
        temp.tab <- data.frame(temp.tab, sample = f)
        temp.orig.data <- data.frame(temp.orig.data, sample = f)
        tab <- rbind(tab, temp.tab)
        orig.data <- rbind(orig.data, temp.orig.data)
    }
    
    m <- cluster_data(tab, col.names, k = num_clusters, algorithm = "clara", sampsize = min(nrow(tab), 1000), samples = num_samples)
    colnames(m) <- gsub("groups", "cellType", colnames(m))
    orig.data <- cbind(orig.data, cellType = m[, "cellType"])
    
    #tab.medians <- ddply(m, ~cellType * sample, colwise(median))
    #pop.size <- ddply(m, ~cellType * sample, nrow)
    #temp <- data.frame(tab.medians, popsize = pop.size[tab.medians$cellType, "V1"], check.names = F, stringsAsFactors = FALSE)
    
    temp <- get_stats_by_sample(m)
    
    colnames(temp) <- gsub("^X", "", colnames(temp))
    m <- data.frame(m, check.names = F)
    orig.data <- data.frame(orig.data, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(orig.data) <- gsub("^X", "", colnames(orig.data))
    colnames(m) <- gsub("^X", "", colnames(m))
    
    write.table(temp, paste(f, ".clustered.txt", sep = ""), row.names = F, sep = "\t", quote = F)
    my_save(m, paste(f, ".clustered.all_events.RData", sep = ""))
    #my_save(orig.data, paste(f, ".clustered.all_events.orig_data.RData", sep = ""))
}



process_file <- function(f, wd, col.names, num_clusters, num_samples, asinh.cofactor)
{
    setwd(wd)
    
    cluster_data <- function(tab, col.names, k, algorithm = "", ...)
    {
        m <- as.matrix(tab[, col.names])
        
        if(algorithm == "clara")
        {
            print("Performing clara clustering")
            groups <- clara(m, k, ...)$clustering
        }
        
        else if(algorithm == "hierarchical")
        {
            print("Performing hierarchical clustering")
            dend <- hclust(dist(m), ...)
            groups <- cutree(dend, k)
        }
        
        print("Clustering done")
        tab <- cbind(tab, groups)
        return(tab)
    }
    
    fcs.file <- read.FCS(f)
    orig.data <- exprs(fcs.file)
    tab <- convert_fcs(fcs.file, asinh.cofactor)
    colnames(tab) <- pData(parameters(fcs.file))$desc

    if(any(is.na(colnames(tab))))
    {
        w <- is.na(colnames(tab))
        colnames(tab)[w] <- pData(parameters(fcs.file))$name[w]
    }
    
    
    
    tab <- as.matrix(tab)
    tab[tab < 0] <- 0
    tab <- as.data.frame(tab)
    
    

    m <- cluster_data(tab, col.names, k = num_clusters, algorithm = "clara", sampsize = min(nrow(tab), 1000), samples = num_samples)
    colnames(m) <- gsub("groups", "cellType", colnames(m))
    orig.data <- cbind(orig.data, cellType = m[, "cellType"])
    
    tab.medians <- ddply(m, ~cellType, colwise(median))
    
    pop.size <- ddply(m, ~cellType, nrow)
    
    temp <- data.frame(tab.medians, sample = f, popsize = pop.size[tab.medians$cellType, "V1"], check.names = F, stringsAsFactors = FALSE)
    
    colnames(temp) <- gsub("^X", "", colnames(temp))
    m <- data.frame(m, check.names = F)
    orig.data <- data.frame(orig.data, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(orig.data) <- gsub("^X", "", colnames(orig.data))
    colnames(m) <- gsub("^X", "", colnames(m))
    
    write.table(temp, paste(f, ".clustered.txt", sep = ""), row.names = F, sep = "\t", quote = F)
    my_save(m, paste(f, ".clustered.all_events.RData", sep = ""))
    #my_save(orig.data, paste(f, ".clustered.all_events.orig_data.RData", sep = ""))
}


##Testing cluster files together 9/24/16
process_files_together <- function(f, wd, col.names, num_clusters, num_samples, asinh.cofactor, fileIDs, totalCellNumbers)
{
    setwd(wd)
    
    cluster_data <- function(tab, col.names, k, algorithm = "", ...)
    {
        m <- as.matrix(tab[, col.names])
        #Ensure that m is a numeric matrix
        class(m) = "numeric"
        
        if(algorithm == "clara")
        {
            print("Performing clara clustering")
            groups <- clara(m, k, ...)$clustering
        }
        
        else if(algorithm == "hierarchical")
        {
            print("Performing hierarchical clustering")
            dend <- hclust(dist(m), ...)
            groups <- cutree(dend, k)
        }
        
        print("Clustering done")
        tab <- cbind(tab, groups)
        return(tab)
    }
    
    tab = as.data.frame(c())
    orig.data = c()
    
    addNextFile = function(fileName) {
        fcs.file <- read.FCS(fileName)
        orig.data <- exprs(fcs.file)
        tabNextFile <- convert_fcs(fcs.file, asinh.cofactor)
        colnames(tabNextFile) <- pData(parameters(fcs.file))$desc
        
        tabNextFile <- as.matrix(tabNextFile)
        tabNextFile[tabNextFile < 0] <- 0
        
        fileID = rep.int(fileName, length(tabNextFile[,1]))
        tabNextFile = as.data.frame(cbind(tabNextFile, fileID))
        
        tab = rbind(tab, tabNextFile)
        return(tab)
    }
    
    f = as.matrix(f)
    tab = apply(f, 1, addNextFile)
    tab = as.data.frame(do.call("rbind", tab))
    
    m <- cluster_data(tab, col.names, k = num_clusters, algorithm = "clara", sampsize = min(nrow(tab), 1000), samples = num_samples)
    colnames(m) <- gsub("groups", "cellType", colnames(m))
    orig.data <- cbind(orig.data, cellType = m[, "cellType"])
    
    #Compute the medians once on the whole data file. 
    #This must be done because the coordinates control the mapping.
    m2 = as.matrix(m[,which(colnames(m) != "fileID")])
    class(m2) = "numeric"
    m2 = as.data.frame(m2)
    tab.medians <- ddply(m2, ~cellType, colwise(median))
    
    
    
    write_medians_and_files = function(fileID) {
        workingFrame = m[which(m$fileID == fileID),]
        populations = as.matrix(1:num_clusters)
        
        getPopSize = function(pop) {
            return(nrow(workingFrame[which(workingFrame$cellType == pop),]))
        }
        
        pop.size <- as.matrix(apply(populations, 1, getPopSize))
        
        
        temp <- data.frame(tab.medians, sample = fileID, popsize = pop.size, check.names = F, stringsAsFactors = FALSE)
        temp$popsize[is.na(temp$popsize)] = 0
        
        colnames(temp) <- gsub("^X", "", colnames(temp))
        workingFrame = as.matrix(workingFrame[,which(colnames(workingFrame) != "fileID")])
        class(workingFrame) = "numeric"
        workingFrame <- data.frame(workingFrame, check.names = F)
        orig.data <- data.frame(orig.data, stringsAsFactors = FALSE, check.names = FALSE)
        colnames(orig.data) <- gsub("^X", "", colnames(orig.data))
        colnames(workingFrame) <- gsub("^X", "", colnames(workingFrame))
        
        write.table(temp, paste(fileID, ".clustered.txt", sep = ""), row.names = F, sep = "\t", quote = F)
        my_save(workingFrame, paste(fileID, ".clustered.all_events.RData", sep = ""))
        #my_save(orig.data, paste(f, ".clustered.all_events.orig_data.RData", sep = ""))
    }
    
    ##Write medians and clustered.txt files for each sample
    f = as.matrix(f)
    apply(f, 1, write_medians_and_files)
    
}


# cluster_fcs_files_in_dir <- function(wd, num.cores, col.names, num_clusters, num_samples, asinh.cofactor)
# {
#     files.list <- list.files(path = wd, pattern = "*.fcs$")
#     parallel::mclapply(files.list, mc.cores = num.cores, mc.preschedule = FALSE,
#              process_file, wd = wd, col.names = col.names, num_clusters = num_clusters, num_samples = num_samples, asinh.cofactor = asinh.cofactor)
#     return(files.list)
# }

##Testing here for cluster together 9/24/16
cluster_fcs_files_in_dir <- function(wd, num.cores, col.names, num_clusters, num_samples, asinh.cofactor, cluster.together)
{
    files.list <- list.files(path = wd, pattern = "*.fcs$")
    
    if (cluster.together) {
        process_files_together(f = files.list, wd = wd, col.names = col.names, num_clusters = num_clusters, num_samples = num_samples, asinh.cofactor = asinh.cofactor, fileIDs = files.list, totalCellNumbers = totalCellNumbers)
        return(files.list)
    } else {
        parallel::mclapply(files.list, mc.cores = num.cores, mc.preschedule = FALSE,
                           process_file, wd = wd, col.names = col.names, num_clusters = num_clusters, num_samples = num_samples, asinh.cofactor = asinh.cofactor)
        return(files.list)
    }
}


cluster_fcs_files_groups <- function(wd, files.list, num.cores, col.names, num_clusters, num_samples, asinh.cofactor)
{
    lapply(files.list,
                       process_files_groups, wd = wd, col.names = col.names, num_clusters = num_clusters, num_samples = num_samples, asinh.cofactor = asinh.cofactor)
    return(files.list)
}



