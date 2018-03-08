### This function takes in user defined vectors of clusters, Group 1 and Group 2 for which
  # they want to compare the magnitude of change in protein expression, defined by histogram
  # intersection distance. 
  # The user must also define a new directory name & protein list to include

####
generateHistIntersect <- function (nameDirectory, working.directory, stat = TRUE, proteins) {
  codeTracker = timestamp()
  print("Starting Histogram Intersection Distance")
  
  if(nameDirectory == "Histogram Intersection") {
    info = gsub(" ","_",gsub(":","-",Sys.time()))
    newDirectory = paste(info,nameDirectory,sep = " ")
    file.rename(paste(working.directory,nameDirectory,sep = "/"), paste(working.directory,newDirectory,sep = "/"))
    nameDirectory = newDirectory
  }
  
  setwd(working.directory)
  Rfiles = as.matrix(list.files(pattern = ".RData"))
  v1 = read.csv(paste(working.directory,nameDirectory,"Vector1.csv",sep = "/"))
    vector1 = v1$x
  v2 = read.csv(paste(working.directory,nameDirectory,"Vector2.csv",sep = "/"))
    vector2 = v2$x

  getCells <- function(currentFile) {
    dat = as.data.frame(readRDS(currentFile))
    g1Dat = dat[which(dat$cellType %in% vector1),]
    g2Dat = dat[which(dat$cellType %in% vector2),]
    g1Dat$Sample = rep(as.character(currentFile), length=nrow(g1Dat))
    g2Dat$Sample = rep(as.character(currentFile), length=nrow(g2Dat))
    bigDat = rbind(g1Dat,g2Dat)
    return(bigDat)
  }
  bigList = apply(Rfiles,1,getCells)
  bigMatrix = do.call(rbind,bigList)

  dat1 = bigMatrix[which(bigMatrix$cellType %in% vector1),]
  dat2 = bigMatrix[which(bigMatrix$cellType %in% vector2),]

  remove_dat1 = which(!(colnames(dat1)[1:55] %in% proteins))
    if(length(remove_dat1) > 0) {dat1 = dat1[,-remove_dat1]}
  remove_dat2 = which(!(colnames(dat2)[1:55] %in% proteins))
    if(length(remove_dat2) > 0) {dat2 = dat2[,-remove_dat2]}
    if("beadDist" %in% colnames(dat1)) {
      dat1 = dat1[,-which(colnames(dat1) == "beadDist")]
      dat2 = dat2[,-which(colnames(dat2) == "beadDist")]
    }
  
  # Order & Untransform data  
  dat1 = dat1[,c(ncol(dat1),(ncol(dat1)-1),1:(ncol(dat1)-2))]
  dat2 = dat2[,c(ncol(dat2),(ncol(dat2)-1),1:(ncol(dat2)-2))]
  dat1[,3:ncol(dat1)] = sinh(dat1[,3:ncol(dat1)]) * 5
  dat2[,3:ncol(dat2)] = sinh(dat2[,3:ncol(dat2)]) * 5
  
  codeTracker = append(codeTracker, "Samples in G1:")
    codeTracker = append(codeTracker,unique(dat1$Sample))
    codeTracker = append(codeTracker, paste("Number of cells in G1:",as.character(nrow(dat1))))
  codeTracker = append(codeTracker,"Samples in G2:")
    codeTracker = append(codeTracker,unique(dat2$Sample))
    codeTracker = append(codeTracker, paste("Number of cells in G2:",as.character(nrow(dat2))))
  
  if(!(all(dat1$Sample %in% dat2$Sample) && all(dat2$Sample %in% dat1$Sample))) {
    codeTracker = append(codeTracker, "Different samples represented in G1 vs G2 - some samples have empty clusters.")
    warning("BA: Different samples represented in G1 vs G2 - some samples have empty clusters.")
  }
  
  if(unique(dat1$Sample) < 2 || unique(dat2$Sample) < 2) {
    codeTracker = append(codeTracker,"Less than 2 samples per group - very rare populations. Cannot run statistics.")
    stat = FALSE
    warning("BA: Less than 2 samples per group - very rare populations. Cannot run statistics.")
  }
  
    
# Initializing MoC functions
  
  # MASTER FUNCTION to run Magnitude of Change
  expressChange <- function(dat1, dat2, runName = "CurrentCells", codeTracker,stat) {

    if(nrow(dat1) < 50 || nrow(dat2) < 50) {
      codeTracker = append(codeTracker,"Less than 50 cells per group. Cannot run statistics.")
      stat = FALSE
      warning("BA: Less than 50 cells per group. Cannot run statistics.")
    } else if(nrow(dat1) < 10 || nrow(dat2) < 10) {
      stop("Cluster(s) too rare for comparison.")
    }

    print("Starting MOC")
    # # # Initializes required variables
    Protein = names(dat1[3:ncol(dat1)])
    output = as.data.frame(Protein)
    Percent_Diff = data.frame()
    outCorrectHID = data.frame()
    Median1 = numeric()
    Median2= numeric()
    HistInt = numeric()
    
    # # # Calculates Magnitude of Change by histogram intersection distance for each protein
    for(i in 3:length(dat1)) {
      g1 = dat1[[i]]
      g2 = dat2[[i]]
      
      Median1 = append(Median1, median(g1))
      Median2 = append(Median2, median(g2))
      
      if (max(g1) > max(g2)){setBreak = as.integer(max(g1))+1
      }else{setBreak = as.integer(max(g2))+1}
      g1.hist = hist(g1, breaks = seq(0,setBreak, l=100), plot = F)
      g2.hist = hist(g2, breaks = seq(0,setBreak, l=100), plot = F)
      HistInt = append(HistInt, Hist.Intersect.Dist(g1.hist, g2.hist)) # BA Histogram Intersection distance Master Dependency 1 *
    }
    if (stat == TRUE) {HisInt_SD = calculateSD(dat1,dat2) # Master Dependency 2 **
    }else {HisInt_SD = rep(0, length = nrow(output))}
    
    if(all(HisInt_SD == 0)) {codeTracker = append(codeTracker, "SD calculation failed - Too few cells per mouse")}
    
    
    # # # Runs correction, iterating over subsets of dataset #1
    for(k in 1:100) {
      correct1 = sample_n(dat1, size = nrow(dat1)/4)
      correct2 = sample_n(dat1, size = nrow(dat1)/4)
      HistInt_corr = numeric()
      
      for(m in 3:length(correct1)) {
        cor1 = correct1[[m]]
        cor2 = correct2[[m]]
        
        if (max(cor1) > max(cor2)){setBreak = as.integer(max(cor1))+1
        }else{setBreak = as.integer(max(cor2))+1}
        cor1.hist = hist(cor1, breaks = seq(0,setBreak, l=100), plot = F)
        cor2.hist = hist(cor2, breaks = seq(0,setBreak, l=100), plot = F)
        HistInt_corr = append(HistInt_corr, Hist.Intersect.Dist(cor1.hist, cor2.hist))
      }
      if(k == 1) {outCorrectHID = as.data.frame(HistInt_corr)
      } else {outCorrectHID = rowMeans(cbind(outCorrectHID, as.data.frame(HistInt_corr)))}
    }
    
    # # # Runs SAM analysis
    print(stat)
    if(stat == TRUE) {SAMresult = runStats(dat1,dat2,runName)} # Master Dependency 3 ***
    
    print("Generating Outputs...")
    # # # Generating output files
    output = cbind(output, Group1 = Median1, Group2 = Median2)
    
    # Calculating Percent Diff
    output$Percent_Diff = NA
    for (i in 1:nrow(output)){
      if(output$Group1[i] == 0 || output$Group2[i] == 0) {output$Percent_Diff[i] = NA
      } else {output$Percent_Diff[i] = (output$Group2[i] / output$Group1[i])*100}
    }
    # Categorizing protein expression change
    output$TypeOfChange = NA
    if(stat == TRUE) {output$TypeOfChange[which(is.na(SAMresult$SAMqval))] = "Unchanged"}
    output$TypeOfChange[which(output$Group1 < 2 & output$Group2 < 2)] = "Not Expressed"
    #output$TypeOfChange[which(is.na(output$TypeOfChange) & output$Group1 < 2 & output$Group2 >= 2)] = "Turns ON"
    #output$TypeOfChange[which(is.na(output$TypeOfChange) & output$Group1 >= 2 & output$Group2 < 2)] = "Turns OFF"
    output$TypeOfChange[which(is.na(output$TypeOfChange) & output$Group1 > output$Group2)] = "Decreases_in_G2"
    output$TypeOfChange[which(is.na(output$TypeOfChange) & output$Group1 < output$Group2)] = "Increases_in_G2"
    
    if(stat == TRUE) {output = cbind(output, SAMresult)}
    output = cbind(output, HistInt_Dist = HistInt, HistInt_Correct = outCorrectHID, HisInt_SD)
    output = output[order(output$HistInt_Dist, decreasing = T),]
    output = do.call(data.frame,lapply(output, function(x) replace(x, is.infinite(x),NA)))
    write.csv(output, file = paste(runName,"SummaryOfChanges.csv", sep="/"),
              row.names=FALSE)
   
    # # # Generating summary plots
    return(plotChange(output, runName, codeTracker, stat)) # Master Dependency 4 ****
  }
  
  
# Master Dependency 1 *
### BA version of Histogram Intersection Distance, relying on percent per bin not absolute count
Hist.Intersect.Dist <- function (h1,h2) {
  stopifnot(inherits(h1, "histogram"), inherits(h2, "histogram"))
  stopifnot(all(h1$breaks == h2$breaks))
  h1.count = (h1$counts/sum(h1$counts))*100
  h2.count = (h2$counts/sum(h2$counts))*100
  intersect.counts <- pmin(h1.count, h2.count)
  h3 = BuildHistogram(h1$breaks, intersect.counts)
  return(1 - (sum(h3$counts) / sum(h2.count)))	
}

BuildHistogram <- function(breaks, counts, xname="") {
  # Returns a histogram object from the given list of breaks and counts.
  
  stopifnot(is.numeric(breaks), is.numeric(counts))
  stopifnot(length(breaks) > 1)
  stopifnot(length(breaks) == (length(counts) + 1))
  hist <- list(breaks = breaks,
               counts = counts,
               density = counts / (sum(counts) * diff(breaks)),
               mids = (head(breaks, -1) + tail(breaks, -1)) / 2,
               xname = xname,
               equidist = BreaksAreEquidistant(breaks))
  class(hist) <- "histogram"
  return(hist)
}

BreaksAreEquidistant <- function(breaks) {
  # Check if breaks are equally spaced.
  diffs <- diff(breaks)
  all(abs(diffs - diffs[1]) < .Machine$double.eps^0.5 * max(diffs))
}

# Master Dependency 2 **
### Calculates Standard Div of histogram intersection distance between mice
calculateSD = function(dat1,dat2) {
  Mouse_SD_HID = as.data.frame(colnames(dat1)[3:ncol(dat1)])
  colnames(Mouse_SD_HID) = "Protein"
  HID_iterations = Mouse_SD_HID
  files = as.character(unique(dat1$Sample))
    if(!(all(files %in% dat2$Sample))) {files = files[-which(!(files %in% dat2$Sample))]}
  
  for(i in 1:length(files)){
    setA = as.data.frame(dat1[which(dat1$Sample == files[i]),3:ncol(dat1)])
    setB = as.data.frame(dat2[which(dat2$Sample == files[i]),3:ncol(dat2)])
    
    if(nrow(setA) > 100 && nrow(setA) > 100) {
        Mouse_HID = numeric()
        for(j in 1:length(setA)) {
          A = setA[[j]]
          B = setB[[j]]
          
          if (max(A) > max(B)){setBreak = as.integer(max(A))+1
          }else{setBreak = as.integer(max(B))+1}
          A.hist = hist(A, breaks = seq(0,setBreak, l=100), plot = F)
          B.hist = hist(B, breaks = seq(0,setBreak, l=100), plot = F)	
          Mouse_HID = append(Mouse_HID, Hist.Intersect.Dist(A.hist, B.hist))
        }
        Mouse_SD_HID$newMouse = Mouse_HID
        colnames(Mouse_SD_HID)[ncol(Mouse_SD_HID)] = files[i]
    }
  }	
  if (ncol(Mouse_SD_HID) > 2) {
    HisInt_SD = apply(Mouse_SD_HID[,2:length(Mouse_SD_HID)],1, function(x) sd(x))
  } else {
    HisInt_SD = rep(0, length = nrow(Mouse_SD_HID))
    warning("BA: SD calculation failed - Too few cells per mouse")
  }
  return(HisInt_SD)
}


# Master Dependency 3 ***
### Runs SAM analysis on each protein across samples
# # # internal dependencies
samTable = function(dat1, dat2) {
  sumStats = data.frame(matrix(nrow=ncol(dat1)-2))
  rownames(sumStats) = names(dat1[3:ncol(dat1)])
  
  files = unique(dat1$Sample)
  for(i in 1:length(files)) {
    x = dat1[which(dat1$Sample == files[i]), 3:ncol(dat1)]
    sumStats[i] = apply(x, 2, function(x) median(x))
    names(sumStats)[i] = paste(files[i], "G1", sep ="_")
  }
  
  for(i in 1:length(files)) {	
    y = dat2[which(dat2$Sample == files[i]), 3:ncol(dat2)]
    sumStats = cbind(sumStats, apply(y, 2, function(x) median(x)))
    names(sumStats)[ncol(sumStats)] = paste(files[i], "G2", sep ="_")
  }
  return(sumStats)
}

getSampleIDs = function(names) {
  sampleID = rep.int(0, length(names))
  sampleID[grep("G1", names)] = 1
  sampleID[grep("G2", names)] = 2
  if (!all(sampleID!=0)) {print("There is an unassigned file")}
  return(sampleID)
}

# # # Execution
runStats = function(dat1, dat2, runName) {
  sumStats = samTable(dat1, dat2)
  write.csv(sumStats, file = paste(runName,"MedianSignalIntensity.csv", sep="/"))
  sampleID = getSampleIDs(colnames(sumStats))
  fullMatrix = as.data.frame(sumStats)
  
  row_without_values = apply(fullMatrix, 1, function(row) all(row==0))
  row_without_values_indecies = which(row_without_values == TRUE)
  if (length(row_without_values_indecies) > 0) {statMatrix = fullMatrix[-row_without_values_indecies,]
  } else {statMatrix = fullMatrix}
  
  samResults = SAM(x=statMatrix,y=sampleID,resp.type="Two class unpaired",
                   genenames=rownames(statMatrix), geneid=rownames(statMatrix), nperms=10000)
  
  qStat = as.data.frame(rep(NA, length=nrow(sumStats)))
  rownames(qStat) = rownames(sumStats)
  colnames(qStat) = "SAMqval"
  qStat$SAMfoldChange = rep(NA, length=nrow(sumStats))
  UP = as.data.frame(samResults$siggenes.table$genes.up, stringsAsFactors = FALSE)
  DOWN = as.data.frame(samResults$siggenes.table$genes.lo, stringsAsFactors = FALSE)
  for(p in 1:nrow(UP)) {
    qStat$SAMqval[which(rownames(qStat) == UP$`Gene ID`[p])] = UP$`q-value(%)`[p]
    qStat$SAMfoldChange[which(rownames(qStat) == UP$`Gene ID`[p])] = UP$`Fold Change`[p]
  }
  for(p in 1:nrow(DOWN)) {
    qStat$SAMqval[which(rownames(qStat) == DOWN$`Gene ID`[p])] = DOWN$`q-value(%)`[p]
    qStat$SAMfoldChange[which(rownames(qStat) == DOWN$`Gene ID`[p])] = DOWN$`Fold Change`[p]
  }
  return(qStat)
}


# Master Dependency 4 ****
### plotChange takes the Summary output file and runName variable from the expressChange function
# Returns a barplot saved as PDF to the working directory that displays magnitude of change by histogram intersection distance
plotChange <- function(result, runName, codeTracker,stat) {
  plotDat = as.data.frame(result$Protein)
  colnames(plotDat) = "Protein"
  
  # Gets relevant parameters
  plotDat = cbind(plotDat, value = result$HistInt_Dist, sd = result$HisInt_SD, NullDist = result$HistInt_Correct,
                  G2_Expression_Change = result$TypeOfChange)
  
  # Cleaning the data
  if(any(plotDat$G2_Expression_Change == "Unchanged")) {plotDat = plotDat[-which(plotDat$G2_Expression_Change == "Unchanged"),]}
  if(any(plotDat$G2_Expression_Change == "Not Expressed")) {plotDat = plotDat[-which(plotDat$G2_Expression_Change == "Not Expressed"),]}
  plotDat$value[which(plotDat$G2_Expression_Change == "Decreases_in_G2")] = plotDat$value[which(plotDat$G2_Expression_Change == "Decreases_in_G2")]*-1
  #plotDat$value[which(plotDat$G2_Expression_Change == "Turns OFF")] = plotDat$value[which(plotDat$G2_Expression_Change == "Turns OFF")]*-1
  
  graphTitle = paste(runName,"Not Tested For Significance",sep="__")
  if(stat == TRUE) {graphTitle = paste("Significant",graphTitle,sep="__")}
  
  plt = ggplot(plotDat, aes(x = reorder(Protein, -value), y = value, fill = G2_Expression_Change)) + 
    geom_bar(stat="identity") + theme_classic() + 
    scale_fill_manual(values = c("Increases_in_G2"="red","Decreases_in_G2"="blue"))+
    theme(axis.text.x = element_text(angle = 90)) + xlab("") + scale_y_continuous(limits=c(-1,1)) +
    ggtitle(graphTitle)+ ylab(paste("Magnitude of Change: Histogram Intersection Distance")) + ggtitle(runName) +
    geom_hline(yintercept=median(plotDat$NullDist), linetype="dashed", color = "darkgrey") +
    geom_hline(yintercept=-median(plotDat$NullDist), linetype="dashed", color = "darkgrey")
  if(!all(plotDat$sd == 0)) {
    plt = plt + geom_errorbar(stat="identity", position="identity",aes(ymin = value-sd, ymax = value+sd))
  }
  ggsave(filename = paste(runName,"Hist-Intersect-Dist.pdf", sep="/"),
         plot=plt,width=11,height=8.5)
  
  codeTracker = append(codeTracker, timestamp())
  write.table(codeTracker,file = paste(runName,"ProcessingNotes.txt", sep="/"),row.names=FALSE)
  return(plt)
}

  ##### Calling MoC
  return(expressChange(dat1,dat2,nameDirectory,codeTracker,stat))
}

