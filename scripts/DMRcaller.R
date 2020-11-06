message(paste("Start:",Sys.time(),sep=" "))
suppressMessages(library(DMRcaller))
library(betareg)
library(RColorBrewer)
args = commandArgs(trailingOnly = T)
thread = as.numeric(args[1])
binsize = as.numeric(args[2])
args[-c(1:2)] -> args
args[c(TRUE,FALSE)] -> genotype
args[c(FALSE,TRUE)] -> replicates
p <- sum(as.numeric(replicates))
data.frame(rep("x",p), stringsAsFactors = FALSE) -> group
p <- 0
for (n in 1:length(genotype)){
  for (k in 1:replicates[n]){
    p = p + 1
    paste(genotype[n],k,sep="_") -> rownames(group)[p]
    genotype[n] -> group[p,1]
  }
}
colnames(group) <-"Type"
diff <- c(0.4,0.2,0.1)
context <- c("CG","CHG","CHH")
files <- paste(rownames(group),"CX_report","txt",sep=".")
message("Loading...")
for(i in 1:p){
  assign(rownames(group)[i], readBismark(files[i]))
}
cmd <- paste(paste(rownames(group),'=',rownames(group),sep=""),collapse = ',')
methylationDataList <- eval(parse(text=paste("GRangesList(",cmd,")",sep="")))
pdf("MethylCProfile.pdf",28,4)
for(i in 1:p){
  message(paste("Plotting", rownames(group)[i], "...", sep=" "))
  plotMethylationProfileFromData(methylationDataList[[i]], conditionsNames = rownames(group)[i], windowSize = 10000, autoscale = FALSE, context = context, labels = LETTERS)
}
dev.off()
control = which(group[,1]==genotype[1])
controlList = methylationDataList[[control[1]]]
for(k in 2:length(control)){
  controlList = joinReplicates(controlList, methylationDataList[[control[k]]], usecomplete = TRUE)
}
for(j in 2:length(genotype)){
  message(paste("Compare",genotype[j],"and",genotype[1],"...",sep=" "))
  treatment = which(group[,1]==genotype[j])
  condition = c(group[group[,1]==genotype[1],],group[group[,1]==genotype[j],])
  compareList <- joinReplicates(controlList, methylationDataList[[treatment[1]]], usecomplete = TRUE)
  for(k in 2:length(treatment)){
    compareList <- joinReplicates(compareList, methylationDataList[[treatment[k]]], usecomplete = TRUE)
  }
  for(l in 1:3){
    message(paste("Calculating",context[l],"...",sep=" "))
    DMRsReplicatesRegions <- computeDMRsReplicates(compareList, condition = condition, context = context[l], method = "bins", minProportionDifference = diff[l], binSize = binsize, cores = thread)
    write.csv(DMRsReplicatesRegions, paste(genotype[j],'vs',genotype[1],context[l],'csv',sep="."), quote = FALSE)
  }
}
message(paste("End:",Sys.time(),sep=" "))
