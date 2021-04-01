message(paste("Start:",Sys.time(),sep=" "))
suppressMessages(library(DMRcaller))
library(betareg)
library(RColorBrewer)
args = commandArgs(trailingOnly = T)
thread = as.numeric(args[1])
args[-1] -> args
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
files <- paste(rownames(group),"CX_report","txt","gz",sep=".")
message("Loading...")
for(i in 1:p){
  assign(rownames(group)[i], readBismark(files[i]))
}
cmd <- paste(paste(rownames(group),'=',rownames(group),sep=""),collapse = ',')
methylationDataList <- eval(parse(text=paste("GRangesList(",cmd,")",sep="")))
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
    DMRsReplicatesRegions <- computeDMRsReplicates(compareList, condition = condition, context = context[l], method = "bins", minProportionDifference = diff[l], cores = thread)
    write.csv(DMRsReplicatesRegions, paste(genotype[j],'vs',genotype[1],context[l],'csv',sep="."), quote = FALSE)
  }
}
message(paste("End:",Sys.time(),sep=" "))
