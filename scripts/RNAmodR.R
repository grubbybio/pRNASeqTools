options(warn=-1)
args = commandArgs(trailingOnly = T)
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(RNAmodR.RiboMethSeq)
  library(RNAmodR.Data)
})
options(ucscChromosomeNames=FALSE)
sequences <- paste(args[1], "fa", sep=".")
annotation <- GFF3File(paste(args[1], "gff", sep="."))
file <- c(treated = paste(args[1],"filtered","bam",sep="."))
mrms <- ModRiboMethSeq(file, annotation = annotation, sequences = sequences)
settings(mrms) <- list(minScoreMean = 0.93, minScoreRMS = 0.5, minCoverage = as.integer(args[2]))
modify(mrms,force = TRUE) -> mrms
sink(paste(args[1],"RMS","wig",sep="."))
for(i in 1:length(ranges(mrms))){
  cat(paste("variableStep chrom=",as.character(seqnames(ranges(mrms)))[i],"\n",sep=""))
  values(getDataTrack(mrms,i,type = c("scoreRMS", "scoreMean"))$scoreRMS) -> tmp
  for(j in 1:width(ranges(ranges(mrms))[[i]])){
    if(tmp[,j]>0){
      cat(paste(j,"\t",tmp[,j],"\n",sep=""))
    }
  }
}
sink()
sink(paste(args[1],"Mean","wig",sep="."))
for(i in 1:length(ranges(mrms))){
  cat(paste("variableStep chrom=",as.character(seqnames(ranges(mrms)))[i],"\n",sep=""))
  values(getDataTrack(mrms,i,type = c("scoreRMS", "scoreMean"))$scoreMean) -> tmp
  for(j in 1:width(ranges(ranges(mrms))[[i]])){
    if(tmp[,j]>0){
      cat(paste(j,"\t",tmp[,j],"\n",sep=""))
    }
  }
}
sink()
