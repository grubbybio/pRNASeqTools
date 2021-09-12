options(warn=-1)
args = commandArgs(trailingOnly = T)
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(RNAmodR.RiboMethSeq)
})
options(ucscChromosomeNames=FALSE)
sequences <- paste(args[1], "fa", sep=".")
annotation <- GFF3File(paste(args[1], "gff", sep="."))
file <- c(treated = paste(args[1],"filtered","bam",sep="."))
mrms <- ModRiboMethSeq(file, annotation = annotation, sequences = sequences, minScoreMean = 0.93, minScoreA = 0.5, minCoverage = as.integer(args[2]),flankingRegion=as.integer(2))
write.csv(mrms@modifications,paste(args[1],"csv",sep = "."), quote = F, row.names = F)
sink(paste(args[1],"ScoreA2","wig",sep="."))
for(i in 1:length(ranges(mrms))){
  cat(paste("variableStep chrom=",as.character(seqnames(ranges(mrms)))[i],"\n",sep=""))
  values(getDataTrack(mrms,i,type = c("scoreA", "scoreMean", "scoreRMS"))$scoreA) -> tmp
  for(j in 1:width(ranges(ranges(mrms))[[i]])){
    if(tmp[,j]>0){
      cat(paste(j,"\t",tmp[,j],"\n",sep=""))
    }
  }
}
sink()
sink(paste(args[1],"ScoreMean","wig",sep="."))
for(i in 1:length(ranges(mrms))){
  cat(paste("variableStep chrom=",as.character(seqnames(ranges(mrms)))[i],"\n",sep=""))
  values(getDataTrack(mrms,i,type = c("scoreA", "scoreMean", "scoreRMS"))$scoreMean) -> tmp
  for(j in 1:width(ranges(ranges(mrms))[[i]])){
    if(tmp[,j]>0){
      cat(paste(j,"\t",tmp[,j],"\n",sep=""))
    }
  }
}
sink()
sink(paste(args[1],"ScoreRMS","wig",sep="."))
for(i in 1:length(ranges(mrms))){
  cat(paste("variableStep chrom=",as.character(seqnames(ranges(mrms)))[i],"\n",sep=""))
  values(getDataTrack(mrms,i,type = c("scoreA", "scoreMean", "scoreRMS"))$scoreRMS) -> tmp
  for(j in 1:width(ranges(ranges(mrms))[[i]])){
    if(tmp[,j]>0){
      cat(paste(j,"\t",tmp[,j],"\n",sep=""))
    }
  }
}
sink()