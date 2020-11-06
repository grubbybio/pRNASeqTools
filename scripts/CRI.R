options(warn=-1)
args = commandArgs(trailingOnly = T)
pdf('CRI.pdf',6,6)
for(i in 1:length(args)){
  read.table(paste(args[i],"_CRI.txt",sep=""),as.is=T,row.names=1, fill = TRUE, sep = "\t") -> tmp
  plot(density(tmp[,4]),main=args[i],xlab = "CRI",xlim=c(-2,2))
  abline(v=median(tmp[,4]),lty=2,col="grey80")
}
dev.off()
