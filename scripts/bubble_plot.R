options(warn=-1)
args = commandArgs(trailingOnly = T)
if(length(args) %% 2 != 0){
  stop("The arguments are not correct!", call. = FALSE)
}else{
  args[1:(length(args)/2)] -> genotype
  args[1:(length(args)/2)] -> replicates
  for(m in 1:length(args)){
    if(m %% 2 == 0){
      args[m] -> replicates[m/2]
    }else{
      args[m] -> genotype[(m+1)/2]
    }
  }
}
message("Samples OK! loading...")
p <- sum(as.numeric(replicates))
rep("x",p) -> a
q <- 0
for (n in 1:length(genotype)){
  for (k in 1:replicates[n]){
    q = q + 1
    paste(genotype[n],"_",k,".out",sep="") -> a[q]
  }
}
for (i in 1:length(a)){
  assign(a[i],read.table(a[i],as.is = T))
}
unique(eval(parse(text=a[1]))[,1]) -> llist
for (m in 1:length(llist)){
  pdf(paste(llist[m],"pdf",sep = "."), length(genotype)*8, max(as.numeric(replicates))*8)
  par(mfcol=c(max(as.numeric(replicates)),length(genotype)))
  for (i in 1:length(a)){
    eval(parse(text=a[i])) -> tmp
    tmp[tmp[,1]==llist[m],] ->tmp2
    max(tmp2[,4])+1 -> bb
    sum(tmp2[,4]) -> ss
    for(k in 1:nrow(tmp2)){
      sqrt(tmp2[k,4] / bb) -> tmp2[k,4]
    }
    cc <- rainbow(nrow(tmp2),start = 0.05, end = 0.95)
    try({
      symbols(x=tmp2[,2] * 3, y=tmp2[,3] * 3, circles = tmp2[,4], ylim=c(-2,26), xlim=c(-2,26),xlab = "Tailing", ylab="Truncation", inches = F, bg = cc, xaxt="n",yaxt="n", main=paste(a[i],ss,sep="\n"))
    }, silent = T)
    axis(1,at=(0:10)*3,labels = 0:10)
    axis(2,at=(0:10)*3,labels = 0:-10)
  }
  dev.off()
}
