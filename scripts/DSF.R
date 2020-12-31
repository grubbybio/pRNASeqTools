options(warn=-1)
args = commandArgs(trailingOnly = T)
method <- args[1]
pvalueoo <- as.numeric(args[2])
foldchangeoo <- as.numeric(args[3])
args[-c(1:3)] -> args
args[c(TRUE,FALSE)] -> genotype
args[c(FALSE,TRUE)] -> replicates
message("Loading...")
p <- sum(as.numeric(replicates))
data.frame(rep("x",p), stringsAsFactors = FALSE) -> b
colnames(b) <- "condition"
myList <- list()
p <- 0
for (n in 1:length(genotype)){
  for (k in 1:replicates[n]){
    p = p + 1
    paste(genotype[n],"_",k,sep="") -> rownames(b)[p]
    genotype[n] -> b[p,]
    myList[[rownames(b)[p]]] <- read.table(paste(genotype[n],"_",k,".feature.count",sep = ""), as.is = T, row.names = 1)[,1:8]
  }
}
b -> ss
for(n in 1:p){
  read.table(paste(rownames(b)[n],"nf",sep="."),as.is=T) -> nf
  nf[nf[,1]==method,2] -> ss[n,]
}
as.numeric(ss[,1])/1000000 -> ss
suppressMessages(library(DESeq2))
suppressMessages(library(NMF))
nmf.options(grid.patch=TRUE)
library(pheatmap)
library(RColorBrewer)
message("Loading completed.")
message("Processing all length...")
features <- unique(unlist(lapply(myList, rownames)))
myListAll <- list()
a <- matrix(0, nrow = length(features), ncol = length(myList), dimnames = list(features, names(myList)))
for(q in 1:p){
  myListAll[[names(myList)[q]]] <- as.matrix(apply(myList[[q]],1,sum))
}
for (z in seq_along(myListAll)) {
  a[rownames(myListAll[[z]]),z] <- myListAll[[z]]
}
a[apply(a, 1, mean) >= 1,] -> a
try({
  suppressMessages(DESeqDataSetFromMatrix(countData = a, colData = b, design = ~ condition) -> dds)
  sizeFactors(dds) <- ss
  suppressMessages(DESeq(dds) -> dds)
  select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:1000]
  nt <- normTransform(dds)
  log2.norm.counts <- assay(nt)[select,]
  df <- as.data.frame(colData(dds)[,c("condition")])
  rownames(b) -> rownames(df)
  colnames(df) <- "Genotype"
  rld <- rlog(dds)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- rownames(b)
  colors <- colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
  pdf(file=paste(method,"gene","all","pdf",sep="."),6,6)
  aheatmap(log2.norm.counts, Rowv=NA, Colv=NA, annCol=df)
  message("Top 1000 completed!")
  print(plotPCA(rld,ntop=1000,intgroup=c("condition")))
  message("PCA completed!")
  pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,show_colnames = F,col=colors)
  message("Dist completed!")
  dev.off()
  for(j in 2:length(genotype)){
    results(dds,contrast = c("condition",genotype[j],genotype[1])) -> res
    cbind(as.data.frame(res),counts(dds,norm = T)) -> out
    subset(out, pvalue < pvalueoo & log2FoldChange >= log2(foldchangeoo)) -> resup
    subset(out, pvalue < pvalueoo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
    write.csv(out,paste(genotype[j],"vs",genotype[1],method,"feature","all","csv",sep="."),quote=F)
    write.csv(resup,paste(genotype[j],"vs",genotype[1],method,"feature","all","hyper","csv",sep="."),quote=F)
    write.csv(resdown,paste(genotype[j],"vs",genotype[1],method,"feature","all","hypo","csv",sep="."),quote=F)
  }
}, silent = T)

for(i in 2:8){
  message(paste("Processing ",i+17,"nt...",sep=""))
  a <- matrix(0, nrow = length(features), ncol = length(myList), dimnames = list(features, names(myList)))
  for(q in 1:p){
    myListAll[[names(myList)[q]]] <- myList[[q]][i]
  }
  for (z in seq_along(myListAll)) {
    a[rownames(myListAll[[z]]), z] <- as.matrix(myListAll[[z]])
  }
  a[apply(a, 1, mean) >= 1,] -> a
  try({
    suppressMessages(DESeqDataSetFromMatrix(countData = a, colData = b, design = ~ condition) -> dds)
    sizeFactors(dds) <- ss
    suppressMessages(DESeq(dds) -> dds)
    select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:1000]
    nt <- normTransform(dds)
    log2.norm.counts <- assay(nt)[select,]
    df <- as.data.frame(colData(dds)[,c("condition")])
    rownames(b) -> rownames(df)
    colnames(df) <- "Genotype"
    rld <- rlog(dds)
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- rownames(b)
    colors <- colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
    pdf(file=paste(method,"gene",i+17,"pdf",sep="."),6,6)
    aheatmap(log2.norm.counts, Rowv=NA, Colv=NA, annCol=df)
    message("Top 1000 completed!")
    print(plotPCA(rld,ntop=1000,intgroup=c("condition")))
    message("PCA completed!")
    pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,show_colnames = F,col=colors)
    message("Dist completed!")
    dev.off()
    for(j in 2:length(genotype)){
      results(dds,contrast = c("condition",genotype[j],genotype[1])) -> res
      cbind(as.data.frame(res),counts(dds,norm = T)) -> out
      subset(out, pvalue < pvalueoo & log2FoldChange >= log2(foldchangeoo)) -> resup
      subset(out, pvalue < pvalueoo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
      write.csv(out,paste(genotype[j],"vs",genotype[1],method,"feature",i+17,"csv",sep="."),quote=F)
      write.csv(resup,paste(genotype[j],"vs",genotype[1],method,"feature",i+17,"hyper","csv",sep="."),quote=F)
      write.csv(resdown,paste(genotype[j],"vs",genotype[1],method,"feature",i+17,"hypo","csv",sep="."),quote=F)
    }
  }, silent = T)
}
