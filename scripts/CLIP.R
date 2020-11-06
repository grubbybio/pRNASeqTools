options(warn=-1)
args = commandArgs(trailingOnly = T)
pvalueoo <- as.numeric(args[1])
foldchangeoo <- as.numeric(args[2])
args[-c(1:2)] -> args
args[1:(length(args)/2)] -> genotype
args[1:(length(args)/2)] -> replicates
for(m in 1:length(args)){
  if(m %% 2 == 0){
    args[m] -> replicates[m/2]
  }else{
    args[m] -> genotype[(m+1)/2]
  }
}
message("Samples OK! loading...")
suppressMessages(library(DESeq2))
suppressMessages(library(NMF))
library(RColorBrewer)
library(pheatmap)
nmf.options(grid.patch=TRUE)
p <- sum(as.numeric(replicates))
data.frame(rep("x",p), stringsAsFactors = FALSE) -> b
colnames(b) <- "genotype"
p <- 0
for (n in 1:length(genotype)){
  for (k in 1:replicates[n]){
    p = p + 1
    paste(genotype[n],"_",k,sep="") -> rownames(b)[p]
    assign(rownames(b)[p], read.table(paste(rownames(b)[p],".txt",sep = ""), as.is = T, row.names = 1))
    genotype[n] -> b[p,1]
  }
}
b -> ss
for(n in 1:p){
  read.table(paste(rownames(b)[n],"nf",sep="."),as.is=T)[1,2] -> ss[n,1]
}
as.numeric(ss[,1])/1000000 -> ss
eval(parse(text=rownames(b)[1]))[,1] -> a
for(q in 2:p){
  cbind(a,eval(parse(text=rownames(b)[q]))[,1]) -> a
}
rownames(a) <- rownames(eval(parse(text=rownames(b)[1])))
colnames(a) <- rownames(b)
a[apply(a,1,mean)>=2,] -> a
message("Loading completed.")
DESeqDataSetFromMatrix(countData = a, colData = b, design = ~ genotype) -> dds
sizeFactors(dds) <- ss
DESeq(dds) -> dds
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:1000]
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("genotype")])
rownames(b)[1:p] -> rownames(df)
colnames(df) <- "Genotype"
rld <- rlog(dds)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pdf(file="total.pdf",6,6)
aheatmap(log2.norm.counts, Rowv=NA, Colv=NA, annCol=df)
message("Top 1000 completed!")
plotPCA(rld,ntop=1000,intgroup=c("genotype"))
message("PCA completed!")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,show_colnames = F,col=colors)
message("Dist completed!")
dev.off()
for(j in 2:length(genotype)){
  results(dds,contrast = c("genotype",genotype[j],genotype[1])) -> res
  cbind(as.data.frame(res),counts(dds,norm = T)) -> out
  subset(out, pvalue < pvalueoo & log2FoldChange >= log2(foldchangeoo)) -> resup
  subset(out, pvalue < pvalueoo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
  write.csv(out,paste(genotype[j],"vs",genotype[1],"csv",sep="."),quote=F)
  write.csv(resup,paste(genotype[j],"vs",genotype[1],"hyper","csv",sep="."),quote=F)
  write.csv(resdown,paste(genotype[j],"vs",genotype[1],"hypo","csv",sep="."),quote=F)
}
