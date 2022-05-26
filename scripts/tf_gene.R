options(warn=-1)
args = commandArgs(trailingOnly = T)
method <- args[1]
pvalueoo <- as.numeric(args[2])
foldchangeoo <- as.numeric(args[3])
args[-c(1:3)] -> args
m <- length(args)/2
args[1] -> genotype1
args[2:m] -> genotype1_tags
args[m+1] -> genotype2
args[(m+2):(2*m)] -> genotype2_tags
genotype1_tags[c(FALSE, TRUE)] -> genotype1_reps
genotype1_tags[c(TRUE, FALSE)] -> genotype1_tags
genotype2_tags[c(FALSE, TRUE)] -> genotype2_reps
genotype2_tags[c(TRUE, FALSE)] -> genotype2_tags
suppressMessages(library(DESeq2))
suppressMessages(library(NMF))
library(RColorBrewer)
library(pheatmap)
nmf.options(grid.patch=TRUE)
message("Parameters OK! loading...")

p <- sum(as.numeric(genotype1_reps)+as.numeric(genotype2_reps))
data.frame(rep("x",p), stringsAsFactors = FALSE) -> b
cbind(b,b) -> b
b -> ss
colnames(b) <- c("genotype","design")
q <- 0
for (n in 1:length(genotype1_tags)){
  for (k in 1:genotype1_reps[n]){
    q = q + 1
    paste(genotype1,genotype1_tags[n],k,sep="_") -> rownames(b)[q]
    genotype1 -> b[q,1]
    genotype1_tags[n] -> b[q,2]
    assign(paste(genotype1,genotype1_tags[n],k,sep="_"), read.table(paste(genotype1,"_",genotype1_tags[n],"_",k,".gene.count",sep = ""), header = T, as.is = T, row.names = 1))
  }
}
for (n in 1:length(genotype2_tags)){
  for (k in 1:genotype2_reps[n]){
    q = q + 1
    paste(genotype2,genotype2_tags[n],k,sep="_") -> rownames(b)[q]
    genotype2 -> b[q,1]
    genotype2_tags[n] -> b[q,2]
    assign(paste(genotype2,genotype2_tags[n],k,sep="_"), read.table(paste(genotype2,"_",genotype2_tags[n],"_",k,".gene.count",sep = ""), header = T, as.is = T, row.names = 1))
  }
}
for(n in 1:p){
  read.table(paste(rownames(b)[n],"nf",sep="."),as.is=T) -> nf
  nf[nf[,1]==method,2] -> ss[n,]
}
as.numeric(ss[,1])/1000000 -> ss
message("Loading completed.")

for(i in 2:8){
  message(paste("Processing ",i+17,"nt...",sep=""))
  eval(parse(text=rownames(b)[1]))[,i] -> a
  for(q in 2:p){
    cbind(a,eval(parse(text=rownames(b)[q]))[,i]) -> a
  }
  rownames(a) <- rownames(eval(parse(text=rownames(b)[1])))
  colnames(a) <- rownames(b)
  a[apply(a, 1, mean) >= 1,] -> a
  try({
    suppressMessages(DESeqDataSetFromMatrix(countData = a, colData = b, design = ~ genotype * design) -> dds)
    sizeFactors(dds) <- ss
    suppressMessages(DESeq(dds) -> dds)
    select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:1000]
    nt <- normTransform(dds)
    log2.norm.counts <- assay(nt)[select,]
    df <- as.data.frame(colData(dds)[,c("genotype","design")])
    rownames(b)[1:p] -> rownames(df)
    colnames(df) <- c("Genotype","Design")
    rld <- rlog(dds)
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- rownames(b)
    colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
    pdf(file=paste(method,"gene",i+17,"tf","pdf",sep="."),6,6)
    aheatmap(log2.norm.counts, Rowv=NA, Colv=NA, annCol=df)
    message("Top 1000 completed!")
    plotPCA(rld,ntop=1000,intgroup=c("genotype","design"))
    message("PCA completed!")
    pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,show_colnames = F,col=colors)
    message("Dist completed!")
    dev.off()
    results(dds,name = resultsNames(dds)[4]) -> res
    cbind(as.data.frame(res),counts(dds,norm = T)) -> out
    subset(out, pvalue < pvalueoo & log2FoldChange >= log2(foldchangeoo)) -> resup
    subset(out, pvalue < pvalueoo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
    write.csv(out,paste(genotype2,"vs",genotype1,method,"gene",i+17,"tf","csv",sep="."),quote=F)
    write.csv(resup,paste(genotype2,"vs",genotype1,method,"gene",i+17,"tf","hyper","csv",sep="."),quote=F)
    write.csv(resdown,paste(genotype2,"vs",genotype1,method,"gene",i+17,"tf","hypo","csv",sep="."),quote=F)
  }, silent = T)
}
