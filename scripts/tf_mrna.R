options(warn=-1)
args = commandArgs(trailingOnly = T)
norm <- args[1]
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
colnames(b) <- c("genotype","design")
q <- 0
for (n in 1:length(genotype1_tags)){
  for (k in 1:genotype1_reps[n]){
    q = q + 1
    paste(genotype1,genotype1_tags[n],k,sep="_") -> rownames(b)[q]
    genotype1 -> b[q,1]
    genotype1_tags[n] -> b[q,2]
    assign(paste(genotype1,genotype1_tags[n],k,sep="_"), read.table(paste(genotype1,"_",genotype1_tags[n],"_",k,".txt",sep = ""), header = T, as.is = T, row.names = 1))
  }
}
for (n in 1:length(genotype2_tags)){
  for (k in 1:genotype2_reps[n]){
    q = q + 1
    paste(genotype2,genotype2_tags[n],k,sep="_") -> rownames(b)[q]
    genotype2 -> b[q,1]
    genotype2_tags[n] -> b[q,2]
    assign(paste(genotype2,genotype2_tags[n],k,sep="_"), read.table(paste(genotype2,"_",genotype2_tags[n],"_",k,".txt",sep = ""), header = T, as.is = T, row.names = 1))
  }
}
message("Loading completed.")

eval(parse(text=rownames(b)[1]))[,1] -> a
for(i in 2:p){
  cbind(a,eval(parse(text=rownames(b)[i]))[,1])-> a
}
rownames(b) -> colnames(a)
rownames(a) <- rownames(eval(parse(text=rownames(b)[i])))
a[apply(a,1,mean)>=1,] -> a
DESeqDataSetFromMatrix(countData = a, colData = b, design = ~ genotype * design) -> dds
if(norm == "RPM"){
  sizeFactors(dds) <- apply(a,2,sum)/1000000
}
DESeq(dds) -> dds
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
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
pdf(file="tf.pdf",6,6)
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
write.csv(out,paste(genotype2,"vs",genotype1,"tf","csv",sep="."),quote=F)
write.csv(resup,paste(genotype2,"vs",genotype1,"tf","hyper","csv",sep="."),quote=F)
write.csv(resdown,paste(genotype2,"vs",genotype1,"tf","hypo","csv",sep="."),quote=F)
