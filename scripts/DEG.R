options(warn=-1)
args = commandArgs(trailingOnly = T)
norm = args[1]
pvalueoo <- as.numeric(args[2])
fdroo <- as.numeric(args[3])
foldchangeoo <- as.numeric(args[4])
try(
  {read.table(paste(args[5],"/reference/",args[6],".BIN",sep=""), as.is = T, sep="\t") -> binref
  as.matrix(table(binref[,2])) -> refbin},
  silent=TRUE
)
args[-c(1:6)] -> args
args[c(TRUE, FALSE)] -> genotype
args[c(FALSE, TRUE)] -> replicates
message("Parameters OK! loading...")
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
    genotype[n] -> b[p,]
    assign(rownames(b)[p], read.table(paste(genotype[n],"_",k,".txt",sep = ""), as.is = T, row.names = 1,header = T))
  }
}
message("Loading completed.")
eval(parse(text=rownames(b)[1]))[,1] -> a
for(q in 2:p){
  cbind(a,eval(parse(text=rownames(b)[q]))[,1]) -> a
}
rownames(a) <- rownames(eval(parse(text=rownames(b)[1])))
colnames(a) <- rownames(b)
apply(a,2,sum)/1000000 -> rtotal
a[apply(a,1,mean)>=1,] -> a
DESeqDataSetFromMatrix(countData = a, colData = b, design = ~ genotype) -> dds
if(norm == "RPM"){
  sizeFactors(dds) <- rtotal
}
DESeq(dds) -> dds
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
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
  cbind(as.data.frame(res),t(t(counts(dds))/rtotal)) -> out
  if(fdroo < 1){
    subset(out, padj < fdroo & log2FoldChange >= log2(foldchangeoo)) -> resup
    subset(out, padj < fdroo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
  }else{
    subset(out, pvalue < pvalueoo & log2FoldChange >= log2(foldchangeoo)) -> resup
    subset(out, pvalue < pvalueoo & log2FoldChange <= -log2(foldchangeoo)) -> resdown
  }
  write.csv(out,paste(genotype[j],"vs",genotype[1],"total","csv",sep="."),quote=F)
  write.csv(resup,paste(genotype[j],"vs",genotype[1],"total","upregulated","csv",sep="."),quote=F)
  write.csv(resdown,paste(genotype[j],"vs",genotype[1],"total","downregulated","csv",sep="."),quote=F)
  try(
    {binref[binref[,1] %in% rownames(res),2] -> tmp1
      as.matrix(table(tmp1)) -> tmp1
      binref[binref[,1] %in% rownames(resup),2] -> tmp2
      as.matrix(table(tmp2)) -> tmp2
      binref[binref[,1] %in% rownames(resdown),2] -> tmp3
      as.matrix(table(tmp3)) -> tmp3
      merge(x=refbin, y=tmp1, all = T, by = 0) -> tmp
      rownames(tmp) <- tmp[,1]
      tmp[,-1] -> tmp
      merge(x=tmp, y=tmp2, all = T, by = 0) -> tmp
      rownames(tmp) <- tmp[,1]
      tmp[,-1] -> tmp
      merge(x=tmp, y=tmp3, all = T, by = 0) -> tmp
      rownames(tmp) <- tmp[,1]
      tmp[,-1] -> tmp
      sum(tmp[!is.na(tmp[,2]),2]) -> sumref
      sum(tmp[!is.na(tmp[,3]),3]) -> sumup
      sum(tmp[!is.na(tmp[,4]),4]) -> sumdown
      tmp[,3:4] -> tmp[,5:6]
      for(k in 1:nrow(tmp)){
        if(!is.na(tmp[k,3])){
          matrix(c(tmp[k,2],sumref-tmp[k,2],tmp[k,3],sumup-tmp[k,3]),nrow=2) -> fis
          fisher.test(fis)$p.value -> tmp[k,5]
        }
        if(!is.na(tmp[k,4])){
          matrix(c(tmp[k,2],sumref-tmp[k,2],tmp[k,4],sumdown-tmp[k,4]),nrow=2) -> fis
          fisher.test(fis)$p.value -> tmp[k,6]
        }
      }
      p.adjust(tmp[,5],"fdr") -> tmp[,7]
      p.adjust(tmp[,6],"fdr") -> tmp[,8]
      colnames(tmp) <- c("GENOME","DETECTED","UP","DOWN","UP P","DOWN P","UP FDR","DOWN FDR")
      write.table(tmp,paste(genotype[j],"vs",genotype[1],"total","bin","txt",sep="."),quote=F,sep="\t")},
    silent = TRUE
  )
}
