suppressMessages(library(riboWaltz))
options(warn=-1)
args = commandArgs(trailingOnly = T)
genome <- args[1]
tags <- args[-1]
ann <- create_annotation(gtfpath = paste(genome,'gtf',sep="."))
reads_list <- bamtolist(bamfolder = '.', annotation = ann)
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")
reads_psite_list <- psite_info(reads_list, psite_offset)
#coverage_dt <- codon_coverage(reads_psite_list, ann, psite = FALSE)
#psite_cds <- psite_per_cds(reads_psite_list, ann)

for(i in 1:length(tags)){
  write.csv(eval(parse(text=paste("reads_psite_list$",tags[i]))),paste(tags[i],"csv",sep="."),quote = F,row.names = F)
  example_length_dist_zoom <- rlength_distr(reads_list, sample = tags[i], cl = 99)
  example_ends_heatmap <- rends_heat(reads_list, ann, sample = tags[i], cl = 85, utr5l = 25, cdsl = 40, utr3l = 25)
  example_psite_region <- region_psite(reads_psite_list, ann, sample = tags[i])
  example_frames_stratified <- frame_psite_length(reads_psite_list, sample = tags[i], region = "all", cl = 90)
  example_frames <- frame_psite(reads_psite_list, sample = tags[i], region = "all")

  pdf(file=paste(tags[i],"pdf",sep="."),12,6)
  print(example_length_dist_zoom[[paste("plot",tags[i],sep = "_")]])
  print(example_ends_heatmap[["plot"]])
  print(example_psite_region[["plot"]])
  print(example_frames_stratified[["plot"]])
  print(example_frames[["plot"]])
  for (j in c(19:21,28:32)){
    try({
      example_metaprofile <- metaprofile_psite(reads_psite_list, ann, sample = tags[i], utr5l = 20, cdsl = 40, utr3l = 20, plot_title = "sample.transcript.length_range", length_range = j)
      print(example_metaprofile[[paste("plot",tags[i],sep = "_")]])
    }, silent = TRUE)
  }
  #codon_usage_barplot <- codon_usage_psite(reads_psite_list, ann, sample = tags[i], fastapath = paste(genome,'fa',sep="."), fasta_genome = F)
  dev.off()
}
