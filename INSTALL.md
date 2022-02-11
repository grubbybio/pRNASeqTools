# INSTALL
1. Create a new conda environment for `pRNASeqTools` with Python 3.

2. Add the path into the bash profile file (.bash_profile or .bashrc or .zshrc if you are using zsh).
```bash
export PATH=/path/to/pRNASeqTools:$PATH
```
3. You may ask the author for the `reference` files, which should be placed in the `pRNASeqTools` folder.

4. Re-login or `source` the modified bash profile.

5. Run `pRNASeqTools`, see if there are any error messages. Usually it will tell you which dependency is missing.

6. Using `conda` to install neccesary dependencies, until `pRNASeqTools` shows the command menu.

   1. `cutadapt`
      Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. **_EMBnet. journal_**, 2011, 17(1): pp. 10-12.

   2. `ShortStack v3.x`
      Johnson N R, Yeoh J M, Coruh C, _et al_. Improved Placement of Multi-Mapping Small RNAs. **_G3: Genes| Genomes| Genetics_**, 2016: g3. 116.030452.

   3. `samtools v1.x`
      Li H, Handsaker B, Wysoker A, _et al_. The sequence alignment map format and SAMtools. **_Bioinformatics_**, 2009, 25(16): 2078-2079.

   4. [SRA Toolkit](https://github.com/ncbi/sra-tools/)

   5. `R` and `R` packages `DESeq2`, `NMF`, `DMRcaller`, `riboWaltz`, `Seurat`, `RNAmodR.RiboMethSeq`, and `pheatmap`
      R Core Team (2016). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria

      Love M I, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. **_Genome biology_**, 2014, 15(12): 1-21.

      Raivo Kolde (2015). pheatmap: Pretty Heatmaps. R package version 1.0.8.
      **Note: Most neccesary R packages can be installed by `BiocManager`. Please use the branch `devel` of _NMF_ which is compatible with the new internal grid of R.**

   6. `bedtools`
      Quinlan A R, Hall I M. BEDTools: a flexible suite of utilities for comparing genomic features. **_Bioinformatics_**, 2010, 26(6): 841-842.

   7. `bowtie` and `bowtie2`
      Langmead B, Trapnell C, Pop M, _et al._ Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. **_Genome biology_**, 2009, 10(3): R25.
      Langmead B, Salzberg S L. Fast gapped-read alignment with Bowtie 2. **_Nature methods_**, 2012, 9(4): 357-359.

   8. `STAR`
      Dobin A, Davis CA, Schlesinger F, _et al._ STAR: ultrafast universal RNA-seq aligner. **_Bioinformatics_**, 2013, 29(1): 15-21.

   9. `CLIPper`
      Lovci MT, Ghanem D, Marr H, _et al._ Rbfox proteins regulate alternative mRNA splicing through evolutionarily conserved RNA bridges. **_Nature structural & molecular biology_**. 2013;20(12): 1434-1442.

   10. `featureCounts`
      Liao Y, Smyth GK and Shi W. Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014: an efficient general-purpose program for assigning sequence reads to genomic features. ***Bioinformatics***, 30(7):923-30, 2014

   11. [gffread](https://github.com/gpertea/gffread)

   12. `UMItools`
       Smith T S, Heger A and Sudbery I. UMI-tools: Modelling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. ***Genome Res.*** gr.209601.116

7. Specifically, the CLIP-seq analysis needs a customized GFF file located in the `data` folder of `clipper`. The example GFF of _Arabidopsis thaliana_ is available upon request.
