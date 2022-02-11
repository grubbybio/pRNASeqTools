message("\nChecking R packages...")

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}

packages <- c("DESeq2", "pheatmap", "DMRcaller", "RNAmodR.RiboMethSeq")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if(!suppressMessages(require(x, character.only=T, quietly = TRUE))){
      message(paste("Installing ",x,"...",sep=""))
      BiocManager::install(x)
    }
  }
)

if(!suppressMessages(require("riboWaltz", character.only=T, quietly = TRUE))){
  message("Installing riboWaltz...")
  devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
}

if(!suppressMessages(require("NMF", character.only=T, quietly = TRUE))){
  message("Installing NMF...")
  devtools::install_github("renozao/NMF@devel", dependencies = TRUE)
}

if(!suppressMessages(require("Seurat", character.only=T, quietly = TRUE))){
  message("Installing Seurat...")
  remotes::install_github("satijalab/seurat")
  remotes::install_github("jlmelville/uwot")
}

message("All R packages are installed!")
