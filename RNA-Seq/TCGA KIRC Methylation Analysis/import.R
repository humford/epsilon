#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#IMPORT

mapper <- NULL

import <- function()
{
  #cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
  cancer.names <- c("KIRC")

  for(cancer in cancer.names)
  {
    setwd(paste("~/Documents/", cancer, sep = ""))
    exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  
    skews <- NULL
    for(i in 1:dim(exprMatrix)[[1]])
    {
      skews[i] <-moment((exprMatrix[i,]), 3)
    }

    skews <- skews[!is.na(skews)]
    names(skews) <- rownames(exprMatrix)   
  
    MvalMatrix <- read.table(paste(cancer, "_Methylation_Processed", sep = ""))
    mapper <- read.table("Methylation_Map")
  }
}