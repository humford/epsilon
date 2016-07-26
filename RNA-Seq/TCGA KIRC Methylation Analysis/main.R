#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#MAIN

source("globals.R")
source("import.R")

main <- function() 
{
  #cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
  cancer.names <- c("KIRC")
  
  for(cancer in cancer.names)
  {
    setwd(paste("~/Documents/", cancer, sep = ""))
    exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
    
    MvalMatrix <- read.table(paste(cancer, "_Methylation_Processed", sep = ""))
    mapper <- read.table("Methylation_Map")
    
    splitMatrix <- exprMatrix
    
    for (symbol in rownames(splitMatrix)) 
    {
      splitMatrix[, symbol] <- split(exprMatrix[, symbol])
    }
    
    
    for (symbol in rownames(splitMatrix))
    {
      geneProbes <- probes(symbol)
      
      TMatrix <- NULL
      NTMatrix <- NULL
      
      for (patient in colnames(splitMatrix)) 
      {
        mvalues <- MvalMatrix[geneProbes, patient]
        
        if (splitMatrix[patient, symbol]  == 1) 
        {
          TMatrix <- cbind(TMatrix, mvalues)
        } 
        else if (splitMatrix[patient, symbol] == 0)
        {
          NTMatrix <- cbind(NTMatrix, mvalues)
        }
      }
      
      rownames(TMatrix) <- geneProbes
      rownames(NTMatrix) <- geneProbes
      
      pvalues <- NULL
      
      for (probe in geneProbes)
      {
        pvalues[probe] <- wilcox.test(TMatrix[probe], NTMatrix[probe], correct = FALSE)
      }
      
      names(pvales) <- geneProbes
      
      adjpvalues <- p.adjust(pvalues, method = "BH")
      sigprobes <- geneProbes[which(adjpvalues < cutoff)]
      write.table(c(gene, length(sigprobes)/length(pvalues), sigprobes), append = TRUE)
    }
  }
}