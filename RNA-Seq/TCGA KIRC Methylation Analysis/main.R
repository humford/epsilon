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
    
    for (symbol in colnames(splitMatrix)) 
    {
      splitMatrix[, symbol] <- split(exprMatrix[, symbol])
    }
    
    for (gene in XXXX)
    {
      for (patient in gene) 
      {
        TMatrix <- NULL
        NTMatrix <- NULL
        mvalues <- MValues(bvalues[probes(gene), patient])
        
        if (splitMatrix[patient, symbol]  == 1) 
        {
          TMatrix <- cbind(TMatrix, mvalues)
        } 
        else if (splitMatrix[patient, symbol] == 0)
        {
          NTMatrix <- cbind(NTMatrix, mvalues)
        }
      }
      
      for (probe in rowNames(TMatrix))
      {
        pvalues[probe] <- wilcox.test(TMatrix[probe], NTMatrix[probe], correct = FALSE)
        adjpvalues <- p.adjust(pvalue[probe], method = "BH")
        sigprobes <- names(pvalues[which(adjpvalues < cutoff)])
        
        write.table(c(gene, length/length(pvalues), sigprobes), append = TRUE)
      }
    }
  }
}