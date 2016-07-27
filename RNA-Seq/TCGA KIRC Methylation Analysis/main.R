#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#MAIN

source("~/Documents/Git/epsilon/RNA-Seq/TCGA KIRC Methylation Analysis/globals.R")


#cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
cancer.names <- c("KIRC")

for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  
  
  for (symbol in rownames(exprMatrix))
  {
    setwd("Gene_Methylation")
    
    MvalMatrix <- read.table(symbol)
    geneProbes <- rownames(MvalMatrix)
    
    TMatrix <- NULL
    NTMatrix <- NULL
    
    splits <- splitter(exprMatrix[symbol,])
    names(splits) <- colnames(exprMatrix)
    
    for (patient in intersect(colnames(exprMatrix), colnames(MvalMatrix))) 
    {
      mvalues <- MvalMatrix[,patient]
      
      if (splits[patient]) 
      {
        TMatrix <- cbind(TMatrix, mvalues)
      } 
      else
      {
        NTMatrix <- cbind(NTMatrix, mvalues)
      }
    }
    
    rownames(TMatrix) <- geneProbes
    rownames(NTMatrix) <- geneProbes
    
    pvalues <- NULL
    
    for (probe in geneProbes)
    {
      pvalues[probe] <- wilcox.test(TMatrix[probe,], NTMatrix[probe,], correct = FALSE)$p.value
    }
    
    adjpvalues <- p.adjust(pvalues, method = "BH")
    sigprobes <- geneProbes[which(adjpvalues < cutoff)]
    
    setwd(paste("~/Documents/", cancer, sep = ""))
    write.table(t(c(symbol, length(sigprobes)/length(pvalues), sigprobes)), paste(cancer, "Signifcant_Methylated_Genes", sep = "_"), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
