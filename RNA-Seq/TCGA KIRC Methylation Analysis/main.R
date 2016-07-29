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
  
  setwd("Gene_Methylation")
  genes <- dir()
  
  sigGenes <- NULL
  
  for (symbol in genes)
  {
    GeneMvalMatrix <- read.table(symbol)
    geneProbes <- rownames(GeneMvalMatrix)
    
    TMatrix <- NULL
    NTMatrix <- NULL
    
    splits <- splitter(exprMatrix[symbol,])
    names(splits) <- colnames(exprMatrix)
    
    if(sum(splits) > 0)
    {
      for (patient in intersect(colnames(exprMatrix), colnames(GeneMvalMatrix))) 
      {
        mvalues <- GeneMvalMatrix[,patient]
        
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
      
      if(length(sigprobes) > 0)
      {
        setwd(paste("~/Documents/", cancer, sep = ""))
        write.table(t(c(symbol, length(sigprobes)/length(pvalues), sigprobes)), paste(cancer, "Signifcant_Methylated_Probes", sep = "_"), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
        setwd("Gene_Methylation")
      }
      
      if(length(sigprobes) > 0)
      {
        sigGenes <- rbind(sigGenes, c(symbol, min(adjpvalues), geneProbes[which(adjpvalues == min(adjpvalues))][1], length(sigprobes)))
      }
    }
 }
 colnames(sigGenes) <- c("Symbol", "p-value", "Most_Significant_Probe", "Number_of_Significant_Probes")
 sigGenes <- sigGenes[order(as.numeric(sigGenes[, "p-value"])), ]
 setwd(paste("~/Documents/", cancer, sep = ""))
 write.table(sigGenes, paste(cancer, "Signifcant_Methylated_Genes", sep = "_"), row.names = FALSE, quote = FALSE)
}
