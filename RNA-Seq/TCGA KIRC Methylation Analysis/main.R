#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#MAIN

source("~/Documents/Git/epsilon/RNA-Seq/TCGA KIRC Methylation Analysis/globals.R")


#cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
cancer.names <- c("SKCM", "LGG")

for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  
  setwd("Gene_Methylation")
  genes <- dir()
  
  geneStats <- data.frame(p.value = numeric(0), NumSigProbes = numeric(0), PercentSigProbes = numeric(0), MvalDiff = numeric(0), Skewness = numeric(0))
  
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
      sigProbes <- geneProbes[which(adjpvalues < cutoff)]
      
      if(length(sigProbes) > 0)
      {
        setwd(paste("~/Documents/", cancer, "/Results", sep = ""))
        write.table(t(c(symbol, length(sigProbes)/length(pvalues), sigProbes)), paste(cancer, "Signifcant_Methylated_Probes", sep = "_"), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
        setwd("Gene_Methylation")
      }
      
      geneStats <- addRow(geneStats, c(min(adjpvalues), length(sigProbes), length(sigProbes)/length(geneProbes), 
                                          rowMeans(TMatrix)[order(adjpvalues)[1]] - rowMeans(NTMatrix)[order(adjpvalues)[1]], moment(exprMatrix[symbol, ], 3)))
      rownames(geneStats)[dim(geneStats)[1]] <- symbol
    }
 }
 
 geneStats$p.value <- p.adjust(geneStats$p.value, method = "BH")
 geneStats <- geneStats[order(geneStats$p.value), ]
 significantGenes <- geneStats[which(geneStats$p.value < cutoff), ]
 setwd(paste("~/Documents/", cancer, "/Results", sep = ""))
 write.table(cbind(Symbol = rownames(significantGenes), format(significantGenes, digits = 3)) , paste(cancer, "Signifcant_Methylated_Genes", sep = "_"), row.names = FALSE, quote = FALSE,)
 
 plot(-log(significantGenes$p.value), abs(significantGenes$Skewness))
 plot(significantGenes$MvalDiff[1:100], significantGenes$Skewness[1:100])
}
