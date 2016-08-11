#TCGA KIRC METHYLATION ANALYSIS
#BEN CHURCH AND HENRY WILLIAMS
#MAIN

source("~/Documents/Git/epsilon/RNA-Seq/TCGA KIRC Methylation Analysis/globals.R")


#cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
cancer.names <- c("SKCM")

functions <- c("promoter", "UTR", "body")

for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  
  setwd("Gene_Methylation")
  genes <- dir()
  
  geneStats <- data.frame(p.value = numeric(0), NumSigProbes = numeric(0), PercentSigProbes = numeric(0), MvalDiff = numeric(0), Skewness = numeric(0))
  
  functionStats <- list(data.frame(p.value = numeric(0), NumSigProbes = numeric(0), PercentSigProbes = numeric(0), MvalDiff = numeric(0), Skewness = numeric(0)),
    data.frame(p.value = numeric(0), NumSigProbes = numeric(0), PercentSigProbes = numeric(0), MvalDiff = numeric(0), Skewness = numeric(0)),
    data.frame(p.value = numeric(0), NumSigProbes = numeric(0), PercentSigProbes = numeric(0), MvalDiff = numeric(0), Skewness = numeric(0)))
  
  names(functionStats) <- functions

  
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
        write.table(t(c(symbol, length(sigProbes)/length(pvalues), sigProbes)), paste(cancer, "Significant_Methylated_Probes", sep = "_"), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
        setwd(paste("~/Documents/", cancer, "/Gene_Methylation", sep = ""))
      }
      
      geneStats <- addRow(geneStats, c(min(adjpvalues), length(sigProbes), length(sigProbes)/length(geneProbes), 
                                          mean(rowMeans(TMatrix)[sigProbes] - rowMeans(NTMatrix)[sigProbes]), moment(exprMatrix[symbol, ], 3)))
      rownames(geneStats)[dim(geneStats)[1]] <- symbol
      
      probe_functions <- unlist(lapply(geneProbes, function(x) probe_function(x, symbol)))
      
      for(f in functions)
      {
        function.probes <- geneProbes[which(probe_functions == f)]
        if(length(function.probes) > 0)
        {
          adjpvalues <- p.adjust(pvalues[function.probes], method = "BH")
          sigProbes <- function.probes[which(adjpvalues < cutoff)]
        
          functionStats[[f]] <- addRow(functionStats[[f]], c(min(adjpvalues), length(sigProbes), length(sigProbes)/length(function.probes), 
                                          mean(rowMeans(TMatrix)[sigProbes] - rowMeans(NTMatrix)[sigProbes]), moment(exprMatrix[symbol, ], 3)))
          rownames(functionStats[[f]])[dim(functionStats[[f]])[1]] <- symbol
        }
      }
    }
 }
 
 geneStats$p.value <- p.adjust(geneStats$p.value, method = "BH")
 geneStats <- geneStats[order(geneStats$p.value), ]
 significantGenes <- geneStats[which(geneStats$p.value < cutoff), ]
 setwd(paste("~/Documents/", cancer, "/Results", sep = ""))
 write.table(cbind(Symbol = rownames(significantGenes), format(significantGenes, digits = 3)) , paste(cancer, "Significant_Methylated_Genes", sep = "_"), row.names = FALSE, quote = FALSE,)

 for(f in functions)
 {
    functionStats[[f]]$p.value <- p.adjust(functionStats[[f]]$p.value, method = "BH")
    functionStats[[f]] <- functionStats[[f]][order(functionStats[[f]]$p.value), ]
 }
 
}         

