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
 
 clusterTable <- NULL
 
 par(mfrow = c(2,2))
 
 m <- Mclust(geneStats[1:500, c(4,5)], G = 1:4)
 plot(m, what = "classification", xlab = "Change in M-value", main = FALSE)
 title(main = paste("all probes"))
 
 for(f in functions)
 {
     corStats <- cor.test(functionStats[[f]]$Skewness[1:500], functionStats[[f]]$MvalDiff[1:500])  # Take 500 most significant genes
     write.table(t(c(f, corStats$conf.int[1], corStats$estimate, corStats$conf.int[2])), "corrTable", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
     
     m <- Mclust(functionStats[[f]][1:500, c(4,5)], G = 4)
     plot(m, what = "classification", xlab = "Change in M-value", main = FALSE)
     title(main = paste(f, "probes"))
     
     clusterTable <- rbind(clusterTable, c(length(which(m$classification == 1)), length(which(m$classification == 2)), length(which(m$classification == 3)), length(which(m$classification == 4))))
 }
 
 fTest <- fisher.test(clusterTable[c(1,2),])
 write.table(clusterTable, "clusterTable", col.names = FALSE)

 corSampling <- NULL
 
 for(i in 100:min(1000, length(which(functionStats[[1]]$p.value < cutoff)), length(which(functionStats[[2]]$p.value < cutoff)), length(which(functionStats[[3]]$p.value < cutoff))))
 {
    row = i
    for(f in functions)
    {
      corStats <- cor.test(functionStats[[f]]$Skewness[1:i], functionStats[[f]]$MvalDiff[1:i])
      row <- c(row, corStats$estimate, corStats$conf.int[1], corStats$conf.int[2])
    }
    corSampling <- rbind(corSampling, row)
 }
 colnames(corSampling) <- c("num", "promoter_estimate", "promoter_lower", "promoter_upper", "UTR_estimate", "UTR_lower", "UTR_upper","body_estimate", "body_lower", "body_upper" )
 
 for(i in c(3,4,6,7,9,10))
 {
  corSampling[, i] <- smooth.spline(corSampling[, i], spar = 0.75)$y
 }
 
 chart<-
      ggplot(as.data.frame(corSampling)) +
      geom_smooth(aes(x=num, y=promoter_estimate), size=1.0,colour="red", span =  0.7) + 
      scale_x_continuous('Number of Genes',limits=c(100,1000)) +   
      scale_y_continuous('Correlation') +
      geom_ribbon(aes(x=num, ymin = promoter_lower, 
      ymax=promoter_upper), colour="red", fill="red",alpha=0.1) + 
      
      geom_smooth(aes(x=num, y=UTR_estimate), size=1.0,colour="blue", span =  0.7) + 
      geom_ribbon(aes(x=num, ymin = UTR_lower, 
      ymax=UTR_upper), colour="blue", fill="blue",alpha=0.1)+
      
      geom_smooth(aes(x=num, y=body_estimate), size=1.0,colour="green", span =  0.7) + 
      geom_ribbon(aes(x=num, ymin = body_lower, 
      ymax=body_upper), colour="green", fill="green",alpha=0.1) +
      labs(title = paste("Skewness, M-value Correlation vs. Number of Significant Included Genes in", cancer))
 print(chart) 
}
