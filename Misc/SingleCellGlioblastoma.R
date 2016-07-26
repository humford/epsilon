#Load Files

setwd("~/Documents/Glioblastoma")

exprMatrix <- read.table("GSE57872.csv", sep = ",", header = TRUE)

rownames(exprMatrix) <- exprMatrix[,1]

exprMatrix <- as.matrix(exprMatrix[,-1])


#Analysis of Moment Distributions

skews <- NULL
for(i in 1:dim(exprMatrix)[[1]])
{
  skews[i] <-moment((exprMatrix[i,]), 3)
}

names(skews) <- rownames(exprMatrix)

plot(density(skews), main = "Glioblastoma", xlab= "Skewness")




  lower <- NULL

  for(a in 1:dim(exprMatrix)[[2]])
  {
    skews <- NULL
    for(i in 1:dim(exprMatrix)[[1]])
    {
      skews[i] <- moment(sample(exprMatrix[i,], a), 3)
    }
    lower[a] <- length(which(skews < 0))
  }

  lower[1:3] <- lower[4]

  lower <- lower/length(skews)

  range <- 0.05 *abs(max(lower) - min(lower))
  accept <- which(abs(lower  - lower[length(lower)])  <= range)
  reject <- which(abs(lower  - lower[length(lower)]) > range)


  plot(reject, lower[reject],col = "red", main= "Single Cell Glioblastoma Gene Spliting by Sample Size", xlab = "Samples", ylab = "Fraction of Genes in Lower Mode", 
                                            xlim = c(1, length(lower)), ylim = c(min(lower) - range, max(lower) + range))
  points(accept, lower[accept], col = "black")

  abline(a = (lower[length(lower)] - range), b = 0)
  abline(a = (lower[length(lower)] + range), b = 0)
  text(min(accept), lower[min(accept)] + 2*range, as.character(min(accept)))


  
  
setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("RNASeqControlProcessed.txt"))
setwd("~/Documents/Graphs")

ControlSkews <- NULL

for(i in 1:dim(controlMatrix)[[1]])
{
  ControlSkews[i] <-moment((controlMatrix[i,]), 3)
}

names(ControlSkews) <- rownames(exprMatrix)

cancer.names <- "Single Cell Glio"

library(GOstats)
library(KEGG.db)
library(org.Hs.eg.db)
library(mclust)

controlDiffTable[["Single Cell Glio"]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

for(cancer in cancer.names)
{
  GeneLists[[cancer]] <- list() 
  
  out <- Mclust(controlDiffTable[[cancer]], G = 3, modelNames = "V")
    
  GeneLists[[cancer]][["UP"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 3)] 
  GeneLists[[cancer]][["MIDDLE"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 2)] 
  GeneLists[[cancer]][["DOWN"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 1)]  
  GeneLists[[cancer]][["ALL"]] <- names(controlDiffTable[[cancer]])
}

GeneLists[["Intersect"]] <- GeneLists[["OV"]]

for(cancer in cancer.names)
{
  GeneLists[["Intersect"]][["UP"]] <- intersect(GeneLists[["Intersect"]][["UP"]], GeneLists[[cancer]][["UP"]])
  GeneLists[["Intersect"]][["MIDDLE"]] <- intersect(GeneLists[["Intersect"]][["MIDDLE"]], GeneLists[[cancer]][["MIDDLE"]]) 
  GeneLists[["Intersect"]][["DOWN"]] <- intersect(GeneLists[["Intersect"]][["DOWN"]], GeneLists[[cancer]][["DOWN"]]) 
  GeneLists[["Intersect"]][["ALL"]] <- intersect(GeneLists[["Intersect"]][["ALL"]], GeneLists[[cancer]][["ALL"]])
}

hgCutoff <- 1
GeneStatsResults <- NULL
AdjustedResults <- NULL

for(cancer in cancer.names)
{
  GeneStatsResults[[cancer]] <- list()
  AdjustedResults[[cancer]] <- list()
  
  universelist <- GeneLists[[cancer]][["ALL"]]
  
  universe <- unlist(mget(intersect(universelist, ls(org.Hs.egSYMBOL2EG)), org.Hs.egSYMBOL2EG))

  universe <- unique(universe[!is.na(universe)])
  
  for(group in c("UP", "MIDDLE", "DOWN"))
  {
    
    GeneStatsResults[[cancer]][[group]] <- list()
    AdjustedResults[[cancer]][[group]] <- list()

    genelist <- GeneLists[[cancer]][[group]]
      
    genelist.eg <- unlist(mget(intersect(genelist, ls(org.Hs.egSYMBOL2EG)), org.Hs.egSYMBOL2EG))
  
    genelist.eg <- unique(genelist.eg[!is.na(genelist.eg)])
    
    for(marker in c("BP", "CC", "MF", "KEGG"))
    {
    
      if(marker == "KEGG")
      {
          params <- new("KEGGHyperGParams", geneIds=genelist.eg, universeGeneIds=universe, 
                annotation="org.Hs.eg.db", pvalueCutoff=hgCutoff, testDirection="over") 
      
      }
      else
      {
          params <- new("GOHyperGParams", geneIds=genelist.eg, universeGeneIds=universe, 
                annotation="org.Hs.eg.db", ontology=marker, pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")               
      }
      hgOver <- hyperGTest(params)
 
      GeneStatsResults[[cancer]][[group]][[marker]] <- summary(hgOver)
      adjustedPvalues <- p.adjust(GeneStatsResults[[cancer]][[group]][[marker]]$Pvalue, method = "BH")
      AdjustedResults[[cancer]][[group]][[marker]] <- GeneStatsResults[[cancer]][[group]][[marker]][which(adjustedPvalues < 0.05),]
      AdjustedResults[[cancer]][[group]][[marker]]$Pvalues <- adjustedPvalues[which(adjustedPvalues < 0.05)]
    }
  }
}


for(marker in c("BP", "CC", "MF", "KEGG"))
{
  for(group in c("UP", "MIDDLE", "DOWN"))
  {
    for(cancer in cancer.names) 
    {
      stats <- AdjustedResults[[cancer]][[group]][[marker]]
      if(dim(stats)[1] > 0)
      {
        write.table(t(c(marker, group, cancer)), "RNASeqAdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(cbind(stats[1:min(length(stats[[1]]), 10),8], stats[1:min(length(stats[[1]]), 10),6:7]), 
                                                                  "RNASeqAdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}
