moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

cancer.names <- c("GBM", "OV", "LumA", "AML", "AML NK")


controlDiffTable <- list()  

#Analyze Control

setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("Control_Processed"))
setwd("~/Documents/Graphs")



ControlSkews <- NULL

for(i in 1:dim(controlMatrix)[[1]])
{
  ControlSkews[i] <-moment((controlMatrix[i,]), 3)
}

names(ControlSkews) <- rownames(controlMatrix)


  
  

#Analyze Cancers
  
for(cancer in cancer.names)
{
  setwd(paste("~/Documents/", cancer, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancer, "ProcessedBC", sep = "_")))
  setwd("~/Documents/Graphs")
  
  skews <- NULL
  for(i in 1:dim(exprMatrix)[[1]])
  {
    skews[i] <-moment((exprMatrix[i,]), 3)
  }

  names(skews) <- rownames(exprMatrix)
  skews <- skews[!is.na(skews)]

  controlDiffTable[[cancer]] <- skews[intersect(names(skews), names(ControlSkews))] - ControlSkews[intersect(names(skews), names(ControlSkews))]

}



library(GOstats)
library(KEGG.db)
library(org.Hs.eg.db)
library(mclust)

GeneLists <- NULL

for(cancer in names(controlDiffTable))
{
  GeneLists[[cancer]] <- list() 
  
  out <- Mclust(controlDiffTable[[cancer]], G = 3, modelNames = "V")
    
  GeneLists[[cancer]][["UP"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 3)] 
  GeneLists[[cancer]][["MIDDLE"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 2)] 
  GeneLists[[cancer]][["DOWN"]] <- names(controlDiffTable[[cancer]])[which(out$classification == 1)]  
  GeneLists[[cancer]][["ALL"]] <- names(controlDiffTable[[cancer]])
}

GeneLists[["Intersect"]] <- GeneLists[["OV"]]

for(cancer in names(controlDiffTable))
{
  GeneLists[["Intersect"]][["UP"]] <- intersect(GeneLists[["Intersect"]][["UP"]], GeneLists[[cancer]][["UP"]])
  GeneLists[["Intersect"]][["MIDDLE"]] <- intersect(GeneLists[["Intersect"]][["MIDDLE"]], GeneLists[[cancer]][["MIDDLE"]]) 
  GeneLists[["Intersect"]][["DOWN"]] <- intersect(GeneLists[["Intersect"]][["DOWN"]], GeneLists[[cancer]][["DOWN"]]) 
  GeneLists[["Intersect"]][["ALL"]] <- intersect(GeneLists[["Intersect"]][["ALL"]], GeneLists[[cancer]][["ALL"]])
}

hgCutoff <- 1
GeneStatsResults <- NULL
AdjustedResults <- NULL

for(cancer in names(GeneLists))
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
        write.table(t(c(marker, group, cancer)), "MicroarrayAdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(stats[1:min(length(stats[[1]]), 10), c(1,8,5,6,7)], 
                                                                  "MicroarrayAdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
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
        write.table(t(c(marker, group, cancer)), "MicroarrayAdjustedResultsNoPval.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(stats[1:min(length(stats[[1]]), 10), c(1,6,7)], 
                                                                  "MicroarrayAdjustedResultsNoPval.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}

      


      

      

