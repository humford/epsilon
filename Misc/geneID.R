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
 
      GeneStatsResults[[cancer]][[group]][[marker]] <- hgOver
      adjustedPvalues <- p.adjust(summary(GeneStatsResults[[cancer]][[group]][[marker]])$Pvalue, method = "BH")
      AdjustedResults[[cancer]][[group]][[marker]] <- summary(GeneStatsResults[[cancer]][[group]][[marker]])[which(adjustedPvalues < 0.05),]
      AdjustedResults[[cancer]][[group]][[marker]]$Pvalues <- adjustedPvalues[which(adjustedPvalues < 0.05)]
    }
  }
}


for(marker in c("BP", "CC", "MF", "KEGG"))
{
  for(group in c("UP", "MIDDLE", "DOWN"))
  {
    for(cancer in names(GeneLists)) 
    {
      stats <- AdjustedResults[[cancer]][[group]][[marker]]
      if(dim(stats)[1] > 0)
      {
        write.table(t(c(marker, group, cancer)), "AdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(cbind(stats[,8], stats[,6:7]), "AdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}


for(marker in c("BP", "CC", "MF", "KEGG"))
{
  for(group in c("UP", "MIDDLE", "DOWN"))
  {
    for(cancer in c("BC OV", "BC Glio", "AML_15434", "AML_6891", "LumA")) 
    {
      stats <- AdjustedResults[[cancer]][[group]][[marker]]
      if(dim(stats)[1] > 0)
      {
        write.table(t(c(marker, group, cancer)), "BCAdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(cbind(stats[1:min(length(stats[[1]]), 10),8], stats[1:min(length(stats[[1]]), 10),6:7]), "BCAdjustedResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}



out <- Mclust(controlDiffTable[["BC Glio"]], G = 1:9, modelNames = "V")
plot(out$BIC)
out <- Mclust(controlDiffTable[["BC Glio"]], G = 3, modelNames = "V")


attach(mtcars)
par(mfrow=c(2,2),oma = c(0, 0, 0, 0))

gene <- names(controlDiffTable[["BC Glio"]][which(out$classification == 1)[3]])
plot(density(exprMatrix[gene,]), main = "Control")
plot(density(BatchCorr[gene,]), main = "Glioblastoma")

gene <- names(controlDiffTable[["BC Glio"]][which(out$classification == 3)[12]])
plot(density(exprMatrix[gene,]), main = "Control")
plot(density(BatchCorr[gene,]), main = "Glioblastoma")




exprTable <- matrix(1:15, ncol = 3,nrow = 5)
colnames(exprTable) <- c("Lower", "Upper", "n")
rownames(exprTable) <- c("Glioblastoma", "Ovarian", "AML_1", "AML_2", "LumA", "Control")

par(mai = c(1,2,1,1))
barplot(rbind(exprTable[,"Lower"], exprTable[,"Upper"]), main="Percentage of Genes in Upper and Lower Modes",
        legend = colnames(exprTable)[1:2], beside = FALSE, names = paste(rownames(exprTable), " (n = ", exprTable[,"n"], ")", sep = ""), las = 1, horiz = FALSE)
      

      
      
      
