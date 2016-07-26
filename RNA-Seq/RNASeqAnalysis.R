moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")


controlDiffTable <- list()  

#Analyze Control

setwd("~/Documents/Control")
controlMatrix <- as.matrix(read.table("RNASeqControlProcessed.txt"))
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
  exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  setwd("~/Documents/Tables")
  
  skews <- NULL
  for(i in 1:dim(exprMatrix)[[1]])
  {
    skews[i] <-moment((exprMatrix[i,]), 3)
  }

  skews <- skews[!is.na(skews)]
  names(skews) <- rownames(exprMatrix)

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

hgCutoff <- 1
pvalueCutoff <- 0.01
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
      AdjustedResults[[cancer]][[group]][[marker]] <- summary(GeneStatsResults[[cancer]][[group]][[marker]])[which(adjustedPvalues < pvalueCutoff),]
      AdjustedResults[[cancer]][[group]][[marker]]$Pvalues <- adjustedPvalues[which(adjustedPvalues < pvalueCutoff)]
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
        setwd("~/Documents/Tables")
        write.table(t(c(marker, group, cancer)), "RNASeq_to_Control_TMP.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(stats[1:min(length(stats[[1]]), 10), c(1,8,5,6,7)], "RNASeq_to_Control_TMP.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        
        write.table(t(c(marker, group, cancer)), "RNASeq_to_Control_TMP (NO PVAL).txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        write.table(stats[1:min(length(stats[[1]]), 10), c(1,6,7)], "RNASeq_to_Control_TMP (NO PVAL).txt", append = TRUE, row.names = FALSE, col.names = FALSE) 
        
        if(marker != "KEGG")
        {
          setwd("~/Documents/Tables/GOSlim")
          write.table(stats[, c(1,8)], paste(group, cancer, "to_Control_TPM.txt", sep = "_"), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
      }
    }
  }
}


      

