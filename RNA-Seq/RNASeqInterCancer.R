moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")


DiffTable <- list()  

#Analyze Cancers
  
for(cancerA in cancer.names)
{
  DiffTable[[cancerA]] <- list()
  
  setwd(paste("~/Documents/", cancerA, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancerA, "Processed", sep = "_")))
  setwd("~/Documents/Graphs")
    
  skewsA <- NULL
  for(i in 1:dim(exprMatrix)[[1]])
  {
    skewsA[i] <-moment((exprMatrix[i,]), 3)
  }

  skewsA <- skewsA[!is.na(skewsA)]
  names(skewsA) <- rownames(exprMatrix)
    
  for(cancerB in setdiff(cancer.names, cancerA))
  {
    setwd(paste("~/Documents/", cancerB, sep = ""))
    exprMatrix <- as.matrix(read.table(paste(cancerB, "Processed", sep = "_")))
    setwd("~/Documents/Tables")
    
    skewsB <- NULL
    for(i in 1:dim(exprMatrix)[[1]])
    {
      skewsB[i] <-moment((exprMatrix[i,]), 3)
    }

    skewsB <- skewsB[!is.na(skewsB)]
    names(skewsB) <- rownames(exprMatrix)

    DiffTable[[cancerA]][[cancerB]] <- skewsA[intersect(names(skewsA), names(skewsB))] - skewsB[intersect(names(skewsA), names(skewsB))]
  }
}




library(GOstats)
library(KEGG.db)
library(org.Hs.eg.db)
library(mclust)

GeneLists <- NULL

for(cancerA in cancer.names)
{
  GeneLists[[cancerA]] <- list() 
  
  for(cancerB in setdiff(cancer.names, cancerA))
  {
    GeneLists[[cancerA]][[cancerB]] <- list() 
    
    out <- Mclust(DiffTable[[cancerA]][[cancerB]], G = 3, modelNames = "V")
      
    GeneLists[[cancerA]][[cancerB]][["UP"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 3)] 
    GeneLists[[cancerA]][[cancerB]][["MIDDLE"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 2)] 
    GeneLists[[cancerA]][[cancerB]][["DOWN"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 1)]  
    GeneLists[[cancerA]][[cancerB]][["ALL"]] <- names(DiffTable[[cancerA]][[cancerB]])
  }
}

hgCutoff <- 1
pvalueCutoff <- 0.01
GeneStatsResults <- NULL
AdjustedResults <- NULL

for(cancerA in cancer.names)
{
  GeneStatsResults[[cancerA]] <- list()
  AdjustedResults[[cancerA]] <- list()
  
  for(cancerB in setdiff(cancer.names, cancerA))
  {
    GeneStatsResults[[cancerA]][[cancerB]] <- list()
    AdjustedResults[[cancerA]][[cancerB]] <- list()
    
    universelist <- GeneLists[[cancerA]][[cancerB]][["ALL"]]
    
    universe <- unlist(mget(intersect(universelist, ls(org.Hs.egSYMBOL2EG)), org.Hs.egSYMBOL2EG))

    universe <- unique(universe[!is.na(universe)])
    
    for(group in c("UP", "MIDDLE", "DOWN"))
    {
      
      GeneStatsResults[[cancerA]][[cancerB]][[group]] <- list()
      AdjustedResults[[cancerA]][[cancerB]][[group]] <- list()

      genelist <- GeneLists[[cancerA]][[cancerB]][[group]]
        
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
  
        GeneStatsResults[[cancerA]][[cancerB]][[group]][[marker]] <- hgOver
        adjustedPvalues <- p.adjust(summary(GeneStatsResults[[cancerA]][[cancerB]][[group]][[marker]])$Pvalue, method = "BH")
        AdjustedResults[[cancerA]][[cancerB]][[group]][[marker]] <- summary(GeneStatsResults[[cancerA]][[cancerB]][[group]][[marker]])[which(adjustedPvalues < pvalueCutoff),]
        AdjustedResults[[cancerA]][[cancerB]][[group]][[marker]]$Pvalues <- adjustedPvalues[which(adjustedPvalues < pvalueCutoff)]
      }
    }
  }
}

for(marker in c("BP", "CC", "MF", "KEGG"))
{
  for(group in c("UP", "MIDDLE", "DOWN"))
  {
    for(cancerA in cancer.names) 
    {
      for(cancerB in setdiff(cancer.names, cancerA))
      {
        stats <- AdjustedResults[[cancerA]][[cancerB]][[group]][[marker]]
        if(dim(stats)[1] > 0)
        {
          setwd("~/Documents/Tables")
          write.table(t(c(marker, group, cancerA, cancerB)), "RNASeq_InterCancer_TMP.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
          write.table(stats[1:min(length(stats[[1]]), 10), c(1,8,5,6,7)], "RNASeq_InterCancer_TMP.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
          
          write.table(t(c(marker, group, cancerA, cancerB)), "RNASeq_InterCancer_TMP (NO PVAL).txt", append = TRUE, row.names = FALSE, col.names = FALSE)
          write.table(stats[1:min(length(stats[[1]]), 10), c(1,6,7)], "RNASeq_InterCancer_TMP (NO PVAL).txt", append = TRUE, row.names = FALSE, col.names = FALSE)
          
          if(marker != "KEGG")
          {
            setwd("~/Documents/Tables/GOSlim")
            write.table(stats[, c(1,8)], paste(group, cancerA, "to", cancerB, "TPM.txt", sep = "_"), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
          }
        }
      }
    }
  }
}

      

