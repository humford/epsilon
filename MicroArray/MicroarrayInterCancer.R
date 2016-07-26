moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

cancer.names <- c("GBM", "OV", "LumA", "AML", "AML NK")


      
#CANCER TO CANCER COMPARISON 


DiffTable <- list()  
GeneLists <- list()

for(cancerA in cancer.names)
{
  DiffTable[[cancerA]] <- list() 
  GeneLists[[cancerA]] <- list()
  for(cancerB in setdiff(cancer.names, cancerA))
  {
    DiffTable[[cancerA]][[cancerB]] <-list()
    GeneLists[[cancerA]][[cancerB]] <- list()
  }
}

#Analyze Cancers

cancers.used <- NULL
  
for(cancerA in cancer.names)
{
  cancers.used <- c(cancers.used, cancerA)
  
  setwd(paste("~/Documents/", cancerA, sep = ""))
  exprMatrix <- as.matrix(read.table(paste(cancerA, "ProcessedBC", sep = "_")))
  setwd("~/Documents/Graphs")
    
  skewsA <- NULL
  for(i in 1:dim(exprMatrix)[[1]])
  {
    skewsA[i] <-moment((exprMatrix[i,]), 3)
  } 
  names(skewsA) <- rownames(exprMatrix)
  skewsA <- skewsA[!is.na(skewsA)]
    
  for(cancerB in setdiff(cancer.names, cancers.used))
  {
    setwd(paste("~/Documents/", cancerB, sep = ""))
    exprMatrix <- as.matrix(read.table(paste(cancerB, "ProcessedBC", sep = "_")))
    setwd("~/Documents/Graphs")
    
    skewsB <- NULL
    for(i in 1:dim(exprMatrix)[[1]])
    {
      skewsB[i] <-moment((exprMatrix[i,]), 3)
    }
    
    names(skewsB) <- rownames(exprMatrix)
    skewsB <- skewsB[!is.na(skewsB)]

    DiffTable[[cancerA]][[cancerB]] <- skewsA[intersect(names(skewsA), names(skewsB))] - skewsB[intersect(names(skewsA), names(skewsB))]
    DiffTable[[cancerB]][[cancerA]] <- - DiffTable[[cancerA]][[cancerB]]
  }
}




library(GOstats)
library(KEGG.db)
library(org.Hs.eg.db)
library(mclust)


cancers.used <- NULL

for(cancerA in cancer.names)
{
  cancers.used <- c(cancers.used, cancerA)
  
  for(cancerB in setdiff(cancer.names, cancers.used))
  {
    
    out <- Mclust(DiffTable[[cancerA]][[cancerB]], G = 3, modelNames = "V")
      
    GeneLists[[cancerA]][[cancerB]][["UP"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 3)] 
    GeneLists[[cancerA]][[cancerB]][["MIDDLE"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 2)] 
    GeneLists[[cancerA]][[cancerB]][["DOWN"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 1)]  
    GeneLists[[cancerA]][[cancerB]][["ALL"]] <- names(DiffTable[[cancerA]][[cancerB]])
    
    GeneLists[[cancerB]][[cancerA]][["UP"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 1)] 
    GeneLists[[cancerB]][[cancerA]][["MIDDLE"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 2)] 
    GeneLists[[cancerB]][[cancerA]][["DOWN"]] <- names(DiffTable[[cancerA]][[cancerB]])[which(out$classification == 3)]  
    GeneLists[[cancerB]][[cancerA]][["ALL"]] <- names(DiffTable[[cancerA]][[cancerB]])
  }
}

hgCutoff <- 1
GeneStatsResults <- list()
AdjustedResults <- list()

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
        AdjustedResults[[cancerA]][[cancerB]][[group]][[marker]] <- summary(GeneStatsResults[[cancerA]][[cancerB]][[group]][[marker]])[which(adjustedPvalues < 0.05),]
        AdjustedResults[[cancerA]][[cancerB]][[group]][[marker]]$Pvalues <- adjustedPvalues[which(adjustedPvalues < 0.05)]
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
          write.table(t(c(marker, group, cancerA, cancerB)), "MicroarrayInterCancer.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
          write.table(cbind(stats[1:min(length(stats[[1]]), 10),8], stats[1:min(length(stats[[1]]), 10),6:7]), 
                                                                    "MicroarrayInterCancer.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        }
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
          write.table(t(c(marker, group, cancerA, cancerB)), "MicroarrayInterCancerResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
          write.table(stats[1:min(length(stats[[1]]), 10), c(1,8,5,6,7)], 
                                                                    "MicroarrayInterCancerResults.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        }
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
          write.table(t(c(marker, group, cancerA, cancerB)), "MicroarrayInterCancerResultsNoPval.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
          write.table(stats[1:min(length(stats[[1]]), 10), c(1,6,7)], 
                                                                    "MicroarrayInterCancerResultsNoPval.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
        }
      }
    }
  }
}


      
      
