moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

categorize_tail <- function(exprs)
{
  return(exprs > mean(exprs) + sd(exprs))
}

mapper <- NULL

probes <- function(gene, map = mapper)
{
  return(mapper[which(mapper[,"symbol"] == gene), "probe"])
}

TCGABarcode <- function(fileName)
{
  return(paste(as.list(strsplit(strsplit(fileName, "lvl-3.")[[1]][2], "-")[[1]][1:3]), sep = "", collapse = "-"))
}

import <- function(){
  #cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
  cancer.names <- c("KIRC")

  for(cancer in cancer.names)
  {
    setwd(paste("~/Documents/", cancer, sep = ""))
    exprMatrix <- as.matrix(read.table(paste(cancer, "Processed", sep = "_")))
  
    skews <- NULL
    for(i in 1:dim(exprMatrix)[[1]])
    {
      skews[i] <-moment((exprMatrix[i,]), 3)
    }

    skews <- skews[!is.na(skews)]
    names(skews) <- rownames(exprMatrix)   
  
    MvalMatrix <- read.table(paste(cancer, "_Methylation_Processed", sep = ""))
    mapper <- read.table("Methylation_Map")
  }
}

