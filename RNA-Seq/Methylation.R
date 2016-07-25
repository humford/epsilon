moment <- function(x, n)
{
  s <- sum((x-mean(x))^n)/(length(x)-1)
  return(abs(s)^(1/n) * sign(s) / sd(x))
}

categorize_tail <- function(exprs)
{
  return(exprs > mean(exprs) + sd(exprs))
}

#cancer.names <- c("SKCM", "HNSC", "LGG", "LUSC", "KIRC")
cancer.names <- c("KIRC")

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

  
}
