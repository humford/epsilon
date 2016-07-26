calc_samples <- function(expr, name)
{
  lower <- NULL

  for(a in 1:dim(expr)[[2]])
  {
    skews <- NULL
    for(i in 1:dim(expr)[[1]])
    {
      skews[i] <- moment(sample(expr[i,], a), 3)
    }
    lower[a] <- length(which(skews < 0))
  }
  
  
  lower[1:3] <- lower[4]
  lower <- lower/length(skews)
  
  range <- 0.05 *abs(max(lower) - min(lower))
  accept <- which(abs(lower  - lower[length(lower)])  <= range)
  reject <- which(abs(lower  - lower[length(lower)]) > range)


  plot(reject, lower[reject],col = "red", main= paste(name, " AML Gene Spliting by Sample Size"), xlab = "Samples", ylab = "Fraction of Genes in Lower Mode", 
                                          xlim = c(1, length(lower)), ylim = c(min(lower) - range, max(lower) + range))
  points(accept, lower[accept], col = "black")

  abline(a = (lower[length(lower)] - range), b = 0)
  abline(a = (lower[length(lower)] + range), b = 0)
  text(min(accept), lower[min(accept)] + 2*range, as.character(min(accept)))

  return(lower)
}


setwd("~/Documents/AML")



SamplingList <- list()

exprMatrix <- as.matrix(read.table("AML_Processed"))
SamplingList[["TGCA"]] <- calc_samples(exprMatrix, "TGCA")


exprMatrix <- as.matrix(read.table("GSE6891_Processed"))
SamplingList[["GSE6891"]] <- calc_samples(exprMatrix, "GSE6891")


exprMatrix <- as.matrix(read.table("GSE15434_Processed"))
SamplingList[["GSE15434"]] <- calc_samples(exprMatrix, "GSE15434")


exprMatrix <- as.matrix(read.table("GSE17855_Processed"))
SamplingList[["GSE17855"]] <- calc_samples(exprMatrix, "GSE17855")

